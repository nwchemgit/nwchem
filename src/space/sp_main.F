      subroutine sp_setup(node,iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa)
c
c $Id: sp_main.F,v 1.17 1999-08-26 15:26:59 d3j191 Exp $
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
c
      integer i_qlst,l_qlst
c
      call sp_gagetixv(node,iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa,
     + int_mb(i_ipl))
c
      if(me.eq.node) then
      nwmloc=numwm
      nsaloc=numsa
      call sp_qatoms(isl,numsa)
      if(.not.ma_push_get(mt_int,nsa,'qlst',l_qlst,i_qlst))
     + call errquit('Failed to allocate memory for qlst',0)
      call sp_latoms(isl,int_mb(i_qlst),numsa)
      if(.not.ma_pop_stack(l_qlst))
     + call errquit('Failed to deallocate memory for qlst',0)
      call sp_putix(me,iwl,xw,numwm,isl,xs,numsa)
      endif
c
      return
      end
      subroutine sp_latoms(isl,lst,nums)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      integer isl(msa,mis2),lst(nsa)
c
      integer i
c
      do 1 i=1,nsa
      lst(i)=0
    1 continue
c
      call cf_links(nums,isl(1,lsgan),isl(1,lsdyn),lst)
c
      call ga_igop(msp_21,lst,nsa,'+')
c
      do 2 i=1,nums
      if(iand(isl(i,lsdyn),lquant).eq.0.and.lst(isl(i,lsgan)).gt.0)
     + isl(i,lsdyn)=ior(isl(i,lsdyn),lqlink)
    2 continue
c
      return
      end
      subroutine sp_qatoms(isl,numsa)
c
      implicit none
c
#include "sp_common.fh"
c
      logical cf_quantuma
      external cf_quantuma
c
      integer numsa
      integer isl(msa,mis2)
c
      integer i
c
      do 1 i=1,numsa
      if(cf_quantuma(isl(i,lsatt)))
     + isl(i,lsdyn)=ior(isl(i,lsdyn),lquant)
    1 continue
c
      return
      end
      subroutine sp_getixv(node,iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
c
      call sp_gagetixv(node,iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa,
     + int_mb(i_ipl))
c
      if(me.eq.node) then
      nwmloc=numwm
      nsaloc=numsa
      endif
c
      return
      end
      subroutine sp_putix(node,iwl,xw,numwm,isl,xs,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
c
      call sp_gaputix(node,iwl,xw,numwm,isl,xs,numsa,int_mb(i_ipl))
c
      return
      end
      subroutine sp_update(node,vlatt,
     + iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3),vlatt(3,3)
c
      integer i,j
c
      call sp_gaputixv(node,iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa,
     + int_mb(i_ipl))
c
      do 1 j=1,3
      box(j)=vlat(j,j)
      boxh(j)=half*box(j)
      do 2 i=1,3
      vlat(i,j)=vlatt(i,j)
      vlati(i,j)=vlatt(i,j)
    2 continue
    1 continue
c
      call matinv(vlati,3,3)
c
      if(me.eq.node) then
      nwmloc=numwm
      nsaloc=numsa
      endif
c
      return
      end
      subroutine sp_gagetixv(node,iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa,
     + ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
      integer ipl(mbox,mip2)
      integer il,ih,jl,jh
c
      call ga_distribution(ga_ip,node,il,ih,jl,jh)
      call ga_get(ga_ip,il,ih,jl,jh,ipl,mbox)
      numwm=ipl(1,2)
      numsa=ipl(2,2)
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,node,il,ih,jl,jh)
      call ga_get(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_get(ga_w,il,ih,jl,jl+3*mwa-1,xw,mwm)
      call ga_get(ga_w,il,ih,jl+3*mwa,jl+6*mwa-1,vw,mwm)
      call ga_get(ga_w,il,ih,jl+6*mwa,jl+6*mwa+2,xwcr,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,node,il,ih,jl,jh)
      call ga_get(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_get(ga_s,il,ih,jl,jl+2,xs,msa)
      call ga_get(ga_s,il,ih,jl+3,jl+5,vs,msa)
      endif
c
      return
      end
      subroutine sp_gaputixv(node,iwl,xw,xwcr,vw,numwm,isl,xs,vs,numsa,
     + ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3)
      integer ipl(mbox,mip2)
      integer il,ih,jl,jh
c
      ipl(1,2)=numwm
      ipl(2,2)=numsa
      call ga_distribution(ga_ip,node,il,ih,jl,jh)
      call ga_put(ga_ip,il,ih,jl,jh,ipl,mbox)
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,node,il,ih,jl,jh)
      call ga_put(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_put(ga_w,il,ih,jl,jl+3*mwa-1,xw,mwm)
      call ga_put(ga_w,il,ih,jl+3*mwa,jl+6*mwa-1,vw,mwm)
      call ga_put(ga_w,il,ih,jl+6*mwa,jl+6*mwa+2,xwcr,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,node,il,ih,jl,jh)
      call ga_put(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_put(ga_s,il,ih,jl,jl+2,xs,msa)
      call ga_put(ga_s,il,ih,jl+3,jl+5,vs,msa)
      endif
c
      return
      end
      subroutine sp_gaputix(node,iwl,xw,numwm,isl,xs,numsa,ipl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer node,numwm,numsa
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer ipl(mbox,mip2)
      integer il,ih,jl,jh
c
      ipl(1,2)=numwm
      ipl(2,2)=numsa
      call ga_distribution(ga_ip,node,il,ih,jl,jh)
      call ga_put(ga_ip,il,ih,jl,jh,ipl,mbox)
c
      if(numwm.gt.0) then
      call ga_distribution(ga_iw,node,il,ih,jl,jh)
      call ga_put(ga_iw,il,il+numwm-1,jl,jh,iwl,mwm)
      call ga_distribution(ga_w,node,il,ih,jl,jh)
      ih=il+numwm-1
      call ga_put(ga_w,il,ih,jl,jl+3*mwa-1,xw,mwm)
      endif
c
      if(numsa.gt.0) then
      call ga_distribution(ga_is,node,il,ih,jl,jh)
      call ga_put(ga_is,il,il+numsa-1,jl,jh,isl,msa)
      call ga_distribution(ga_s,node,il,ih,jl,jh)
      ih=il+numsa-1
      call ga_put(ga_s,il,ih,jl,jl+2,xs,msa)
      endif
c
      return
      end
      subroutine sp_owner(xw,iwl,numwm,xs,isl,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer numwm,numsa
c
c     determine solvent ownership
c
      call sp_ownerw(xw,iwl,numwm,dbl_mb(i_boxs),int_mb(i_iown))
c
c     determine solute ownership
c
      call sp_owners(xs,isl,numsa,dbl_mb(i_boxs),int_mb(i_iown))
c
      return
      end
      subroutine sp_ownerw(xw,iwl,numwm,boxsiz,ibownr)
c
      implicit none
c
#include "sp_common.fh"
c
      real*8 xw(mwm,3,mwa)
      real*8 boxsiz(maxbox,3)
      integer iwl(mwm,miw2),ibownr(maxbox,3)
      integer numwm
c
      integer iwm,ibx,iby,ibz,i,k,ipx,ipy,ipz
      real*8 cgx,cgy,cgz,xt(3)
c
c     determine the box and owning node for each solvent molecule
c
      do 1 iwm=1,numwm
      cgx=zero
      cgy=zero
      cgz=zero
      do 2 k=1,nwa
      cgx=cgx+xw(iwm,1,k)
      cgy=cgy+xw(iwm,2,k)
      cgz=cgz+xw(iwm,3,k)
    2 continue
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
      do 3 i=1,nbx-1
      if(xt(1)/nwa+boxh(1).gt.boxsiz(i,1)) ibx=i
    3 continue
      do 4 i=1,nby-1
      if(xt(2)/nwa+boxh(2).gt.boxsiz(i,2)) iby=i
    4 continue
      do 5 i=1,nbz-1
      if(xt(3)/nwa+boxh(3).gt.boxsiz(i,3)) ibz=i
    5 continue
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
      iwl(iwm,lwbox)=(ibz*nby+iby)*nbx+ibx
      iwl(iwm,lwnod)=(ipz*npy+ipy)*npx+ipx
    1 continue
c
      return
      end
      subroutine sp_owners(xs,isl,numsa,boxsiz,ibownr)
c
      implicit none
c
#include "sp_common.fh"
c
      real*8 xs(msa,3)
      real*8 boxsiz(maxbox,3)
      integer ibownr(maxbox,3)
      integer isl(msa,mis2)
      integer numsa
c
      integer isa,isaf,isal,jsa,ibox,inod
      integer ibx,iby,ibz,i,ipx,ipy,ipz
      real*8 xcgx,xcgy,xcgz,factor
c
      isaf=1
      isal=0
      do 5 isa=1,numsa
c
      if(isa.lt.numsa) then
      if(isl(isa+1,lssgm).ne.isl(isaf,lssgm)) isal=isa
      else
      isal=isa
      endif
c
      if(isal.gt.0) then
      xcgx=zero
      xcgy=zero
      xcgz=zero
      do 6 jsa=isaf,isal
      xcgx=xcgx+xs(jsa,1)
      xcgy=xcgy+xs(jsa,2)
      xcgz=xcgz+xs(jsa,3)
    6 continue
      factor=one/dble(isal-isaf+1)
      xcgx=factor*xcgx
      xcgy=factor*xcgy
      xcgz=factor*xcgz
c
      ibx=0
      iby=0
      ibz=0
      do 7 i=1,nbx-1
      if(xcgx+boxh(1).gt.boxsiz(i,1)) ibx=i
    7 continue
      do 8 i=1,nby-1
      if(xcgy+boxh(2).gt.boxsiz(i,2)) iby=i
    8 continue
      do 9 i=1,nbz-1
      if(xcgz+boxh(3).gt.boxsiz(i,3)) ibz=i
    9 continue
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
      ibox=(ibz*nby+iby)*nbx+ibx
      inod=(ipz*npy+ipy)*npx+ipx
      do 10 jsa=isaf,isal
      isl(jsa,lsbox)=ibox
      isl(jsa,lsnod)=inod
   10 continue
c
      isaf=isa+1
      isal=0
      endif
c
    5 continue
c
      return
      end
      subroutine sp_travel(bx,xw,vw,xwcr,iwl,numwm,xs,vs,isl,numsa)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),bx(3)
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer numwm,numsa
c
      integer lenx,i_ndx,l_ndx,i_itmp,l_itmp,i_rtmp,l_rtmp
c
      call ga_sync()
c
      box(1)=bx(1)
      box(2)=bx(2)
      box(3)=bx(3)
c
      call sp_bscale(dbl_mb(i_boxs))
c
      call sp_owner(xw,iwl,nwmloc,xs,isl,nsaloc)
c
      lenx=max(nwm,nsa)
      if(.not.ma_push_get(mt_int,lenx,'ndx',l_ndx,i_ndx))
     + call errquit('Failed to allocate ndx',0)
      if(.not.ma_push_get(mt_int,lenx,'itmp',l_itmp,i_itmp))
     + call errquit('Failed to allocate itmp',0)
      if(.not.ma_push_get(mt_dbl,lenx,'rtmp',l_rtmp,i_rtmp))
     + call errquit('Failed to allocate rtmp',0)
c
      call sp_trvl(xw,vw,xwcr,iwl,xs,vs,isl,
     + dbl_mb(i_boxs),int_mb(i_iown),int_mb(i_ipl),
     + int_mb(i_ndx),int_mb(i_itmp),dbl_mb(i_rtmp),lenx)
c
      call ga_sync()
c
      if(.not.ma_pop_stack(l_rtmp))
     + call errquit('Failed to de-allocate rtmp',0)
      if(.not.ma_pop_stack(l_itmp))
     + call errquit('Failed to de-allocate itmp',0)
      if(.not.ma_pop_stack(l_ndx))
     + call errquit('Failed to de-allocate ndx',0)
c
      call sp_gaputixv(me,iwl,xw,xwcr,vw,nwmloc,isl,xs,vs,nsaloc,
     + int_mb(i_ipl))
c
      numwm=nwmloc
      numsa=nsaloc
c
      return
      end
      subroutine sp_trvl(xw,vw,xwcr,iwl,xs,vs,isl,
     + boxsiz,ibownr,ipl,ndx,itmp,rtmp,lenx)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3)
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer lenx
      real*8 boxsiz(maxbox,3)
      integer ipl(mbox,mip2),ibownr(maxbox,3)
      integer ndx(lenx),itmp(lenx)
      real*8 rtmp(lenx)
      logical lrec(27)
c
      integer i,indexw,indexs,j,k,ibx,iby,ibz,ipx,ipy,ipz
      integer isbox,isnod,nrbox
      integer ilp,ihp,jlp,jhp
      integer il,ih,jl,jh,ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      integer iliw,ihiw,jliw,jhiw
      integer ilis,ihis,jlis,jhis
      integer iwm,iwstay,jwstay,lwstay,nwgo,nwgosm
      integer nwgtsm
      integer nwnew,nwstay,iwmloc,jwmloc,lwmloc,irw
      integer isa,jsa,isstay,jsstay,icsgm,ifsgm,ilsgm
      integer nsnew,nsstay,lsstay,isaloc,jsaloc,lsaloc,irs
      integer nsgo,jnode,iwfr,iwto,isfr,isto
      real*8 factor,xscx,xscy,xscz
      integer itemps
      logical lend
c
      nwstay=0
c
c     order the solvent molecules
c
      if(nwmloc.gt.0) then
      do 1 i=1,nwmloc
      ndx(i)=i
    1 continue
      endif
      if(nwmloc.gt.1) then
      lwmloc=nwmloc/2+1
      irw=nwmloc
    2 continue
      if(lwmloc.gt.1) then
      lwmloc=lwmloc-1
      itemps=ndx(lwmloc)
      else
      itemps=ndx(irw)
      ndx(irw)=ndx(1)
      irw=irw-1
      if(irw.eq.1) then
      ndx(1)=itemps
      goto 3
      endif
      endif
      iwmloc=lwmloc
      jwmloc=lwmloc+lwmloc
    4 continue
      if(jwmloc.le.irw) then
      if(jwmloc.lt.irw) then
      if((iwl(ndx(jwmloc),lwnod).eq.iwl(ndx(jwmloc+1),lwnod).and.
     + iwl(ndx(jwmloc),lwbox).le.iwl(ndx(jwmloc+1),lwbox)).or.
     + ((iwl(ndx(jwmloc),lwnod).eq.me.or.
     + (iwl(ndx(jwmloc),lwnod).ne.me.and.
     + iwl(ndx(jwmloc),lwnod).le.iwl(ndx(jwmloc+1),lwnod))).and.
     + iwl(ndx(jwmloc+1),lwnod).ne.me)) jwmloc=jwmloc+1
      endif
      if((iwl(itemps,lwnod).eq.iwl(ndx(jwmloc),lwnod).and.
     + iwl(itemps,lwbox).le.iwl(ndx(jwmloc),lwbox)).or.
     + ((iwl(itemps,lwnod).eq.me.or. (iwl(itemps,lwnod).ne.me.and.
     + iwl(itemps,lwnod).le.iwl(ndx(jwmloc),lwnod))).and.
     + iwl(ndx(jwmloc),lwnod).ne.me)) then
      ndx(iwmloc)=ndx(jwmloc)
      iwmloc=jwmloc
      jwmloc=jwmloc+jwmloc
      else
      jwmloc=irw+1
      endif
      goto 4
      endif
      ndx(iwmloc)=itemps
      goto 2
    3 continue
c
      do 5 k=1,3
      do 8 i=1,nwmloc
      rtmp(i)=xwcr(i,k)
    8 continue
      do 9 i=1,nwmloc
      xwcr(i,k)=rtmp(ndx(i))
    9 continue
      do 10 j=1,nwa
      do 11 i=1,nwmloc
      rtmp(i)=xw(i,k,j)
   11 continue
      do 12 i=1,nwmloc
      xw(i,k,j)=rtmp(ndx(i))
   12 continue
      do 13 i=1,nwmloc
      rtmp(i)=vw(i,k,j)
   13 continue
      do 14 i=1,nwmloc
      vw(i,k,j)=rtmp(ndx(i))
   14 continue
   10 continue
    5 continue
      do 18 k=1,miw2
      do 19 i=1,nwmloc
      itmp(i)=iwl(i,k)
   19 continue
      do 20 i=1,nwmloc
      iwl(i,k)=itmp(ndx(i))
   20 continue
   18 continue
      endif
c
      if(nwmloc.gt.0) then
      do 21 iwm=1,nwmloc
      if(iwl(iwm,lwnod).eq.me) nwstay=iwm
   21 continue
      endif
c
c     order the solute atoms
c
c     isl(isa,lsbox) : box
c     isl(isa,lsnod) : node
c     isl(isa,lssgm) : segment
c
      nsstay=0
      if(nsaloc.gt.0) then
      do 22 i=1,nsaloc
      ndx(i)=i
   22 continue
      endif
c
      if(nsaloc.gt.1) then
      lsaloc=nsaloc/2+1
      irs=nsaloc
   23 continue
      if(lsaloc.gt.1) then
      lsaloc=lsaloc-1
      itemps=ndx(lsaloc)
      else
      itemps=ndx(irs)
      ndx(irs)=ndx(1)
      irs=irs-1
      if(irs.eq.1) then
      ndx(1)=itemps
      goto 24
      endif
      endif
      isaloc=lsaloc
      jsaloc=lsaloc+lsaloc
   25 continue
      if(jsaloc.le.irs) then
      if(jsaloc.lt.irs) then
      if((isl(ndx(jsaloc),lsnod).eq.isl(ndx(jsaloc+1),lsnod).and.
     + (isl(ndx(jsaloc),lsbox).lt.isl(ndx(jsaloc+1),lsbox).or.
     + (isl(ndx(jsaloc),lsbox).eq.isl(ndx(jsaloc+1),lsbox).and.
     + isl(ndx(jsaloc),lssgm).le.isl(ndx(jsaloc+1),lssgm)))).or.
     + ((isl(ndx(jsaloc),lsnod).eq.me.or.
     + (isl(ndx(jsaloc),lsnod).ne.me.and.
     + isl(ndx(jsaloc),lsnod).le.isl(ndx(jsaloc+1),lsnod))).and.
     + isl(ndx(jsaloc+1),lsnod).ne.me)) jsaloc=jsaloc+1
      endif
      if((isl(itemps,lsnod).eq.isl(ndx(jsaloc),lsnod).and.
     + (isl(itemps,lsbox).lt.isl(ndx(jsaloc),lsbox).or.
     + (isl(itemps,lsbox).eq.isl(ndx(jsaloc),lsbox).and.
     + isl(itemps,lssgm).le.isl(ndx(jsaloc),lssgm)))).or.
     + ((isl(itemps,lsnod).eq.me.or. (isl(itemps,lsnod).ne.me.and.
     + isl(itemps,lsnod).le.isl(ndx(jsaloc),lsnod))).and.
     + isl(ndx(jsaloc),lsnod).ne.me)) then
      ndx(isaloc)=ndx(jsaloc)
      isaloc=jsaloc
      jsaloc=jsaloc+jsaloc
      else
      jsaloc=irs+1
      endif
      goto 25
      endif
      ndx(isaloc)=itemps
      goto 23
   24 continue
c
      do 26 k=1,3
      do 27 i=1,nsaloc
      rtmp(i)=xs(i,k)
   27 continue
      do 28 i=1,nsaloc
      xs(i,k)=rtmp(ndx(i))
   28 continue
      do 29 i=1,nsaloc
      rtmp(i)=vs(i,k)
   29 continue
      do 30 i=1,nsaloc
      vs(i,k)=rtmp(ndx(i))
   30 continue
   26 continue
      do 40 k=1,mis2
      do 41 i=1,nsaloc
      itmp(i)=isl(i,k)
   41 continue
      do 42 i=1,nsaloc
      isl(i,k)=itmp(ndx(i))
   42 continue
   40 continue
      endif
c
      if(nsa.gt.0) then
      do 43 isa=1,nsaloc
      if(isl(isa,lsnod).eq.me) nsstay=isa
   43 continue
      endif
c
c     make packages ready for shipment
c
c     loop over all neighboring nodes
c
      call ga_distribution(ga_iw,me,iliw,ihiw,jliw,jhiw)
      call ga_distribution(ga_w,me,ilw,ihw,jlw,jhw)
      call ga_distribution(ga_is,me,ilis,ihis,jlis,jhis)
      call ga_distribution(ga_s,me,ils,ihs,jls,jhs)
c
      indexw=0
      indexs=0
      nwgosm=0
c
      do 70 i=1,27
      jnode=neighb(i,1)
      if(jnode.ge.0.and.jnode.ne.me) then
c
c     for the solvent
c
      iwfr=0
      iwto=0
      do 71 iwm=nwstay+1,nwmloc
      if(iwl(iwm,lwnod).eq.jnode) then
      if(iwfr.eq.0) iwfr=iwm
      iwto=iwm
      endif
   71 continue
c
c     if molecules need to travel copy coordinates etc into global array
c
      nwgo=iwto-iwfr+1
      if(iwfr.eq.0) nwgo=0
      ipl(1,1)=0
      ipl(1,2)=0
c
      if(nwgo.gt.0) then
      nwgosm=nwgosm+nwgo
      il=iliw+indexw
      ih=il+nwgo-1
      call ga_put(ga_iw,il,ih,jliw,jhiw,iwl(iwfr,1),mwm)
      il=ilw+indexw
      ih=il+nwgo-1
      call ga_put(ga_w,il,ih,jlw,jlw+3*mwa-1,xw(iwfr,1,1),mwm)
      call ga_put(ga_w,il,ih,jlw+3*mwa,jlw+6*mwa-1,vw(iwfr,1,1),mwm)
      call ga_put(ga_w,il,ih,jlw+6*mwa,jlw+6*mwa+2,xwcr(iwfr,1),mwm)
      ipl(1,1)=indexw+1
      ipl(1,2)=indexw+nwgo
      indexw=indexw+nwgo
      endif
c
c     for the solute
c
      isfr=0
      isto=0
      do 72 isa=nsstay+1,nsaloc
      if(isl(isa,lsnod).eq.jnode) then
      if(isfr.eq.0) isfr=isa
      isto=isa
      endif
   72 continue
      nsgo=isto-isfr+1
      if(isfr.eq.0) nsgo=0
      ipl(1,3)=0
      ipl(1,4)=0
      if(nsgo.gt.0) then
      il=ilis+indexs
      ih=il+nsgo-1
      call ga_put(ga_is,il,ih,jlis,jhis,isl(isfr,1),msa)
      call ga_put(ga_s,il,ih,jls,jls+2,xs(isfr,1),msa)
      call ga_put(ga_s,il,ih,jls+3,jls+5,vs(isfr,1),msa)
      ipl(1,3)=indexs+1
      ipl(1,4)=indexs+nsgo
      indexs=indexs+nsgo
      endif
c
c     inform other node of number of molecules to get
c
      if(ipl(1,1).gt.0.or.ipl(1,3).gt.0) then
      call ga_distribution(ga_ip,jnode,ilp,ihp,jlp,jhp)
      ilp=ilp+2+i
      call ga_put(ga_ip,ilp,ilp,jlp,jhp,ipl,mbox)
      endif
      endif
      lrec(i)=.false.
   70 continue
c
      call ga_sync()
c
c     receive molecules from other nodes
c
      nwgtsm=0
c
      call ga_distribution(ga_ip,me,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ilp+30,jlp,jhp,ipl,mbox)
c
      do 74 i=1,27
      jnode=neighb(i,2)
      if(jnode.ge.0.and.jnode.ne.me.and..not.lrec(i)) then
c
      iwfr=ipl(3+i,1)
      iwto=ipl(3+i,2)
      isfr=ipl(3+i,3)
      isto=ipl(3+i,4)
c
      nwnew=iwto-iwfr+1
      nsnew=isto-isfr+1
c
      if(iwfr.eq.0) nwnew=0
      if(isfr.eq.0) nsnew=0
c
      if(nwstay+nwnew.gt.mwm)
     + call errquit('Travel: Increase mwm',nwstay+nwnew)
      if(nsstay+nsnew.gt.msa)
     + call errquit('Travel: Increase msa',nsstay+nsnew)
c
      lrec(i)=.true.
c
      if(iwfr.gt.0) then
      nwgtsm=nwgtsm+nwnew
      iwto=ipl(3+i,2)
      call ga_distribution(ga_iw,jnode,iliw,ihiw,jliw,jhiw)
      call ga_distribution(ga_w,jnode,ilw,ihw,jlw,jhw)
c
c     get data for additional molecules
c
      il=iliw+iwfr-1
      ih=iliw+iwto-1
      call ga_get(ga_iw,il,ih,jliw,jhiw,iwl(nwstay+1,1),mwm)
      call ga_get(ga_w,il,ih,jlw,jlw+3*mwa-1,xw(nwstay+1,1,1),mwm)
      call ga_get(ga_w,il,ih,jlw+3*mwa,jlw+6*mwa-1,vw(nwstay+1,1,1),mwm)
      call ga_get(ga_w,il,ih,jlw+6*mwa,jlw+6*mwa+2,xwcr(nwstay+1,1),mwm)
c
c     update number of local solvent molecules
c
      nwstay=nwstay+nwnew
c
      endif
c
c     for the solute
c
      if(isfr.gt.0) then
      call ga_distribution(ga_is,jnode,ilis,ihis,jlis,jhis)
      call ga_distribution(ga_s,jnode,ils,ihs,jls,jhs)
      il=ilis+isfr-1
      ih=ilis+isto-1
      jl=jlis
      jh=jhis
      call ga_get(ga_is,il,ih,jlis,jhis,isl(nsstay+1,1),msa)
      call ga_get(ga_s,il,ih,jls,jls+2,xs(nsstay+1,1),msa)
      call ga_get(ga_s,il,ih,jls+3,jls+5,vs(nsstay+1,1),msa)
c
      nsstay=nsstay+nsnew
      endif
c
      endif
c
c     reset the pointers to zero
c
      ipl(3+i,1)=0
      ipl(3+i,2)=0
      ipl(3+i,3)=0
      ipl(3+i,4)=0
c
   74 continue
c
c     reset ipl in global array
c
      call ga_put(ga_ip,ilp,ilp+30,jlp,jhp,ipl,mbox)
c
c     order the solvent molecules according to subbox and
c     store indices into ip
c
c     ip(1,1)    : number of boxes on this node
c     ip(1,2)    : number of solvent molecules on this node
c     ip(2,2)    : number of solute atoms on this node
c
c     ip(3+i,1)  : index for solvents to be moved to the i-th neighbor
c
c     ip(30+i,1) : number of i-th box on this node
c     ip(30+i,2) : index to first solvent in i-th box
c     ip(30+i,3) : index to lasst solvent in i-th box
c
      if(nwstay.gt.0.and.(nwgosm.gt.0.or.nwgtsm.gt.0)) then
      do 81 i=1,nwstay
      ndx(i)=i
   81 continue
      if(nwstay.gt.1) then
      lwstay=nwstay/2+1
      irw=nwstay
   82 continue
      if(lwstay.gt.1) then
      lwstay=lwstay-1
      itemps=ndx(lwstay)
      else
      itemps=ndx(irw)
      ndx(irw)=ndx(1)
      irw=irw-1
      if(irw.eq.1) then
      ndx(1)=itemps
      goto 83
      endif
      endif
      iwstay=lwstay
      jwstay=lwstay+lwstay
   84 continue
      if(jwstay.le.irw) then
      if(jwstay.lt.irw) then
      if(iwl(ndx(jwstay),lwbox).le.iwl(ndx(jwstay+1),lwbox))
     + jwstay=jwstay+1
      endif
      if(iwl(itemps,lwbox).le.iwl(ndx(jwstay),lwbox)) then
      ndx(iwstay)=ndx(jwstay)
      iwstay=jwstay
      jwstay=jwstay+jwstay
      else
      jwstay=irw+1
      endif
      goto 84
      endif
      ndx(iwstay)=itemps
      goto 82
   83 continue
c
      do 85 k=1,3
      do 88 i=1,nwstay
      rtmp(i)=xwcr(i,k)
   88 continue
      do 89 i=1,nwstay
      xwcr(i,k)=rtmp(ndx(i))
   89 continue
      do 90 j=1,nwa
      do 91 i=1,nwstay
      rtmp(i)=xw(i,k,j)
   91 continue
      do 92 i=1,nwstay
      xw(i,k,j)=rtmp(ndx(i))
   92 continue
      do 93 i=1,nwstay
      rtmp(i)=vw(i,k,j)
   93 continue
      do 94 i=1,nwstay
      vw(i,k,j)=rtmp(ndx(i))
   94 continue
   90 continue
   85 continue
      do 98 k=1,miw2
      do 99 i=1,nwstay
      itmp(i)=iwl(i,k)
   99 continue
      do 100 i=1,nwstay
      iwl(i,k)=itmp(ndx(i))
  100 continue
   98 continue
c
      endif
      endif
c
c     order the solute according to segment
c
      if(nsstay.gt.0) then
      do 122 i=1,nsstay
      ndx(i)=i
  122 continue
      if(nsstay.gt.1) then
      lsstay=nsstay/2+1
      irs=nsstay
  123 continue
      if(lsstay.gt.1) then
      lsstay=lsstay-1
      itemps=ndx(lsstay)
      else
      itemps=ndx(irs)
      ndx(irs)=ndx(1)
      irs=irs-1
      if(irs.eq.1) then
      ndx(1)=itemps
      goto 124
      endif
      endif
      isstay=lsstay
      jsstay=lsstay+lsstay
  125 continue
      if(jsstay.le.irs) then
      if(jsstay.lt.irs) then
      if(isl(ndx(jsstay),lssgm).le.isl(ndx(jsstay+1),lssgm))
     + jsstay=jsstay+1
      endif
      if(isl(itemps,lssgm).le.isl(ndx(jsstay),lssgm)) then
      ndx(isstay)=ndx(jsstay)
      isstay=jsstay
      jsstay=jsstay+jsstay
      else
      jsstay=irs+1
      endif
      goto 125
      endif
      ndx(isstay)=itemps
      goto 123
  124 continue
      endif
c
c     for each segment : 1. determine box number
c                        2. assign box number to each atom
c                        3. when box not owned by node:
c                           a. assign box number
c                           b. assign correct node number
c
      icsgm=isl(ndx(1),lssgm)
      ifsgm=1
      ilsgm=1
      do 126 isa=2,nsstay+1
c
c     if isa is first atom of a new segment or very last atom
c
      
      if(isa.le.nsstay) then
      lend=isl(ndx(isa),lssgm).ne.icsgm
      else
      lend=.true.
      endif
      if(lend) then
      if(isa.gt.nsstay) ilsgm=nsstay
      if(ifsgm.gt.0.and.ilsgm.ge.ifsgm) then
      xscx=zero
      xscy=zero
      xscz=zero
      do 127 jsa=ifsgm,ilsgm
      xscx=xscx+xs(ndx(jsa),1)
      xscy=xscy+xs(ndx(jsa),2)
      xscz=xscz+xs(ndx(jsa),3)
  127 continue
      factor=one/dble(ilsgm-ifsgm+1)
      xscx=factor*xscx
      xscy=factor*xscy
      xscz=factor*xscz
c
c     determine the box number
c
      ibx=0
      iby=0
      ibz=0
      do 128 i=1,nbx-1
      if(xscx+boxh(1).gt.boxsiz(i,1)) ibx=i
  128 continue
      do 129 i=1,nby-1
      if(xscy+boxh(2).gt.boxsiz(i,2)) iby=i
  129 continue
      do 1130 i=1,nbz-1
      if(xscz+boxh(3).gt.boxsiz(i,3)) ibz=i
 1130 continue
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
      isbox=(ibz*nby+iby)*nbx+ibx
      isnod=(ipz*npy+ipy)*npx+ipx
c
c     assign box and node numbers
c
      do 1131 jsa=ifsgm,ilsgm
      isl(ndx(jsa),lsbox)=isbox
      isl(ndx(jsa),lsnod)=isnod
 1131 continue
c
      endif
      if(isa.le.nsstay) icsgm=isl(ndx(isa),lssgm)
      ifsgm=isa
      else
      ilsgm=isa
      endif
  126 continue
c
c     order solute according to box, segment, charge group, atom number
c
      if(nsstay.gt.1) then
      lsstay=nsstay/2+1
      irs=nsstay
  132 continue
      if(lsstay.gt.1) then
      lsstay=lsstay-1
      itemps=ndx(lsstay)
      else
      itemps=ndx(irs)
      ndx(irs)=ndx(1)
      irs=irs-1
      if(irs.eq.1) then
      ndx(1)=itemps
      goto 133
      endif
      endif
      isstay=lsstay
      jsstay=lsstay+lsstay
  134 continue
      if(jsstay.le.irs) then
      if(jsstay.lt.irs) then
      if(isl(ndx(jsstay),lsbox).lt.isl(ndx(jsstay+1),lsbox).or.
     + (isl(ndx(jsstay),lsbox).eq.isl(ndx(jsstay+1),lsbox).and.
     + (isl(ndx(jsstay),lssgm).lt.isl(ndx(jsstay+1),lssgm).or.
     + (isl(ndx(jsstay),lssgm).eq.isl(ndx(jsstay+1),lssgm).and.
     + (isl(ndx(jsstay),lsgrp).lt.isl(ndx(jsstay+1),lsgrp).or.
     + (isl(ndx(jsstay),lsgrp).eq.isl(ndx(jsstay+1),lsgrp).and.
     + isl(ndx(jsstay),lsgan).le.isl(ndx(jsstay+1),lsgan)))))))
     + jsstay=jsstay+1
      endif
      if(isl(itemps,lsbox).lt.isl(ndx(jsstay),lsbox).or.
     + (isl(itemps,lsbox).eq.isl(ndx(jsstay),lsbox).and.
     + (isl(itemps,lssgm).lt.isl(ndx(jsstay),lssgm).or.
     + (isl(itemps,lssgm).eq.isl(ndx(jsstay),lssgm).and.
     + (isl(itemps,lsgrp).lt.isl(ndx(jsstay),lsgrp).or.
     + (isl(itemps,lsgrp).eq.isl(ndx(jsstay),lsgrp).and.
     + isl(itemps,lsgan).le.isl(ndx(jsstay),lsgan))))))) then
      ndx(isstay)=ndx(jsstay)
      isstay=jsstay
      jsstay=jsstay+jsstay
      else
      jsstay=irs+1
      endif
      goto 134
      endif
      ndx(isstay)=itemps
      goto 132
  133 continue
      endif
c
      do 135 k=1,3
      do 136 i=1,nsstay
      rtmp(i)=xs(i,k)
  136 continue
      do 137 i=1,nsstay
      xs(i,k)=rtmp(ndx(i))
  137 continue
      do 138 i=1,nsstay
      rtmp(i)=vs(i,k)
  138 continue
      do 139 i=1,nsstay
      vs(i,k)=rtmp(ndx(i))
  139 continue
  135 continue
      do 149 k=1,mis2
      do 150 i=1,nsstay
      itmp(i)=isl(i,k)
  150 continue
      do 151 i=1,nsstay
      isl(i,k)=itmp(ndx(i))
  151 continue
  149 continue
c
      endif
c
      do 200 i=1,ipl(1,1)
      ipl(30+i,2)=0
      ipl(30+i,3)=0
      ipl(30+i,4)=0
      ipl(30+i,5)=0
  200 continue
c
      do 201 i=1,ipl(1,1)
      nrbox=ipl(30+i,1)
      if(nwstay.gt.0) then
      do 202 iwm=1,nwstay
      if(iwl(iwm,lwbox).eq.nrbox) then
      if(ipl(30+i,2).eq.0) ipl(30+i,2)=iwm
      ipl(30+i,3)=iwm
      endif
  202 continue
      endif
      if(nsstay.gt.0) then
      do 203 isa=1,nsstay
      if(isl(isa,lsbox).eq.nrbox) then
      if(ipl(30+i,4).eq.0) ipl(30+i,4)=isa
      ipl(30+i,5)=isa
      endif
  203 continue
      endif
  201 continue
c
      nwmloc=nwstay
      ipl(1,2)=nwmloc
      nsaloc=nsstay
      ipl(2,2)=nsaloc
c
      return
      end
      subroutine sp_nbbl(n)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer n
c
      n=nbbl
      if(lqmd) n=0
c
      return
      end
      subroutine sp_gethdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,isfr,isto
      logical lself
c
      call sp_gthdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto,
     + int_mb(i_bb),int_mb(i_ipl))
c
      return
      end
      subroutine sp_gthdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto,
     + lbbl,ipl)
c
      implicit none
c
#include "sp_common.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,isfr,isto
      integer lbbl(mbbl,mbb2)
      integer ipl(mbox,mip2)
      logical lself
c
      integer ndx,i,ibox
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call errquit('Index to lbbl out of range',ibbl)
c
      ibox=lbbl(ibbl,3)
      lhandl=lbbl(ibbl,4)
      lself=ibox.eq.lbbl(ibbl,2)
c
      ndx=0
      do 1 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.ibox) ndx=30+i
    1 continue
      if(ndx.eq.0) call errquit('Local node box not found',0)
      iwfr=ipl(ndx,2)
      iwto=ipl(ndx,3)
      isfr=ipl(ndx,4)
      isto=ipl(ndx,5)
c
      return
      end
      subroutine sp_getxbl(ibbl,lhandl,
     + iwl,xw,iwfr,iwto,jwfr,jwto,isl,xs,isfr,isto,jsfr,jsto,lself)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      logical lself
c
      call sp_gtxbl(ibbl,lhandl,
     + iwl,xw,iwfr,iwto,jwfr,jwto,isl,xs,isfr,isto,jsfr,jsto,
     + int_mb(i_bb),int_mb(i_ipl),int_mb(i_jpl),lself)
c
      return
      end
      subroutine sp_gtxbl(ibbl,lhandl,
     + iwl,xw,iwfr,iwto,jwfr,jwto,isl,xs,isfr,isto,jsfr,jsto,
     + lbbl,ipl,jpl,lself)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer ibbl,lhandl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer lbbl(mbbl,mbb2)
      integer ipl(mbox,mip2),jpl(mbox,mip2)
      logical lself
c
      integer i,ibox,jbox,jnode,il,ih,jl,jh,ndx,nwnon,nsnon
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call errquit('Index to lbbl out of range',ibbl)
c
      jnode=lbbl(ibbl,1)
      jbox=lbbl(ibbl,2)
      ibox=lbbl(ibbl,3)
      lhandl=lbbl(ibbl,4)
c
      lself=ibox.eq.jbox
c
      if(ibbl.eq.1) then
      jcnode=-1
      jcbox=-1
      icbox=-1
      endif
c
      if(ibox.ne.icbox) then
      ndx=0
      do 1 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.ibox) ndx=30+i
    1 continue
      if(ndx.eq.0) call errquit('Local node box not found',0)
      icbox=ibox
      iwfr=ipl(ndx,2)
      iwto=ipl(ndx,3)
      isfr=ipl(ndx,4)
      isto=ipl(ndx,5)
      endif
c
      if(jnode.eq.me) then
c
      if(jbox.ne.jcbox) then
      ndx=0
      do 2 i=1,ipl(1,1)
      if(ipl(30+i,1).eq.jbox) ndx=30+i
    2 continue
      if(ndx.eq.0) call errquit('Local node box not found',0)
      jcbox=jbox
      jwfr=ipl(ndx,2)
      jwto=ipl(ndx,3)
      jsfr=ipl(ndx,4)
      jsto=ipl(ndx,5)
      endif
c
      else
c
      if(jnode.ne.jcnode) then
      call ga_distribution(ga_ip,jnode,il,ih,jl,jh)
      call ga_get(ga_ip,il,ih,jl,jh,jpl,mbox)
      jcnode=jnode
      endif
c
      if(jbox.ne.jcbox) then
      ndx=0
      do 3 i=1,jpl(1,1)
      if(jpl(30+i,1).eq.jbox) ndx=30+i
    3 continue
      if(ndx.eq.0) call errquit('Remote node box not found',0)
      jcbox=jbox
      jwfr=jpl(ndx,2)
      jwto=jpl(ndx,3)
      jsfr=jpl(ndx,4)
      jsto=jpl(ndx,5)
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(nwmloc+nwnon.gt.mwm)
     + call errquit('Dimension mwm too small',nwmloc+nwnon)
      call ga_distribution(ga_iw,jnode,il,ih,jl,jh)
      call ga_get(ga_iw,il+jwfr-1,il+jwto-1,jl,jh,iwl(nwmloc+1,1),mwm)
      call ga_distribution(ga_w,jnode,il,ih,jl,jh)
      call ga_get(ga_w,il+jwfr-1,il+jwto-1,jl,jl+3*mwa-1,
     + xw(nwmloc+1,1,1),mwm)
      jwfr=nwmloc+1
      jwto=nwmloc+1+jpl(ndx,3)-jpl(ndx,2)
      endif
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(nsaloc+nsnon.gt.msa)
     + call errquit('Dimension msa too small',nsaloc+nsnon)
      call ga_distribution(ga_is,jnode,il,ih,jl,jh)
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl,jh,isl(nsaloc+1,1),msa)
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_get(ga_s,il+jsfr-1,il+jsto-1,jl,jl+2,
     + xs(nsaloc+1,1),msa)
      jsfr=nsaloc+1
      jsto=nsaloc+1+jpl(ndx,5)-jpl(ndx,4)
      endif
      endif
c
      endif
c
      return
      end
      subroutine sp_accfbl(ibbl,lhandl,fw,fs,lpair,iwz,isz)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      call sp_acfbl(ibbl,lhandl,fw,fs,int_mb(i_bb),int_mb(i_jpl),
     + lpair,iwz,isz)
c
      return
      end
      subroutine sp_acfbl(ibbl,lhandl,fw,fs,lbbl,jpl,lpair,iwz,isz)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
c
      integer ibbl,lhandl
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      integer lbbl(mbbl,mbb2)
      integer jpl(mbox,mip2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      integer i,j,k,l,m,ibox,jbox,jnode,il,ih,jl,jh,ndx,nwnon,nsnon
      integer jwfr,jwto,jsfr,jsto
c
      if(ibbl.le.0.or.ibbl.gt.nbbl)
     + call errquit('Index to lbbl out of range',ibbl)
c
      jnode=lbbl(ibbl,1)
      jbox=lbbl(ibbl,2)
      ibox=lbbl(ibbl,3)
      lbbl(ibbl,4)=lhandl
c
      if(jnode.eq.me) return
      if(ibbl.lt.nbbl) then
      if(lbbl(ibbl+1,1).eq.jnode.and.lbbl(ibbl+1,2).eq.jbox) return
      endif
c
      ndx=0
      do 3 i=1,jpl(1,1)
      if(jpl(30+i,1).eq.jbox) ndx=30+i
    3 continue
      if(ndx.eq.0) call errquit('Remote node box not found',0)
      jwfr=jpl(ndx,2)
      jwto=jpl(ndx,3)
      jsfr=jpl(ndx,4)
      jsto=jpl(ndx,5)
      if(jwfr.gt.0) then
      nwnon=jwto-jwfr+1
      if(nwmloc+nwnon.gt.mwm)
     + call errquit('Dimension mwm too small',nwmloc+nwnon)
      call ga_distribution(ga_w,jnode,il,ih,jl,jh)
      call ga_acc(ga_w,il+jwfr-1,il+jwto-1,jl+6*mwa+3,jl+9*mwa+2,
     + fw(nwmloc+1,1,1,1),mwm,one)
      if(llong) call ga_acc(ga_w,il+jwfr-1,il+jwto-1,jl+9*mwa+3,
     + jl+12*mwa+2,fw(nwmloc+1,1,1,2),mwm,one)
      if(lpair) then
      call ga_distribution(ga_iwz,jnode,il,ih,jl,jh)
      call ga_acc(ga_iwz,il+jwfr-1,il+jwto-1,jl,jl,iwz(nwmloc+1),mwm,1)
      endif
      endif
      if(jsfr.gt.0) then
      nsnon=jsto-jsfr+1
      if(nsaloc+nsnon.gt.msa)
     + call errquit('Dimension msa too small',nsaloc+nsnon)
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_acc(ga_s,il+jsfr-1,il+jsto-1,jl+6,jl+8,
     + fs(nsaloc+1,1,1),msa,one)
      if(llong) call ga_acc(ga_s,il+jsfr-1,il+jsto-1,jl+9,jl+11,
     + fs(nsaloc+1,1,2),msa,one)
      if(lpair) then
      call ga_distribution(ga_isz,jnode,il,ih,jl,jh)
      call ga_acc(ga_isz,il+jsfr-1,il+jsto-1,jl,jl,isz(nsaloc+1),msa,1)
      endif
      endif
c
      m=1
      if(llong) m=2
c
      do 9 l=1,m
      do 4 k=1,mwa
      do 5 j=1,3
      do 6 i=nwmloc+1,mwm
      fw(i,j,k,l)=zero
    6 continue
    5 continue
    4 continue
c
      do 7 j=1,3
      do 8 i=nsaloc+1,msa
      fs(i,j,l)=zero
    8 continue
    7 continue
    9 continue
c
      return
      end
      subroutine sp_bscale(boxsiz)
c
      implicit none
c
#include "sp_common.fh"
c
      real*8 boxsiz(maxbox,3)
c
      integer ibx,iby,ibz
      real*8 boxscl
c
c     scale the subbox boundaries
c
      boxscl=box(1)/boxsiz(nbx,1)
      do 1 ibx=1,nbx
      boxsiz(ibx,1)=boxscl*boxsiz(ibx,1)
    1 continue
      boxscl=box(2)/boxsiz(nby,2)
      do 2 iby=1,nby
      boxsiz(iby,2)=boxscl*boxsiz(iby,2)
    2 continue
      boxscl=box(3)/boxsiz(nbz,3)
      do 3 ibz=1,nbz
      boxsiz(ibz,3)=boxscl*boxsiz(ibz,3)
    3 continue
c
      return
      end
      subroutine sp_balanc(stime,syntim)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
c
      logical sp_diffbb
      external sp_diffbb
c
      real*8 stime,syntim
c
      integer l_dlb,i_dlb,i_nod,l_nod,i_nen,l_nen,i_lst,l_lst
c
c      if(loadb.eq.0) return
c
      if(.not.ma_verify_allocator_stuff())
     + call errquit('Problems in sp_balanc at 1',me)
      if(.not.ma_push_get(mt_int,np,'nod',l_nod,i_nod))
     + call errquit('Failed to allocate nod',me)
      if(.not.ma_push_get(mt_int,2*np,'lst',l_lst,i_lst))
     + call errquit('Failed to allocate lst',me)
      if(.not.ma_push_get(mt_int,mbbl*np,'nen',l_nen,i_nen))
     + call errquit('Failed to allocate nen',me)
      if(.not.ma_push_get(mt_dbl,np,'dlb',l_dlb,i_dlb))
     + call errquit('Failed to allocate dlb',me)
      if(.not.ma_verify_allocator_stuff())
     + call errquit('Problems in sp_balanc at 2',me)
c
      call sp_dldbal(stime,syntim,int_mb(i_nod),int_mb(i_lst),
     + int_mb(i_nen),int_mb(i_bb),dbl_mb(i_dlb))
      if(.not.ma_verify_allocator_stuff())
     + call errquit('Problems in sp_balanc at 3',me)
c
      if(.not.ma_pop_stack(l_dlb))
     + call errquit('Failed to deallocate dlb',me)
      if(.not.ma_pop_stack(l_nen))
     + call errquit('Failed to deallocate nen',me)
      if(.not.ma_pop_stack(l_lst))
     + call errquit('Failed to deallocate lst',me)
      if(.not.ma_pop_stack(l_nod))
     + call errquit('Failed to deallocate nod',me)
      if(.not.ma_verify_allocator_stuff())
     + call errquit('Problems in sp_balanc at 4',me)
c
      if(sp_diffbb(dbl_mb(i_boxs),int_mb(i_rng))) then
      call sp_numbb(int_mb(i_iown),dbl_mb(i_boxs))
      call sp_listbb(int_mb(i_iown),dbl_mb(i_boxs),int_mb(i_bb))
      endif
      if(.not.ma_verify_allocator_stuff())
     + call errquit('Problems in sp_balanc at 5',me)
c
      return
      end
      subroutine sp_dldbal(stime,syntim,nod,nlst,nen,lbbl,dlb)
c
      implicit none
c
#include "sp_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nod(np),nlst(np,2),nen(np,mbbl)
      integer lbbl(mbbl,mbb2)
      real*8 stime,syntim,synsum,dlb(np),factor,facx,facy,facz
c
      integer i,j,k,n,ibusy
c
c     broad cast times
c
      do 1 i=1,np
      dlb(i)=0.0d0
      nod(i)=i-1
    1 continue
      dlb(me+1)=syntim
      synsum=syntim
c
      call ga_dgop(msp_14,dlb,np,'+')
      call ga_dgop(msp_15,synsum,1,'+')
cx
      if(me.eq.0) then
      write(lfnsyn,3000) stime,synsum
 3000 format('synchronization',/,2f12.6)
      endif
cx
      if(synsum.lt.1.0d-5) return
c
      facx=1.0d0
      facy=1.0d0
      facz=1.0d0
c
c     count load balancing steps with increasing synsum
c
      if(nldup.lt.0) then
      nldup=0
      tsyncp=synsum
      elseif(tsyncp.lt.synsum) then
      nldup=nldup+1
      else
      nldup=0
      tsyncp=synsum
      endif
c
c     find busiest node
c
      ibusy=1
      do 2 i=2,np
      if(dlb(i).lt.dlb(ibusy)) ibusy=i
    2 continue
c
c     order node list
c
      do 3 i=1,np-1
      do 4 j=i+1,np
      if(dlb(nod(i)+1).gt.dlb(nod(j)+1)) then
      n=nod(i)
      nod(i)=nod(j)
      nod(j)=n
      endif
    4 continue
    3 continue
c
c     resize boxes
c     ------------
c
      if(loadb.eq.1.or.(loadb.eq.3.and.nldup.gt.lbpair)) then
      nldup=0
      tsyncp=synsum
      factor=(1.0d0-factld*(synsum/dble(np)-dlb(ibusy))/synsum)**third
      call sp_resize(ibusy-1,int_mb(i_iown),dbl_mb(i_boxs),factor,
     + facx,facy,facz)
      if(me.eq.0) then
      do 16 i=1,np
      if(i.eq.1) then
      write(lfnsyn,3001) nod(i),dlb(nod(i)+1),facx,facy,facz
 3001 format(i5,f12.6,12x,3f12.6)
      else
      write(lfnsyn,3002) nod(i),dlb(nod(i)+1)
 3002 format(i5,f12.6)
      endif
   16 continue
      endif
      return
      endif
c
c     get active neighbors
c
      n=0
      do 5 i=1,np
      do 6 j=1,mbbl
      nen(i,j)=0
    6 continue
      nlst(i,1)=0
      nlst(i,2)=0
    5 continue
      do 7 i=1,nbbl
      if(lbbl(i,1).eq.me) goto 7
      nlst(me+1,2)=nlst(me+1,2)+1
      do 8 j=1,n
      if(lbbl(i,1)+1.eq.nen(me+1,j)) goto 7
    8 continue
      n=n+1
      nen(me+1,n)=lbbl(i,1)+1
    7 continue
      nlst(me+1,1)=nbbl
c
      call ga_igop(msp_16,nen,mbbl*np,'+')
      call ga_igop(msp_17,nlst,2*np,'+')
      if(me.eq.0) then
      do 3031 i=1,np
      write(lfnsyn,3030) i-1,nlst(i,1),nlst(i,2),
     + dlb(i),(nen(i,j)-1,j=1,16)
 3030 format(i3,2i5,f12.6,20i3)
 3031 continue
      endif
c
c     determine the least busy neighbor of each node
c
      do 9 i=1,np
      nlst(i,1)=-1
    9 continue
      do 10 i=1,np-1
      do 11 j=np,i+1,-1
      do 12 k=1,np
      if(nlst(k,1).eq.nod(j)) goto 11
   12 continue
      do 13 k=1,mbbl
      if(nod(j)+1.eq.nen(nod(i)+1,k)) then
      nlst(i,1)=nod(j)
      goto 10
      endif
   13 continue
   11 continue
   10 continue
c
c     setup transfer list
c
      do 14 i=1,np
      nen(nod(i)+1,1)=nlst(i,1)
   14 continue
c
      if(me.eq.0) then
      do 15 i=1,np
      write(lfnsyn,3003) nod(i),dlb(nod(i)+1),nlst(nod(i),2),
     + nen(nod(i)+1,1)
 3003 format(i5,f12.6,i7,i5)
   15 continue
      endif
c
      call ga_brdcst(msp_18,nen,np*ma_sizeof(mt_int,1,mt_byte),0)
c
c     cascade box-box pairs
c
      if(loadb.ge.2) then
      call sp_cascad(nen,int_mb(i_bb))
      endif
c
      return
      end
      subroutine sp_resize(ipmin,ibownr,boxsiz,factor,facx,facy,facz)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ipmin,ibownr(maxbox,3)
      real*8 factor,facx,facy,facz,boxsiz(maxbox,3)
c
      integer ibx,iby,ibz,ipminx,ipminy,ipminz
      real*8 bxi,bxr,byi,byr,bzi,bzr,rfact
c
c     load balancing based on resizing of node domains
c     ------------------------------------------------
c
      ipminx=mod(ipmin,npx)
      ipminy=mod((ipmin-ipminx)/npx,npy)
      ipminz=((ipmin-ipminx)/npx-ipminy)/npy
c
c     resizing boxes in x-dimension
c
      if(npx.gt.1) then
      do 7 ibx=nbx,2,-1
      boxsiz(ibx,1)=boxsiz(ibx,1)-boxsiz(ibx-1,1)
    7 continue
      bxi=zero
      bxr=zero
      do 8 ibx=1,nbx
      if(ibownr(ibx,1).eq.ipminx) then
      facx=1.0d0/boxsiz(ibx,1)
      boxsiz(ibx,1)=factor*boxsiz(ibx,1)
      if(boxsiz(ibx,1).lt.bxmin) boxsiz(ibx,1)=bxmin
      facx=facx*boxsiz(ibx,1)
      bxi=bxi+boxsiz(ibx,1)
      else
      bxr=bxr+boxsiz(ibx,1)
      endif
    8 continue
      rfact=(box(1)-bxi)/bxr
      do 9 ibx=1,nbx
      if(ibownr(ibx,1).ne.ipminx) boxsiz(ibx,1)=rfact*boxsiz(ibx,1)
      if(ibx.gt.1) boxsiz(ibx,1)=boxsiz(ibx,1)+boxsiz(ibx-1,1)
    9 continue
      endif
c
c     resizing boxes in y-dimension
c
      if(npy.gt.1) then
      do 10 iby=nby,2,-1
      boxsiz(iby,2)=boxsiz(iby,2)-boxsiz(iby-1,2)
   10 continue
      byi=zero
      byr=zero
      do 11 iby=1,nby
      if(ibownr(iby,2).eq.ipminy) then
      facy=1.0d0/boxsiz(iby,2)
      boxsiz(iby,2)=factor*boxsiz(iby,2)
      if(boxsiz(iby,2).lt.bymin) boxsiz(iby,2)=bymin
      facy=1.0d0/boxsiz(iby,2)
      byi=byi+boxsiz(iby,2)
      else
      byr=byr+boxsiz(iby,2)
      endif
   11 continue
      rfact=(box(2)-byi)/byr
      do 12 iby=1,nby
      if(ibownr(iby,2).ne.ipminy) boxsiz(iby,2)=rfact*boxsiz(iby,2)
      if(iby.gt.1) boxsiz(iby,2)=boxsiz(iby,2)+boxsiz(iby-1,2)
   12 continue
      endif
c
c     resizing boxes in z-dimension
c
      if(npz.gt.1) then
      do 13 ibz=nbz,2,-1
      boxsiz(ibz,3)=boxsiz(ibz,3)-boxsiz(ibz-1,3)
   13 continue
      bzi=zero
      bzr=zero
      do 14 ibz=1,nbz
      if(ibownr(ibz,3).eq.ipminz) then
      facz=1.0d0/boxsiz(ibz,3)
      boxsiz(ibz,3)=factor*boxsiz(ibz,3)
      if(boxsiz(ibz,3).lt.bzmin) boxsiz(ibz,3)=bzmin
      facz=1.0d0/boxsiz(ibz,3)
      bzi=bzi+boxsiz(ibz,3)
      else
      bzr=bzr+boxsiz(ibz,3)
      endif
   14 continue
      rfact=(box(3)-bzi)/bzr
      do 15 ibz=1,nbz
      if(ibownr(ibz,3).ne.ipminz) boxsiz(ibz,3)=rfact*boxsiz(ibz,3)
      if(ibz.gt.1) boxsiz(ibz,3)=boxsiz(ibz,3)+boxsiz(ibz-1,3)
   15 continue
      endif
c
      call ga_brdcst
     + (msp_19,boxsiz,ma_sizeof(mt_dbl,3*maxbox,mt_byte),0)
c
      return
      end
      subroutine sp_cascad(ltran,lbbl)
c
      implicit none
c
#include "sp_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ltran(np,3),lbbl(mbbl,mbb2)
c
      integer i,j,ltemp,ntargt
c
c     load balancing based on cascade box-list redistribution
c     -------------------------------------------------------
c
      do 1 i=1,np
      ltran(i,2)=0
      ltran(i,3)=0
    1 continue
c
      ntargt=ltran(me+1,1)
      ltran(me+1,2)=-1
      ltran(me+1,3)=-1
c
      if(ntargt.ge.0) then
      do 2 i=1,nbbl
      if(lbbl(i,1).eq.ntargt) then
      ltran(me+1,2)=lbbl(i,2)
      ltran(me+1,3)=lbbl(i,3)
      do 3 j=i,nbbl
      lbbl(j,1)=lbbl(j+1,1)
      lbbl(j,2)=lbbl(j+1,2)
      lbbl(j,3)=lbbl(j+1,3)
    3 continue
      nbbl=nbbl-1
      goto 4
      endif
    2 continue
      endif
c
    4 continue
      call ga_igop(msp_20,ltran(1,2),2*np,'+')
c
      do 5 i=1,np
      if(ltran(i,1).eq.me.and.ltran(i,2).ge.0) then
      nbbl=nbbl+1
      lbbl(nbbl,1)=i-1
      lbbl(nbbl,2)=ltran(i,3)
      lbbl(nbbl,3)=ltran(i,2)
      goto 7
      endif
    5 continue
    7 continue
c
      nbbloc=0
      do 8 i=1,nbbl-1
      do 9 j=i+1,nbbl
      if((lbbl(i,1).ne.me.and.lbbl(j,1).eq.me).or.
     + (lbbl(i,1).gt.lbbl(j,1).and.lbbl(i,1).ne.me).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).gt.lbbl(j,2)).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).eq.lbbl(j,2).and.
     + lbbl(i,3).gt.lbbl(j,3))) then
      ltemp=lbbl(i,1)
      lbbl(i,1)=lbbl(j,1)
      lbbl(j,1)=ltemp
      ltemp=lbbl(i,2)
      lbbl(i,2)=lbbl(j,2)
      lbbl(j,2)=ltemp
      ltemp=lbbl(i,3)
      lbbl(i,3)=lbbl(j,3)
      lbbl(j,3)=ltemp
      endif
    9 continue
      if(lbbl(i,1).eq.me) nbbloc=i
    8 continue
      if(lbbl(nbbl,1).eq.me) nbbloc=nbbl
c
      return
      end
      logical function sp_diffbb(boxsiz,ibxrng)
c
      implicit none
c
#include "sp_common.fh"
c
      integer ibxrng(maxbox,3,2)
      real*8 boxsiz(maxbox,3)
      integer i,ib,j,k,num
      real*8 rd
c
c     temporarily skip for non-rectangular boxes
c
      if(nbxtyp.eq.1) then
      sp_diffbb=.true.
      return
      endif
c
c     initialize
c
      if(nbbdif.lt.0) then
      do 1 k=1,2
      do 2 j=1,3
      do 3 i=1,maxbox
      ibxrng(i,j,k)=-1
    3 continue
    2 continue
    1 continue
      endif
c
      nbbdif=0
c
      do 4 i=1,nbx
      num=0
      rd=zero
      do 5 j=1,nbx
      ib=i+1
      if(ib.gt.nbx) ib=ib-nbx
      rd=rd+boxsiz(ib,1)
      if(rd.gt.rlong) goto 6
      num=num+1
    5 continue
    6 continue
      if(num.ne.ibxrng(i,1,1)) then
      ibxrng(i,1,1)=num
      nbbdif=nbbdif+1
      endif
      num=0
      rd=zero
      do 7 j=1,nbx
      ib=i-1
      if(ib.lt.1) ib=ib+nbx
      rd=rd+boxsiz(ib,1)
      if(rd.gt.rlong) goto 8
      num=num+1
    7 continue
    8 continue
      if(num.ne.ibxrng(i,1,2)) then
      ibxrng(i,1,2)=num
      nbbdif=nbbdif+1
      endif
    4 continue
c
      do 9 i=1,nby
      num=0
      rd=zero
      do 10 j=1,nby
      ib=i+1
      if(ib.gt.nby) ib=ib-nby
      rd=rd+boxsiz(ib,2)
      if(rd.gt.rlong) goto 11
      num=num+1
   10 continue
   11 continue
      if(num.ne.ibxrng(i,2,1)) then
      ibxrng(i,2,1)=num
      nbbdif=nbbdif+1
      endif
      num=0
      rd=zero
      do 12 j=1,nby
      ib=i-1
      if(ib.lt.1) ib=ib+nby
      rd=rd+boxsiz(ib,2)
      if(rd.gt.rlong) goto 13
      num=num+1
   12 continue
   13 continue
      if(num.ne.ibxrng(i,2,2)) then
      ibxrng(i,2,2)=num
      nbbdif=nbbdif+1
      endif
    9 continue
c
      do 14 i=1,nbz
      num=0
      rd=zero
      do 15 j=1,nbz
      ib=i+1
      if(ib.gt.nbz) ib=ib-nbz
      rd=rd+boxsiz(ib,3)
      if(rd.gt.rlong) goto 16
      num=num+1
   15 continue
   16 continue
      if(num.ne.ibxrng(i,3,1)) then
      ibxrng(i,3,1)=num
      nbbdif=nbbdif+1
      endif
      num=0
      rd=zero
      do 17 j=1,nbz
      ib=i-1
      if(ib.lt.1) ib=ib+nbz
      rd=rd+boxsiz(ib,3)
      if(rd.gt.rlong) goto 18
      num=num+1
   17 continue
   18 continue
      if(num.ne.ibxrng(i,3,2)) then
      ibxrng(i,3,2)=num
      nbbdif=nbbdif+1
      endif
   14 continue
c
      sp_diffbb=nbbdif.ne.0
c
      return
      end
