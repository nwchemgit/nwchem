      subroutine cf_weight(numw,nums,isat,isml)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums,isat(msa),isml(msa)
c
      call cf_wbox(numw,nums,int_mb(i_iwa),isat,isml,
     + dbl_mb(i_wsm),dbl_mb(i_mas))
c
      return
      end
      subroutine cf_volume()
c
      implicit none
c
#include "cf_common.fh"
c
      integer i
c
c     recalculate half boxsize and volume
c
c     reactangular box type
c
      if(nbxtyp.eq.0) volume=box(1)*box(2)*box(3)
c
c     non-rectangular box type
c
      if(nbxtyp.eq.1) then
      do 1 i=1,3
      box(i)=sqrt(vlat(1,i)*vlat(1,i)+vlat(2,i)*vlat(2,i)+
     + vlat(3,i)*vlat(3,i))
    1 continue
      volume=vlat(1,1)*(vlat(2,2)*vlat(3,3)-vlat(2,3)*vlat(3,2))+
     +       vlat(1,2)*(vlat(2,3)*vlat(3,1)-vlat(2,1)*vlat(3,3))+
     +       vlat(1,3)*(vlat(2,1)*vlat(3,2)-vlat(2,2)*vlat(3,1))
      endif
c
      if(nbxtyp.eq.2) then
      volume=half*box(1)*box(2)*box(3)
      endif
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      return
      end
      subroutine cf_wbox(numw,nums,iwat,isat,isml,wsm,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer numw,nums,iwat(mwa),isat(msa),isml(msa)
      real*8 wsm(msm),wgt(mat,6)
c
      integer i
c
      wbox=zero
c
      if(numw.gt.0) then
      do 1 i=1,nwa
      wbox=wbox+wgt(iwat(i),iset)
    1 continue
      wbox=dble(numw)*wbox
      endif
c
      do 2 i=1,msm
      wsm(i)=zero
    2 continue
c
      if(nums.gt.0) then
      do 3 i=1,nums
      wbox=wbox+wgt(isat(i),iset)
      wsm(isml(i))=wsm(isml(i))+wgt(isat(i),iset)
    3 continue
      endif
c
      if(np.gt.1) then
      call ga_dgop(mcf_30,wbox,1,'+')
      if(msm.gt.0) call ga_dgop(mcf_32,wsm,msm,'+')
      endif
c
      return
      end
      subroutine cf_cenmas(numw,xw,xwm,nums,isat,isml,xs,xsm)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      real*8 xw(mwm,3,mwa),xwm(mwm,3)
      real*8 xs(msa,3),xsm(msm,3)
      integer isat(msa),isml(msa)
c
      if(.not.ma_verify_allocator_stuff()) stop '00'
      call cf_comw(xw,xwm,1,numw)
      if(.not.ma_verify_allocator_stuff()) stop '01'
      call cf_coms(isat,isml,xs,xsm,1,nums)
      if(.not.ma_verify_allocator_stuff()) stop '02'
c
      return
      end
      subroutine cf_comw(xw,xwm,iwfr,iwto)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer iwfr,iwto
      real*8 xw(mwm,3,mwa),xwm(mwm,3)
c
      call cf_cmw(xw,xwm,iwfr,iwto,int_mb(i_iwa),dbl_mb(i_mas))
c
      return
      end
      subroutine cf_cmw(xw,xwm,iwfr,iwto,iwa,wgt)
c
      implicit none
c
#include "cf_common.fh"
c
      integer iwfr,iwto,iwa(mwa)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),wgt(mat,6)
c
      integer i,j,k
      real*8 factor
c
      if(iwfr.eq.0.or.iwfr.gt.iwto) return
c
      factor=zero
      do 1 i=1,mwa
      factor=factor+wgt(iwa(i),iset)
    1 continue
      factor=one/factor
c
      do 2 j=1,3
      do 3 i=iwfr,iwto
      xwm(i,j)=zero
      do 4 k=1,mwa
      xwm(i,j)=xwm(i,j)+wgt(iwa(k),iset)*xw(i,j,k)
    4 continue
      xwm(i,j)=factor*xwm(i,j)
    3 continue
    2 continue
c
      return
      end
      subroutine cf_coms(isat,isml,xs,xsm,isfr,isto)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer isfr,isto
      integer isat(msa),isml(msa)
      real*8 xs(msa,3),xsm(msm,3)
c
      call cf_cms(isat,isml,xs,xsm,isfr,isto,dbl_mb(i_wsm),
     + dbl_mb(i_mas))
c
      return
      end
      subroutine cf_cms(isat,isml,xs,xsm,isfr,isto,wsm,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer isfr,isto,isat(msa),isml(msa)
      real*8 xs(msa,3),xsm(msm,3),wsm(msm),wgt(mat,6)
c
      integer i,j
c
      do 1 i=1,msm
      xsm(i,1)=zero
      xsm(i,2)=zero
      xsm(i,3)=zero
    1 continue
c
      if(isfr.gt.0.and.isfr.le.isto) then
c
      do 2 j=1,3
      do 3 i=isfr,isto
      xsm(isml(i),j)=xsm(isml(i),j)+wgt(isat(i),iset)*xs(i,j)
    3 continue
      do 4 i=1,msm
      xsm(i,j)=xsm(i,j)/wsm(i)
    4 continue
    2 continue
c
      endif
c
      if(np.gt.0) call ga_dgop(mcf_33,xsm,3*msm,'+')
c
      return
      end
      subroutine cf_comb(numw,iwat,xwm,nums,isat,xs,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer numw,nums
      real*8 xwm(mwm,3),ww,xs(msa,3)
      integer iwat(msa),isat(msa)
      real*8 wgt(mat,6)
c
      integer iwm,iwa,isa
c
      xbm(1)=zero
      xbm(2)=zero
      xbm(3)=zero
      if(numw.gt.0) then
      ww=zero
      do 1 iwa=1,nwa
      ww=ww+wgt(iwat(iwa),iset)
    1 continue
      do 2 iwm=1,numw
      xbm(1)=xbm(1)+ww*xwm(iwm,1)
      xbm(2)=xbm(2)+ww*xwm(iwm,2)
      xbm(3)=xbm(3)+ww*xwm(iwm,3)
    2 continue
      endif
      if(nums.gt.0) then
      do 3 isa=1,nums
      ww=wgt(isat(isa),iset)
      xbm(1)=xbm(1)+ww*xs(isa,1)
      xbm(2)=xbm(2)+ww*xs(isa,2)
      xbm(3)=xbm(3)+ww*xs(isa,3)
    3 continue
      endif
      xbm(1)=xbm(1)/wbox
      xbm(2)=xbm(2)/wbox
      xbm(3)=xbm(3)/wbox
c
c     globally accumulate box center of mass coordinates
c
      if(np.gt.1) call ga_dgop(mcf_31,xbm,3,'+')
c
      return
      end
      subroutine cf_mdstep(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      integer iwdt(mwm),isdt(msa),isat(msa)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3),fs(msa,3)
c
      call cf_tscale()
c
      if(mdalgo.eq.1) then
      call cf_lpfrog(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,int_mb(i_iwa),
     + dbl_mb(i_mas))
      elseif(mdalgo.eq.2) then
      else
      call errquit('Unknown integration algorithm',mdalgo)
      endif
c
      return
      end
      subroutine cf_tscale()
c
      implicit none
c
#include "cf_common.fh"
c
      scalew=one
      scales=one
      if(itscal.eq.1) then
      if(temp.gt.small) scalew=sqrt(one+(tmpext/temp-one)*tstep/tmwrlx)
      scales=scalew
      else
      if(tempw.gt.small)
     + scalew=sqrt(one+(tmpext/tempw-one)*tstep/tmwrlx)
      if(temps.gt.small)
     + scales=sqrt(one+(tmpext/temps-one)*tstep/tmsrlx)
      endif
c
      return
      end
      subroutine cf_lpfrog(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,iwat,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      integer iwdt(mwm),isdt(msa),isat(msa)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3),fs(msa,3)
      integer iwat(mwa)
      real*8 wgt(mat,6)
c
      integer iwa,iwm,isa
      real*8 factor
c
      if(numw.gt.0) then
      do 1 iwa=1,mwa
      factor=tstep/wgt(iwat(iwa),iset)
      do 2 iwm=1,numw
      yw(iwm,1,iwa)=xw(iwm,1,iwa)
      yw(iwm,2,iwa)=xw(iwm,2,iwa)
      yw(iwm,3,iwa)=xw(iwm,3,iwa)
      if(iand(iwdt(iwm),mfixed).eq.lfixed) then
      vw(iwm,1,iwa)=zero
      vw(iwm,2,iwa)=zero
      vw(iwm,3,iwa)=zero
      vwt(iwm,1,iwa)=zero
      vwt(iwm,2,iwa)=zero
      vwt(iwm,3,iwa)=zero
      else
      vwt(iwm,1,iwa)=vw(iwm,1,iwa)
      vwt(iwm,2,iwa)=vw(iwm,2,iwa)
      vwt(iwm,3,iwa)=vw(iwm,3,iwa)
      vw(iwm,1,iwa)=scalew*(vw(iwm,1,iwa)+fw(iwm,1,iwa)*factor)
      vw(iwm,2,iwa)=scalew*(vw(iwm,2,iwa)+fw(iwm,2,iwa)*factor)
      vw(iwm,3,iwa)=scalew*(vw(iwm,3,iwa)+fw(iwm,3,iwa)*factor)
      xw(iwm,1,iwa)=yw(iwm,1,iwa)+vw(iwm,1,iwa)*tstep
      xw(iwm,2,iwa)=yw(iwm,2,iwa)+vw(iwm,2,iwa)*tstep
      xw(iwm,3,iwa)=yw(iwm,3,iwa)+vw(iwm,3,iwa)*tstep
      endif
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      do 3 isa=1,nums
      factor=tstep/wgt(isat(isa),iset)
      ys(isa,1)=xs(isa,1)
      ys(isa,2)=xs(isa,2)
      ys(isa,3)=xs(isa,3)
      if(iand(isdt(isa),mfixed).eq.lfixed) then
      vs(isa,1)=zero
      vs(isa,2)=zero
      vs(isa,3)=zero
      vst(isa,1)=zero
      vst(isa,2)=zero
      vst(isa,3)=zero
      else
      vst(isa,1)=vs(isa,1)
      vst(isa,2)=vs(isa,2)
      vst(isa,3)=vs(isa,3)
      vs(isa,1)=scales*(vs(isa,1)+fs(isa,1)*factor)
      vs(isa,2)=scales*(vs(isa,2)+fs(isa,2)*factor)
      vs(isa,3)=scales*(vs(isa,3)+fs(isa,3)*factor)
      xs(isa,1)=xs(isa,1)+vs(isa,1)*tstep
      xs(isa,2)=xs(isa,2)+vs(isa,2)*tstep
      xs(isa,3)=xs(isa,3)+vs(isa,3)*tstep
      endif
    3 continue
      endif
c
      return
      end
      subroutine cf_fold(numw,xw,xwm,nums,isat,isml,xs,xsm)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3),xsm(msm,3)
      integer isat(msa),isml(msa)
c
      if(.not.lscr) then
      if(.not.ma_push_get(mt_dbl,3*mscr,'s3r1',l_s3r1,i_s3r1))
     + call errquit('Failed to allocate scratch array s3r1',0)
      endif
c
      call cf_fld(numw,xw,xwm,nums,isat,isml,xs,xsm,
     + dbl_mb(i_s3r1))
c
      if(.not.lscr) then
      if(.not.ma_pop_stack(l_s3r1))
     + call errquit('Failed to deallocate scratch array s3r1',0)
      endif
c
      return
      end
      subroutine cf_fld(numw,xw,xwm,nums,isat,isml,xs,xsm,rwx)
c
      implicit none
c
#include "cf_common.fh"
c
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3),xsm(msm,3)
      real*8 rwx(mscr,3)
      integer isat(msa),isml(msa)
c
      integer iwm,iwa,isa,ism,ix,numw,nums
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
c     vacuo simulations do not need folding
c
      if(npbtyp.eq.0) return
c
      if(numw.gt.0) then
      call cf_pbc(0,xwm,mwm,rwx,mscr,0,1,numw)
      do 1 ix=1,3
      do 2 iwm=1,numw
      xwm(iwm,ix)=xwm(iwm,ix)-rwx(iwm,ix)
    2 continue
      do 3 iwa=1,nwa
      do 4 iwm=1,numw
      xw(iwm,ix,iwa)=xw(iwm,ix,iwa)-rwx(iwm,ix)
    4 continue
    3 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      call cf_pbc(0,xsm,nsm,rwx,mscr,0,1,nsm)
      do 5 ism=1,nsm
      do 6 ix=1,3
      xsm(ism,ix)=xsm(ism,ix)-rwx(ism,ix)
      do 7 isa=1,nums
      if(isml(isa).eq.ism) then
      xs(isa,ix)=xs(isa,ix)-rwx(ism,ix)
      endif
    7 continue
    6 continue
    5 continue
      call cf_pbc(0,xs,msa,rwx,mscr,0,1,nums)
      do 8 isa=1,nums
      if(isml(isa).eq.0) then
      xs(isa,1)=xs(isa,1)-rwx(isa,1)
      xs(isa,2)=xs(isa,2)-rwx(isa,2)
      xs(isa,3)=xs(isa,3)-rwx(isa,3)
      endif
    8 continue
      endif
c
      return
      end
      subroutine cf_final(xw,xwm,yw,vw,vwt,numw,
     + xs,xsm,ys,vs,vst,isat,ismol,isdyn,zs,nums,
     + bx,prs,prsscl,tmp)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      real*8 prs,prsscl,pscalm
      real*8 xw(mwm,3,mwa),xwm(mwm,3),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa)
      real*8 xs(msa,3),xsm(msm,3),ys(msa,3),vs(msa,3),vst(msa,3)
      real*8 zs(msf,3,2),bx(3),tmp
      integer isat(msa),ismol(msa),isdyn(msa)
      integer numw,nums
c
      integer i,j,k
      real*8 z,prssc3,rdata(4),xwcd,xscd
c
      if(ltwin) then
      do 8 i=1,3
      zw(i,1)=zw(i,1)+zw(i,2)
    8 continue
      do 9 i=1,nsf
      do 10 j=1,3
      zs(i,j,1)=zs(i,j,1)+zs(i,j,2)
   10 continue
    9 continue
      endif
c
      z=zw(1,1)+zw(2,1)+zw(3,1)
      do 1 i=1,nsf
      z=z+zs(i,1,1)+zs(i,2,1)+zs(i,3,1)
    1 continue
c
      call ga_dgop(mcf_30,z,1,'+')
c
      volume=box(1)*box(2)*box(3)
c
      pres=(two*(ewkc+eskc)-z)/(three*volume*avoga)
      prs=pres
c
      prsscl=one
      if(ipscal.gt.0) then
c
      prssc3=facpsc*(pres-prsext)+one
      if(prssc3.lt.small) prssc3=small
      prsscl=prssc3**third
c
c     scale the box dimensions
c
      box(1)=prsscl*box(1)
      box(2)=prsscl*box(2)
      box(3)=prsscl*box(3)
c
      if(np.gt.1) then
      rdata(1)=box(1)
      rdata(2)=box(2)
      rdata(3)=box(3)
      rdata(4)=prsscl
      call ga_brdcst(mcf_30,rdata,ma_sizeof(mt_dbl,4,mt_byte),0)
      box(1)=rdata(1)
      box(2)=rdata(2)
      box(3)=rdata(3)
      prsscl=rdata(4)
      endif
c
      call cf_volume()
c
      pscalm=prsscl-one
c
      if(numw.gt.0) then
      do 2 j=1,3
      do 3 i=1,numw
      xwcd=pscalm*xwm(i,j)
      xwm(i,j)=xwm(i,j)+xwcd
      do 4 k=1,nwa
      xw(i,j,k)=xw(i,j,k)+xwcd
    4 continue
    3 continue
    2 continue
      endif
c
      if(nums.gt.0) then
      do 5 j=1,3
      do 6 k=1,nsm
      xscd=pscalm*xsm(k,j)
      xsm(k,j)=xsm(k,j)+xscd
      do 7 i=1,nums
      if(ismol(i).eq.k.and.
     + iand(isdyn(i),mfixed).ne.lfixed) xs(i,j)=xs(i,j)+xscd
    7 continue
    6 continue
    5 continue
      endif
c
      endif
c
      call cf_vw(xw,yw,vw,vwt,numw,int_mb(i_iwa),dbl_mb(i_mas))
      call cf_vs(xs,ys,vs,vst,nums,isat,ismol,dbl_mb(i_mas))
c
      if(np.gt.1) then
      call ga_dgop(mcf_33,ewk,1,'+')
      call ga_dgop(mcf_34,esk,1,'+')
      call ga_dgop(mcf_35,ewkc,1,'+')
      endif
c
      tempw=factmw*ewk
      temps=factms*esk
      temp=factmp*(ewk+esk)
c
      tmp=temp
c
      return
      end
      subroutine cf_vw(xw,yw,vw,vwt,numw,iwat,wgt)
c
      implicit none
c
#include "cf_common.fh"
c
      integer numw
      integer iwat(mwa)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),vw(mwm,3,mwa),vwt(mwm,3,mwa)
      real*8 wgt(mat,6)
c
      integer i,j,k
      real*8 vvsum,ww
c
      ewk=zero
      if(numw.gt.0) then
      do 1 k=1,mwa
      vvsum=zero
      if(mdalgo.eq.1) then
      do 2 j=1,3
      do 3 i=1,numw
      vwt(i,j,k)=vw(i,j,k)
      vw(i,j,k)=(xw(i,j,k)-yw(i,j,k))*tstepi
      vwt(i,j,k)=half*(vw(i,j,k)+vwt(i,j,k))
      vvsum=vvsum+vwt(i,j,k)*vwt(i,j,k)
    3 continue
    2 continue
      else
      do 4 j=1,3
      do 5 i=1,numw
      vw(i,j,k)=(xw(i,j,k)-yw(i,j,k))*tstepi
      vvsum=vvsum+vw(i,j,k)*vw(i,j,k)
    5 continue
    4 continue
      endif
      ewk=ewk+wgt(iwat(k),iset)*vvsum
    1 continue
      ewk=half*ewk
      endif
c
      ewkc=zero
      if(numw.gt.0) then
      ww=zero
      do 6 i=1,mwa
      ww=ww+wgt(iwat(i),iset)
    6 continue
      if(mdalgo.eq.1) then
      do 7 j=1,3
      do 8 i=1,numw
      vvsum=zero
      do 9 k=1,mwa
      vvsum=vvsum+wgt(iwat(k),iset)*vwt(i,j,k)
    9 continue
      ewkc=ewkc+vvsum*vvsum
    8 continue
    7 continue
      else
      do 10 j=1,3
      do 11 i=1,numw
      vvsum=zero
      do 12 k=1,mwa
      vvsum=vvsum+wgt(iwat(k),iset)*vw(i,j,k)
   12 continue
      ewkc=ewkc+vvsum*vvsum
   11 continue
   10 continue
      endif
      ewkc=half*ewkc/ww
      endif
c
      return
      end
      subroutine cf_vs(xs,ys,vs,vst,nums,isat,ismol,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      integer isat(msa),ismol(msa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3),wgt(mat,6)
c
      integer i,j,k
      real*8 vv(4)
c
      esk=zero
      if(nums.gt.0) then
      do 1 j=1,3
      if(mdalgo.eq.1) then
      do 2 i=1,nums
      vst(i,j)=vs(i,j)
      vs(i,j)=(xs(i,j)-ys(i,j))*tstepi
      vst(i,j)=half*(vs(i,j)+vst(i,j))
      esk=esk+wgt(isat(i),iset)*vst(i,j)*vst(i,j)
    2 continue
      else
      do 3 i=1,nums
      vs(i,j)=(xs(i,j)-ys(i,j))*tstepi
      esk=esk+wgt(isat(i),iset)*vs(i,j)*vs(i,j)
    3 continue
      endif
    1 continue
      esk=half*esk
      endif
c
      eskc=zero
      do 4 k=1,nsm
      vv(1)=zero
      vv(2)=zero
      vv(3)=zero
      vv(4)=zero
      if(nums.gt.0) then
      do 5 i=1,nums
      if(ismol(i).eq.k) then
      vv(4)=vv(4)+wgt(isat(i),iset)
      if(mdalgo.eq.1) then
      do 6 j=1,3
      vv(1)=vv(1)+vst(i,j)
    6 continue
      else
      do 7 j=1,3
      vv(1)=vv(1)+vs(i,j)
    7 continue
      endif
      endif
    5 continue
      endif
      if(np.gt.1) call ga_dgop(mcf_34,vv,4,'+')
      eskc=eskc+half*vv(4)*(vv(1)*vv(1)+vv(2)*vv(2)+vv(3)*vv(3))
    4 continue
c
      return
      end
