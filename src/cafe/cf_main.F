      subroutine cf_weight(numw,nums,isat,isml,wbx)
c $Id: cf_main.F,v 1.133 2002-02-27 14:44:54 d3j191 Exp $
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums,isat(msa),isml(msa)
      real*8 wbx
c
      call cf_wbox(numw,nums,int_mb(i_iwa),isat,isml,
     + dbl_mb(i_wsm),dbl_mb(i_mas))
c
      wbx=wbox
c
      return
      end
      subroutine cf_mass(ww,ws,isat,nums)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer nums,isat(msa)
      real*8 ww(mwa),ws(msa)
c
      call cf_m(ww,ws,isat,nums,int_mb(i_iwa),dbl_mb(i_mas))
c
      return
      end
      subroutine cf_m(ww,ws,isat,nums,iwat,wgt)
c
      implicit none
c
#include "cf_common.fh"
c
      integer nums,isat(msa)
      real*8 ww(mwa),ws(msa)
      integer iwat(msa)
      real*8 wgt(mat,6)
c
      integer i
c
      do 1 i=1,nwa
      ww(i)=wgt(iwat(i),iset)
    1 continue
      do 2 i=1,nums
      ws(i)=wgt(isat(i),iset)
    2 continue
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
      call ga_dgop(mcf_34,wbox,1,'+')
      if(msm.gt.0) call ga_dgop(mcf_35,wsm,msm,'+')
      endif
c
      return
      end
      subroutine cf_cenmas(numw,xw,xwm,nums,isat,isml,xs,xsm,gsm)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      real*8 xw(mwm,3,mwa),xwm(mwm,3)
      real*8 xs(msa,3),xsm(msm,3),gsm(msm,4,2)
      integer isat(msa),isml(msa)
c
      call cf_comw(xw,xwm,1,numw)
      call cf_coms(isat,isml,xs,xsm,gsm,1,nums)
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
      subroutine cf_coms(isat,isml,xs,xsm,gsm,isfr,isto)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer isfr,isto
      integer isat(msa),isml(msa)
      real*8 xs(msa,3),xsm(msm,3),gsm(msm,4,2)
c
      call cf_cms(isat,isml,xs,xsm,gsm,isfr,isto,dbl_mb(i_wsm),
     + dbl_mb(i_mas))
c
      return
      end
      subroutine cf_cms(isat,isml,xs,xsm,gsm,isfr,isto,wsm,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer isfr,isto,isat(msa),isml(msa)
      real*8 xs(msa,3),xsm(msm,3),gsm(msm,4,2),wsm(msm),wgt(mat,6)
c
      integer i,j
      real*8 dist2
c
      do 1 i=1,msm
      xsm(i,1)=zero
      xsm(i,2)=zero
      xsm(i,3)=zero
      gsm(i,1,1)=zero
      gsm(i,2,1)=zero
      gsm(i,3,1)=zero
      gsm(i,4,1)=zero
      gsm(i,1,2)=zero
      gsm(i,2,2)=zero
      gsm(i,3,2)=zero
      gsm(i,4,2)=zero
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
      if(np.gt.0) call ga_dgop(mcf_36,xsm,3*msm,'+')
c
      if(isfr.gt.0.and.isfr.le.isto) then
c
      do 5 i=isfr,isto
      dist2=zero
      do 6 j=1,3
      dist2=dist2+(xs(i,j)-xsm(isml(i),j))**2
    6 continue
      do 7 j=1,3
      gsm(isml(i),j,1)=gsm(isml(i),j,1)+wgt(isat(i),j)*dist2
      gsm(isml(i),j,2)=gsm(isml(i),j,2)+wgt(isat(i),j)
    7 continue
    5 continue
c
      endif
c
      if(np.gt.0) call ga_dgop(mcf_67,gsm,8*msm,'+')
c
      do 8 j=1,4
      do 9 i=1,nsm
      gsm(i,j,2)=sqrt(gsm(i,j,1)/gsm(i,j,2))
    9 continue
    8 continue
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
      if(np.gt.1) call ga_dgop(mcf_37,xbm,3,'+')
c
      return
      end
      subroutine cf_restor(xw,yw,fw,vw,numw,xs,ys,fs,vs,nums)
c
      implicit none
c
#include "cf_common.fh"
c
      integer numw,nums
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),fw(mwm,3,mwa),vw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),fs(msa,3),vs(msa,3)
c
      integer i,j
c
      if(numw.gt.0) then
      do 1 j=1,nwa
      do 2 i=1,numw
      xw(i,1,j)=yw(i,1,j)
      xw(i,2,j)=yw(i,2,j)
      xw(i,3,j)=yw(i,3,j)
      fw(i,1,j)=vw(i,1,j)
      fw(i,2,j)=vw(i,2,j)
      fw(i,3,j)=vw(i,3,j)
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      do 4 i=1,nums
      xs(i,1)=ys(i,1)
      xs(i,2)=ys(i,2)
      xs(i,3)=ys(i,3)
      fs(i,1)=vs(i,1)
      fs(i,2)=vs(i,2)
      fs(i,3)=vs(i,3)
    4 continue
      endif
c
      return
      end
      subroutine cf_mdstep(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,iwgmn,isgan,tmpx)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      integer iwdt(mwm),isdt(msa),isat(msa),isgan(msa),iwgmn(mwm)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3),fs(msa,3)
      real*8 tmpx
c
      tmpext=tmpx
c
      call cf_tscale()
c
      if(mdalgo.eq.1) then
      call cf_lpfrog(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,int_mb(i_iwa),
     + dbl_mb(i_mas),iwgmn,isgan)
      elseif(mdalgo.eq.2) then
      call cf_lpfrog_bc(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,int_mb(i_iwa),
     + dbl_mb(i_mas),iwgmn,isgan)
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
      elseif(itscal.gt.1) then
      if(tempw.gt.small)
     + scalew=sqrt(one+(tmpext/tempw-one)*tstep/tmwrlx)
      if(temps.gt.small)
     + scales=sqrt(one+(tmpext/temps-one)*tstep/tmsrlx)
      endif
c
      return
      end
      subroutine cf_lpfrog(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,iwat,wgt,iwgmn,isgan)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      integer iwdt(mwm),isdt(msa),isat(msa),isgan(msa),iwgmn(mwm)
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
      subroutine cf_lpfrog_bc(iwdt,xw,yw,vw,vwt,fw,numw,
     + isdt,isat,xs,ys,vs,vst,fs,nums,iwat,wgt,iwgmn,isgan)
c
c     leap frog variant by Brown and Clarke, Mol.Phys., 51, 1243 (1984)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer numw,nums
      integer iwdt(mwm),isdt(msa),isat(msa),isgan(msa),iwgmn(mwm)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3),fs(msa,3)
      integer iwat(mwa)
      real*8 wgt(mat,6)
c
      integer iwa,iwm,isa
      real*8 factor,chi,ekin(2)
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
      vwt(iwm,1,iwa)=vw(iwm,1,iwa)+half*fw(iwm,1,iwa)*factor
      vwt(iwm,2,iwa)=vw(iwm,2,iwa)+half*fw(iwm,2,iwa)*factor
      vwt(iwm,3,iwa)=vw(iwm,3,iwa)+half*fw(iwm,3,iwa)*factor
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
      vst(isa,1)=vs(isa,1)+half*fs(isa,1)*factor
      vst(isa,2)=vs(isa,2)+half*fs(isa,2)*factor
      vst(isa,3)=vs(isa,3)+half*fs(isa,3)*factor
      endif
    3 continue
      endif
c
      call cf_tw(vwt,numw,iwat,dbl_mb(i_mas),ekin(1))
      call cf_ts(vst,nums,isat,dbl_mb(i_mas),ekin(2))
c
      call ga_dgop(mcf_65,ekin,2,'+')
c
      chi=sqrt(tmpext/(factmp*(ekin(1)+ekin(2))))
c
      if(numw.gt.0) then
      do 4 iwa=1,mwa
      factor=tstep/wgt(iwat(iwa),iset)
      do 5 iwm=1,numw
      if(iand(iwdt(iwm),mfixed).ne.lfixed) then
      vw(iwm,1,iwa)=(two*chi-one)*vw(iwm,1,iwa)+chi*fw(iwm,1,iwa)*factor
      vw(iwm,2,iwa)=(two*chi-one)*vw(iwm,2,iwa)+chi*fw(iwm,2,iwa)*factor
      vw(iwm,3,iwa)=(two*chi-one)*vw(iwm,3,iwa)+chi*fw(iwm,3,iwa)*factor
      xw(iwm,1,iwa)=yw(iwm,1,iwa)+vw(iwm,1,iwa)*tstep
      xw(iwm,2,iwa)=yw(iwm,2,iwa)+vw(iwm,2,iwa)*tstep
      xw(iwm,3,iwa)=yw(iwm,3,iwa)+vw(iwm,3,iwa)*tstep
      endif
    5 continue
    4 continue
      endif
c
      if(nums.gt.0) then
      do 6 isa=1,nums
      factor=tstep/wgt(isat(isa),iset)
      if(iand(isdt(isa),mfixed).ne.lfixed) then
      vs(isa,1)=(two*chi-one)*vs(isa,1)+chi*fs(isa,1)*factor
      vs(isa,2)=(two*chi-one)*vs(isa,2)+chi*fs(isa,2)*factor
      vs(isa,3)=(two*chi-one)*vs(isa,3)+chi*fs(isa,3)*factor
      xs(isa,1)=xs(isa,1)+vs(isa,1)*tstep
      xs(isa,2)=xs(isa,2)+vs(isa,2)*tstep
      xs(isa,3)=xs(isa,3)+vs(isa,3)*tstep
      endif
    6 continue
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
     + xs,xsm,ys,vs,vst,isat,ismol,isdyn,isfrc,zs,esk,nums,
     + bx,vl,prs,tmp,tmpw,tmps)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      real*8 pscalm
      real*8 xw(mwm,3,mwa),xwm(mwm,3),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa)
      real*8 xs(msa,3),xsm(msm,3),ys(msa,3),vs(msa,3),vst(msa,3)
      real*8 zs(msf,3,3,2),bx(3),prs,tmp,tmpw,tmps,esk(msf),vl(3,3)
      integer isat(msa),ismol(msa),isdyn(msa),isfrc(msa)
      integer numw,nums
c
      integer i,j,k
      real*8 prssc3,rdata(4),xwcd,xscd,eskin,ekinc
      real*8 pscx(3)
      real*8 pbuff(4,3,3)
      integer i_vsm,l_vsm
c
      call cf_vw(xw,yw,vw,vwt,numw,int_mb(i_iwa),dbl_mb(i_mas))
      if(.not.ma_push_get(mt_dbl,3*msm,'vsm',l_vsm,i_vsm))
     + call errquit('Failed to allocate scratch array vsm',0)
      call cf_vs(xs,ys,vs,vst,esk,nums,isat,ismol,isfrc,dbl_mb(i_wsm),
     + dbl_mb(i_mas),dbl_mb(i_vsm))
      if(.not.ma_pop_stack(l_vsm))
     + call errquit('Failed to deallocate scratch array vsm',0)
c
      if(ltwin) then
      do 1 i=1,3
      zw(1,i,1)=zw(1,i,1)+zw(1,i,2)
      zw(2,i,1)=zw(2,i,1)+zw(2,i,2)
      zw(3,i,1)=zw(3,i,1)+zw(3,i,2)
    1 continue
      do 2 i=1,nsf
      do 3 j=1,3
      zs(i,1,j,1)=zs(i,1,j,1)+zs(i,1,j,2)
      zs(i,2,j,1)=zs(i,2,j,1)+zs(i,2,j,2)
      zs(i,3,j,1)=zs(i,3,j,1)+zs(i,3,j,2)
    3 continue
    2 continue
      ep2(1)=ep2(1)+ep2(2)
      ep3(1)=ep3(1)+ep3(2)
      if(ntype.eq.3) then
      do 4 i=1,24
      deriv(i,1)=deriv(i,1)+deriv(i,2)
    4 continue
      endif
      endif
c
      z(1,1)=zw(1,1,1)
      z(1,2)=zw(1,2,1)
      z(1,3)=zw(1,3,1)
      z(2,1)=zw(2,1,1)
      z(2,2)=zw(2,2,1)
      z(2,3)=zw(2,3,1)
      z(3,1)=zw(3,1,1)
      z(3,2)=zw(3,2,1)
      z(3,3)=zw(3,3,1)
      do 5 i=1,nsf
      z(1,1)=z(1,1)+zs(i,1,1,1)
      z(2,1)=z(2,1)+zs(i,2,1,1)
      z(3,1)=z(3,1)+zs(i,3,1,1)
      z(1,2)=z(1,2)+zs(i,1,2,1)
      z(2,2)=z(2,2)+zs(i,2,2,1)
      z(3,2)=z(3,2)+zs(i,3,2,1)
      z(1,3)=z(1,3)+zs(i,1,3,1)
      z(2,3)=z(2,3)+zs(i,2,3,1)
      z(3,3)=z(3,3)+zs(i,3,3,1)
    5 continue
c
      z(1,1)=z(1,1)+half*vpme(1)
      z(2,2)=z(2,2)+half*vpme(4)
      z(3,3)=z(3,3)+half*vpme(6)
c
      ekc(1,1)=ewkc(1,1)+eskc(1,1)
      ekc(2,1)=ewkc(2,1)+eskc(2,1)
      ekc(3,1)=ewkc(3,1)+eskc(3,1)
      ekc(1,2)=ewkc(1,2)+eskc(1,2)
      ekc(2,2)=ewkc(2,2)+eskc(2,2)
      ekc(3,2)=ewkc(3,2)+eskc(3,2)
      ekc(1,3)=ewkc(1,3)+eskc(1,3)
      ekc(2,3)=ewkc(2,3)+eskc(2,3)
      ekc(3,3)=ewkc(3,3)+eskc(3,3)
c
      do 18 i=1,3
      do 19 j=1,3
      pbuff(1,i,j)=z(i,j)
      pbuff(2,i,j)=ekc(i,j)
      pbuff(3,i,j)=ewkc(i,j)
      pbuff(4,i,j)=eskc(i,j)
   19 continue
   18 continue
c
      call ga_dgop(mcf_38,pbuff,36,'+')
c
      do 20 i=1,3
      do 21 j=1,3
      z(i,j)=pbuff(1,i,j)
      ekc(i,j)=pbuff(2,i,j)
      ewkc(i,j)=pbuff(3,i,j)
      eskc(i,j)=pbuff(4,i,j)
   21 continue
   20 continue
c
c      call ga_dgop(mcf_38,z,9,'+')
c      call ga_dgop(mcf_70,ekc,9,'+')
c      call ga_dgop(mcf_71,ewkc,9,'+')
c      call ga_dgop(mcf_72,eskc,9,'+')
c
      volume=box(1)*box(2)*box(3)
c
      virial=z(1,1)+z(2,2)+z(3,3)
      ekinc=ekc(1,1)+ekc(2,2)+ekc(3,3)
c
      pres=0.0d0
      prx(1,1)=0.0d0
      prx(2,1)=0.0d0
      prx(3,1)=0.0d0
      prx(1,2)=0.0d0
      prx(2,2)=0.0d0
      prx(3,2)=0.0d0
      prx(1,3)=0.0d0
      prx(2,3)=0.0d0
      prx(3,3)=0.0d0
      if(lpress) then
      pres=(two*ekinc-virial)/(three*volume*avoga)
      prx(1,1)=(two*ekc(1,1)-z(1,1))/(three*volume*avoga)
      prx(2,1)=(two*ekc(2,1)-z(2,1))/(three*volume*avoga)
      prx(3,1)=(two*ekc(3,1)-z(3,1))/(three*volume*avoga)
      prx(1,2)=(two*ekc(1,2)-z(1,2))/(three*volume*avoga)
      prx(2,2)=(two*ekc(2,2)-z(2,2))/(three*volume*avoga)
      prx(3,2)=(two*ekc(3,2)-z(3,2))/(three*volume*avoga)
      prx(1,3)=(two*ekc(1,3)-z(1,3))/(three*volume*avoga)
      prx(2,3)=(two*ekc(2,3)-z(2,3))/(three*volume*avoga)
      prx(3,3)=(two*ekc(3,3)-z(3,3))/(three*volume*avoga)
      endif
c
      pscale=one
      if(ipscal.ne.0) then
c
      prssc3=facpsc*(pres-prsext)+one
      if(prssc3.lt.small) prssc3=small
      pscale=prssc3**third
      pscx(1)=pscale
      pscx(2)=pscale
      pscx(3)=pscale
      if(ipscal.eq.2) then
      do 17 i=1,3
      prssc3=facpsc*(three*prx(i,i)-prsext)+one
      if(prssc3.lt.small) prssc3=small
      pscx(i)=prssc3**third
   17 continue
      endif
      if(ipopt.eq.12) then
      pscx(1)=one
      pscx(2)=one
      endif
      if(ipopt.eq.3) pscx(3)=one
      pscale=pscx(1)*pscx(2)*pscx(3)
c
      if(np.gt.1) then
      rdata(1)=pscx(1)
      rdata(2)=pscx(2)
      rdata(3)=pscx(3)
      rdata(4)=pscale
      call ga_brdcst(mcf_39,rdata,ma_sizeof(mt_dbl,4,mt_byte),0)
      pscx(1)=rdata(1)
      pscx(2)=rdata(2)
      pscx(3)=rdata(3)
      pscale=rdata(4)
      endif
c
c     scale the box dimensions
c
      box(1)=pscx(1)*box(1)
      box(2)=pscx(2)*box(2)
      box(3)=pscx(3)*box(3)
c
      do 6 i=1,3
      do 7 j=1,3
      vlat(i,j)=sqrt(pscx(i)*pscx(j))*vlat(i,j)
    7 continue
    6 continue
c
      call cf_volume()
c
      if(numw.gt.0) then
      do 8 j=1,3
      pscalm=pscx(j)-one
      do 9 i=1,numw
      xwcd=pscalm*xwm(i,j)
      xwm(i,j)=xwm(i,j)+xwcd
      do 10 k=1,nwa
      xw(i,j,k)=xw(i,j,k)+xwcd
   10 continue
    9 continue
    8 continue
      endif
c
      if(nums.gt.0) then
      do 11 j=1,3
      pscalm=pscx(j)-one
      do 12 k=1,nsm
      xscd=pscalm*xsm(k,j)
      xsm(k,j)=xsm(k,j)+xscd
      do 13 i=1,nums
      if(ismol(i).eq.k.and.
     + iand(isdyn(i),mfixed).ne.lfixed) xs(i,j)=xs(i,j)+xscd
   13 continue
   12 continue
   11 continue
      endif
c
      endif
c
      if(np.gt.1) then
      call ga_dgop(mcf_40,ewk,1,'+')
      call ga_dgop(mcf_41,esk,msf,'+')
      call ga_dgop(mcf_42,ewkc,9,'+')
      endif
c
      eskin=zero
      do 14 i=1,nsf
      eskin=eskin+esk(i)
   14 continue
c
      tempw=factmw*ewk
      temps=factms*eskin
      temp=factmp*(ewk+eskin)
c
      prs=pres
      tmp=temp
      tmpw=tempw
      temps=temps
c
      do 15 i=1,3
      bx(i)=box(i)
      boxh(i)=half*box(i)
      do 16 j=1,3
      vl(i,j)=vlat(i,j)
   16 continue
   15 continue
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
      integer i,j,k,l
      real*8 vvsum,ww,visum,vjsum
c
      ewk=zero
      if(numw.gt.0) then
      do 1 k=1,mwa
      vvsum=zero
      if(mdalgo.eq.1) then
      do 2 j=1,3
      do 3 i=1,numw
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
      ewkc(1,1)=zero
      ewkc(2,1)=zero
      ewkc(3,1)=zero
      ewkc(1,2)=zero
      ewkc(2,2)=zero
      ewkc(3,2)=zero
      ewkc(1,3)=zero
      ewkc(2,3)=zero
      ewkc(3,3)=zero
      if(numw.gt.0) then
      ww=zero
      do 6 i=1,mwa
      ww=ww+wgt(iwat(i),iset)
    6 continue
      if(mdalgo.eq.1) then
      do 17 l=1,3
      do 7 j=1,3
      do 8 i=1,numw
      visum=zero
      vjsum=zero
      do 9 k=1,mwa
      visum=visum+wgt(iwat(k),iset)*vwt(i,j,k)
      vjsum=vjsum+wgt(iwat(k),iset)*vwt(i,l,k)
    9 continue
      ewkc(l,j)=ewkc(l,j)+visum*vjsum
    8 continue
    7 continue
   17 continue
      else
      do 110 l=1,3
      do 10 j=1,3
      do 11 i=1,numw
      visum=zero
      vjsum=zero
      do 12 k=1,mwa
      visum=visum+wgt(iwat(k),iset)*vw(i,j,k)
      vjsum=vjsum+wgt(iwat(k),iset)*vw(i,l,k)
   12 continue
      ewkc(l,j)=ewkc(l,j)+visum*vjsum
   11 continue
   10 continue
  110 continue
      endif
      ewkc(1,1)=half*ewkc(1,1)/ww
      ewkc(2,1)=half*ewkc(2,1)/ww
      ewkc(3,1)=half*ewkc(3,1)/ww
      ewkc(1,2)=half*ewkc(1,2)/ww
      ewkc(2,2)=half*ewkc(2,2)/ww
      ewkc(3,2)=half*ewkc(3,2)/ww
      ewkc(1,3)=half*ewkc(1,3)/ww
      ewkc(2,3)=half*ewkc(2,3)/ww
      ewkc(3,3)=half*ewkc(3,3)/ww
      endif
c
      return
      end
      subroutine cf_vs
     + (xs,ys,vs,vst,esk,nums,isat,ismol,isfrc,wsm,wgt,vsm)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      integer isat(msa),ismol(msa),isfrc(msa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),vst(msa,3)
      real*8 wgt(mat,6),vsm(msm,3)
      real*8 wsm(msm)
      real*8 esk(msf)
c
      integer i,j,k
      real*8 w,vv(4)
c
      do 11 i=1,msf
      esk(i)=zero
   11 continue
      if(nums.gt.0) then
      do 1 j=1,3
      if(mdalgo.eq.1) then
      do 2 i=1,nums
      vs(i,j)=(xs(i,j)-ys(i,j))*tstepi
      vst(i,j)=half*(vs(i,j)+vst(i,j))
      esk(isfrc(i))=esk(isfrc(i))+
     + half*wgt(isat(i),iset)*vst(i,j)*vst(i,j)
    2 continue
      else
      do 3 i=1,nums
      vs(i,j)=(xs(i,j)-ys(i,j))*tstepi
      esk(isfrc(i))=esk(isfrc(i))+half*wgt(isat(i),iset)*vs(i,j)*vs(i,j)
    3 continue
      endif
    1 continue
      endif
c
      eskc(1,1)=zero
      eskc(2,1)=zero
      eskc(3,1)=zero
      eskc(1,2)=zero
      eskc(2,2)=zero
      eskc(3,2)=zero
      eskc(1,3)=zero
      eskc(2,3)=zero
      eskc(3,3)=zero
      do 5 i=1,msm
      vsm(i,1)=zero
      vsm(i,2)=zero
      vsm(i,3)=zero
    5 continue
      do 4 i=1,nums
      k=ismol(i)
      w=wgt(isat(i),iset)
      if(mdalgo.eq.1) then
      do 6 j=1,3
      vsm(k,j)=vsm(k,j)+w*vst(i,j)
    6 continue
      else
      do 7 j=1,3
      vsm(k,j)=vsm(k,j)+w*vs(i,j)
    7 continue
      endif
    4 continue
      if(np.gt.1) call ga_dgop(mcf_43,vsm,3*msm,'+')
      if(me.eq.0) then
      do 8 i=1,nsm
      w=half/wsm(i) 
      eskc(1,1)=eskc(1,1)+w*vsm(i,1)*vsm(i,1)
      eskc(2,1)=eskc(2,1)+w*vsm(i,2)*vsm(i,1)
      eskc(3,1)=eskc(3,1)+w*vsm(i,3)*vsm(i,1)
      eskc(1,2)=eskc(1,2)+w*vsm(i,1)*vsm(i,2)
      eskc(2,2)=eskc(2,2)+w*vsm(i,2)*vsm(i,2)
      eskc(3,2)=eskc(3,2)+w*vsm(i,3)*vsm(i,2)
      eskc(1,3)=eskc(1,3)+w*vsm(i,1)*vsm(i,3)
      eskc(2,3)=eskc(2,3)+w*vsm(i,2)*vsm(i,3)
      eskc(3,3)=eskc(3,3)+w*vsm(i,3)*vsm(i,3)
    8 continue
      endif
c
      return
      end
      subroutine cf_tw(vw,numw,iwat,wgt,ekinw)
c
      implicit none
c
#include "cf_common.fh"
c
      integer numw
      integer iwat(mwa)
      real*8 vw(mwm,3,mwa)
      real*8 wgt(mat,6)
c
      integer i,j,k
      real*8 ekinwa,ekinw
c
      ekinw=zero
      if(numw.gt.0) then
      do 1 k=1,mwa
      ekinwa=zero
      do 2 j=1,3
      do 3 i=1,numw
      ekinwa=ekinwa+vw(i,j,k)*vw(i,j,k)
    3 continue
    2 continue
      ekinw=ekinw+wgt(iwat(k),iset)*ekinwa
    1 continue
      ekinw=half*ekinw
      endif
c
      return
      end
      subroutine cf_ts(vs,nums,isat,wgt,ekins)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      integer isat(msa)
      real*8 vs(msa,3),wgt(mat,6)
c
      integer i,j
      real*8 ekins
c
      ekins=zero
      if(nums.gt.0) then
      do 1 j=1,3
      do 2 i=1,nums
      ekins=ekins+wgt(isat(i),iset)*vs(i,j)*vs(i,j)
    2 continue
    1 continue
      ekins=half*ekins
      endif
c
      return
      end
      subroutine cf_proper(vol,tmp,tmpw,tmps,prs,tsw,tss,psc,
     + ewkin,nwwl,nwws,nswl,nsws,nssl,nsss,nshw,nshs,ept2,ept3,
     + edistr,epomef,vir,prest,virt,ekct,ubiasi)
c
      implicit none
c
#include "cf_common.fh"
c
      real*8 vol,tmp,tmpw,tmps,prs,tsw,tss,psc,ewkin,ept2,ept3
      real*8 edistr,epomef,vir,ubiasi
      real*8 prest(3,3),virt(3,3),ekct(3,3)
      integer nwwl,nwws,nswl,nsws,nssl,nsss,nshw,nshs
c
      integer i,j
c
      vol=volume
      tmp=temp
      tmpw=tempw
      tmps=temps
      prs=pres
      tsw=scalew
      tss=scales
      psc=pscale
      vir=virial
c
      ewkin=ewk
c
      nwwl=llww
      nwws=lsww
      nswl=llsw
      nsws=lssw
      nssl=llss
      nsss=lsss
c
      nshw=niterw
      nshs=niters
c
      ept2=ep2(1)
      ept3=ep3(1)
c
      edistr=edrs
      epomef=epmf
c
      ubiasi=ubias
c
      do 1 i=1,3
      do 2 j=1,3
      prest(i,j)=prx(i,j)
      virt(i,j)=z(i,j)
      ekct(i,j)=ekc(i,j)
    2 continue
    1 continue
c
      return
      end
      subroutine cf_fnorm(iwdt,fw,numw,isat,isdt,fs,nums,fnorm,fmax)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      real*8 fw(mwm,3,mwa),fs(msa,3),fnorm,fmax
      integer isat(msa),numw,nums
      integer isdt(msa),iwdt(mwm)
c
      call cf_fnrm(dbl_mb(i_mas),int_mb(i_iwa),iwdt,fw,numw,
     + isat,isdt,fs,nums,fnorm,fmax)
c
      return
      end
      subroutine cf_fnrm(wgt,iwat,iwdt,fw,numw,isat,isdt,fs,nums,
     + fnorm,fmax)
c
      implicit none
c
#include "cf_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      real*8 wgt(mat,6)
      real*8 fw(mwm,3,mwa),fs(msa,3),fnorm,fmax
      integer iwat(mwa),isat(msa),numw,nums
      integer isdt(msa),iwdt(mwm)
c
      integer i,j
      real*8 fm,factor
c
      fnorm=zero
      fmax=zero
c
      if(numw.gt.0) then
      do 1 j=1,nwa
      factor=(one/wgt(iwat(j),iset))*(one/wgt(iwat(j),iset))
      do 2 i=1,numw
      if(iand(iwdt(i),lfixed).eq.0) then
      fm=fw(i,1,j)**2+fw(i,2,j)**2+fw(i,3,j)**2
      fnorm=fnorm+fm
      fmax=max(fmax,factor*fm)
      endif
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      do 3 i=1,nums
      if(iand(isdt(i),lfixed).eq.0) then
      fm=fs(i,1)**2+fs(i,2)**2+fs(i,3)**2
      fnorm=fnorm+fm
      factor=(one/wgt(isat(i),iset))*(one/wgt(isat(i),iset))
      fmax=max(fmax,factor*fm)
      endif
    3 continue
      endif
c
      if(np.gt.1) then
      call ga_dgop(mcf_44,fnorm,1,'+')
      call ga_dgop(mcf_45,fmax,1,'max')
      endif
c
      fnorm=sqrt(fnorm)
      fmax=sqrt(fmax)
c
      return
      end
      subroutine cf_collapse(ncoll,fcoll,nums,numw,ismol,
     + xs,xsm,fs,xwm,xw,fw)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer ncoll,numw,nums
      integer ismol(msa)
      real*8 fcoll
      real*8 xsm(msm,3),xs(msa,3),fs(msa,3)
      real*8 xwm(mwm,3),xw(mwm,3,mwa),fw(mwm,3,mwa)
c
      integer i,j,k,ki,kf
c
      ki=1
      kf=3
      if(ncoll.eq.1) ki=3
      if(ncoll.eq.2) kf=2
c
      do 1 k=ki,kf
      do 2 i=1,nums
      fs(i,k)=fs(i,k)-fcoll*xsm(ismol(i),k)
    2 continue
      do 3 j=1,nwa
      do 4 i=1,numw
      fw(i,k,j)=fw(i,k,j)-fcoll*xwm(i,k)
    4 continue
    3 continue
    1 continue
c
      return
      end
      subroutine cf_extern(stime,nums,fs,isq,numw,fw)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer nums,numw
      integer isq(msa)
      real*8 stime,fs(msa,3),fw(mwm,3,mwa)
c
      call cf_external_field(stime,nums,fs,isq,numw,fw,int_mb(i_iwq),
     + dbl_mb(i_chg))
c
      return
      end
      subroutine cf_external_field(stime,nums,fs,isq,numw,fw,iwq,chg)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer nums,numw
      integer isq(msa),iwq(mwa)
      real*8 stime,fs(msa,3),fw(mwm,3,mwa),chg(mqt,mqp,6)
c
      integer i,j,k
      real*8 efield,q,df
c
      efield=field*cos(2.0d-6*pi*ffreq*stime)
c
      if(numw.gt.0) then
      do 1 k=1,nwa
      q=chg(iwq(k),1,iset)*efield
      do 2 j=1,3
      df=q*fvect(j)
      do 3 i=1,numw
      fw(i,j,k)=fw(i,j,k)+df
    3 continue
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      do 4 j=1,3
      do 5 i=1,nums
      fs(i,j)=fs(i,j)+chg(isq(i),1,iset)*efield*fvect(j)
    5 continue
    4 continue
      endif
c
      return
      end
      subroutine cf_multi(nums,xs,fs,isgan,isatt,isfrc,isdyn,isct1,ess,
     + lfnpmf,lrec,lpmfci)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer nums,lfnpmf
      real*8 xs(msa,3),fs(msa,3),ess(msf,msf,mpe,2)
      integer isgan(msa),isatt(msa),isfrc(msa),isdyn(msa),isct1(msa)
      logical lrec,lpmfci
c
      lpmfc=lpmfci
      edrs=zero
      epmf=zero
      ubias=zero
c
      if(mmult.gt.0) then
      call cf_mlist(nums,xs,fs,isgan,isfrc,isdyn,isct1,
     + mbt(2),int_mb(i_ibnd(2)),mht(2),int_mb(i_iang(2)),
     + mdt(2),int_mb(i_idih(2)),mit(2),int_mb(i_iimp(2)),
     + int_mb(i_ixmul),int_mb(i_imul),dbl_mb(i_xmul),dbl_mb(i_fmul),ess)
      endif
c
      if(ndrs.gt.0) then
      call cf_dlist(nums,xs,fs,isgan,int_mb(i_idrs),dbl_mb(i_rdrs),
     + dbl_mb(i_xdrs))
      endif
c
      if(lpmf) then
      call cf_plist(nums,xs,fs,isgan,isatt,int_mb(i_ipmf),
     + int_mb(i_jpmf),dbl_mb(i_rpmf),dbl_mb(i_xpmf),dbl_mb(i_wpmf),
     + lfnpmf,lrec)
      endif
c
      return
      end
      subroutine cf_plist(nums,xs,fs,isgan,isatt,
     + ipmf,jpmf,rpmf,xpmf,wpmf,lfnpmf,lrec)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 cf_atmass
      external cf_atmass
c
      integer nums,lfnpmf
      real*8 xs(msa,3),fs(msa,3)
      integer isgan(msa),isatt(msa)
      integer ipmf(numpmf,7),jpmf(numpmf,4,npmfa)
      real*8 rpmf(3,6,numpmf),xpmf(4,3,numpmf),wpmf(4,numpmf)
      real*8 et(2),ft(4,3)
      integer ip(5)
      real*8 rv(2),chg(2)
      integer it(2),ign(4),iml(4),idt(4),ich(4)
      logical lrec
c
      integer i,j,k,l,m
      real*8 wm,wmtot
      logical skip
c
      call timer_start(29)
      epmf=zero
      skip=.true.
c
      chg(1)=zero
c
      do 1 i=1,numpmf
      do 2 j=1,4
      xpmf(j,1,i)=zero
      xpmf(j,2,i)=zero
      xpmf(j,3,i)=zero
      wpmf(j,i)=zero
    2 continue
      if(ipmf(i,6).eq.2.and.lpmfc) goto 1
      m=0
c      if(ipmf(i,1).eq.1.or.ipmf(i,1).eq.5.or.ipmf(i,1).eq.6) m=2
      if(ipmf(i,1).eq.1) m=2
      if(ipmf(i,1).eq.2) m=3
      if(ipmf(i,1).eq.3) m=4
      if(ipmf(i,1).eq.4) m=4
      if(ipmf(i,1).eq.5) m=2
      if(ipmf(i,1).eq.6) m=1
      if(m.eq.0) call errquit('Error in pmf xx',ipmf(i,1))
      skip=.false.
      do 3 j=1,m
      wmtot=0.0d0
      do 4 k=1,ipmf(i,j+1)
      do 5 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      wm=cf_atmass(l,isatt)
      wmtot=wmtot+wm
      xpmf(j,1,i)=xpmf(j,1,i)+wm*xs(l,1)
      xpmf(j,2,i)=xpmf(j,2,i)+wm*xs(l,2)
      xpmf(j,3,i)=xpmf(j,3,i)+wm*xs(l,3)
      wpmf(j,i)=wpmf(j,i)+wm
      endif
    5 continue
    4 continue
    3 continue
    1 continue
      call timer_stop(29)
c
      if(skip) return
c
      call timer_start(30)
      call ga_dgop(mcf_59,wpmf,4*numpmf,'+')
      call ga_dgop(mcf_73,xpmf,12*numpmf,'+')
      call timer_stop(30)
c
      call timer_start(31)
      do 6 i=1,numpmf
      if(ipmf(i,6).eq.2.and.lpmfc) goto 6
c
      m=0
      if(ipmf(i,1).eq.1) m=2
      if(ipmf(i,1).eq.2) m=3
      if(ipmf(i,1).eq.3) m=4
      if(ipmf(i,1).eq.4) m=4
      if(ipmf(i,1).eq.5) m=2
      if(ipmf(i,1).eq.6) m=1
      do 222 j=1,m
      wmtot=one/wpmf(j,i)
      wpmf(j,i)=wmtot
      xpmf(j,1,i)=xpmf(j,1,i)*wmtot
      xpmf(j,2,i)=xpmf(j,2,i)*wmtot
      xpmf(j,3,i)=xpmf(j,3,i)*wmtot
  222 continue
      if(ipmf(i,1).eq.5) xpmf(2,3,i)=xpmf(1,3,i)
c
      if(ipmf(i,1).eq.1.or.ipmf(i,1).eq.5) then
      do 7 j=1,2
      ign(j)=j
      iml(j)=1
      idt(j)=0
      ich(j)=1
      do 8 k=1,3
      ft(j,k)=zero
    8 continue
    7 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=0
      call cf_fsb(1,it,1,3,ip,rpmf(1,1,i),rv,2,4,ign,iml,idt,ich,chg,
     + xpmf(1,1,i),ft,et,.false.,.false.,.false.,dbl_mb(i_sti))
      if(me.eq.0) then
      if(lrec.and.lfnpmf.lt.0) write(-lfnpmf,1000) i,rv
 1000 format(i3,f12.6,e12.5)
      if(ipmf(i,6).eq.1) ubias=ubias+rv(2)
      epmf=epmf+rv(2)
      endif
      do 9 j=1,2
      do 10 k=1,ipmf(i,j+1)
      do 11 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      wm=cf_atmass(l,isatt)
      fs(l,1)=fs(l,1)+ft(j,1)*wm*wpmf(j,i)
      fs(l,2)=fs(l,2)+ft(j,2)*wm*wpmf(j,i)
      fs(l,3)=fs(l,3)+ft(j,3)*wm*wpmf(j,i)
      endif
   11 continue
   10 continue
    9 continue
      endif
c
      if(ipmf(i,1).eq.2) then
      do 12 j=1,3
      ign(j)=j
      iml(j)=1
      idt(j)=0
      ich(j)=1
      do 13 k=1,3
      ft(j,k)=zero
   13 continue
   12 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=3
      ip(4)=0
      call cf_fsh(1,it,1,3,ip,rpmf(1,1,i),rv,3,4,ign,iml,idt,ich,chg,
     + xpmf(1,1,i),ft,et,.true.,.true.,.false.,dbl_mb(i_sti))
c      write(*,'(a,f12.6,e12.5,3e12.5)') 'Angle:',rv,(ft(3,j),j=1,3)
      if(me.eq.0) then
      if(lrec.and.lfnpmf.lt.0) write(-lfnpmf,1000) i,rv
      epmf=epmf+rv(2)
      if(ipmf(i,6).eq.1) ubias=ubias+rv(2)
      endif
      do 14 j=1,3
      do 15 k=1,ipmf(i,j+1)
      do 16 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      wm=cf_atmass(l,isatt)
      fs(l,1)=fs(l,1)+ft(j,1)*wm*wpmf(j,i)
      fs(l,2)=fs(l,2)+ft(j,2)*wm*wpmf(j,i)
      fs(l,3)=fs(l,3)+ft(j,3)*wm*wpmf(j,i)
      endif
   16 continue
   15 continue
   14 continue
      endif
c
      if(ipmf(i,1).eq.3.or.ipmf(i,1).eq.4) then
      do 17 j=1,4
      ign(j)=j
      iml(j)=1
      idt(j)=0
      ich(j)=1
      do 18 k=1,3
      ft(j,k)=zero
   18 continue
   17 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=3
      ip(4)=4
      ip(5)=0
      if(ipmf(i,1).eq.3) then
      call cf_fsd(1,it,1,3,ip,rpmf(1,1,i),rv,4,4,ign,iml,idt,
     + xpmf(1,1,i),ft,et,.true.,.true.,.false.,dbl_mb(i_sti))
c      write(*,'(a,f12.6,e12.5,3e12.5)') 'Dih:',rv,(ft(3,j),j=1,3)
      endif
      if(ipmf(i,1).eq.4) then
      call cf_fso(1,it,1,3,ip,rpmf(1,1,i),rv,4,4,ign,iml,idt,
     + xpmf(1,1,i),ft,et,.true.,.true.,.false.,dbl_mb(i_sti))
c      write(*,'(a,f12.6,e12.5,3e12.5)') 'Imp:',rv,(ft(3,j),j=1,3)
      endif
      if(me.eq.0) then
      if(lrec.and.lfnpmf.lt.0) write(-lfnpmf,1000) i,rv
      epmf=epmf+rv(2)
      if(ipmf(i,6).eq.1) ubias=ubias+rv(2)
      endif
      do 19 j=1,4
      do 20 k=1,ipmf(i,j+1)
      do 21 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      wm=cf_atmass(l,isatt)
      fs(l,1)=fs(l,1)+ft(j,1)*wm*wpmf(j,i)
      fs(l,2)=fs(l,2)+ft(j,2)*wm*wpmf(j,i)
      fs(l,3)=fs(l,3)+ft(j,3)*wm*wpmf(j,i)
      endif
   21 continue
   20 continue
   19 continue
      endif
c
c     planar pmf
c
      if(ipmf(i,1).eq.6) then
c
c     copy x and y coordinates of center of mass
c
      xpmf(2,1,i)=xpmf(1,1,i)
      xpmf(2,2,i)=xpmf(1,2,i)
c
c     loop over all atoms on this processor
c
      do 22 j=1,nums
c
c     loop over all atoms in the selected group
c
      do 23 k=1,ipmf(i,2)
      if(isgan(j).eq.jpmf(i,1,k)) then
c
c     copy z coordinate of the atom in xpmf
c
      xpmf(2,3,i)=xs(j,3)
c
c     initialize arrays for fsb
c
      do 24 l=1,2
      ign(l)=l
      iml(l)=1
      idt(l)=0
      ich(l)=1
      do 25 m=1,3
      ft(l,m)=zero
   25 continue
   24 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=0
      call cf_fsb(1,it,1,3,ip,rpmf(1,1,i),rv,2,4,ign,iml,idt,ich,chg,
     + xpmf(1,1,i),ft,et,.true.,.true.,.false.,dbl_mb(i_sti))
c
c     sum ubias on every processor since this is not done in replicated mode
c
      ubias=ubias+rv(2)
      epmf=epmf+rv(2)
c
c     sum the forces
c
      wm=cf_atmass(j,isatt)
      fs(j,3)=fs(j,3)+ft(2,3)*wm*wpmf(1,i)
c
c     note that the reaction force is not used since they sum to zero.
c
      endif
   23 continue
   22 continue
      endif
c
    6 continue
      call timer_stop(31)
c
      return
      end
      subroutine cf_dlist(nums,xs,fs,isgan,idrs,rdrs,xdrs)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nums
      real*8 xs(msa,3),fs(msa,3)
      integer isgan(msa)
      integer idrs(ndrs,2)
      real*8 rdrs(ndrs,6),xdrs(ndrs,6)
c
      integer i,j,k
      real*8 d,di,dd,dn,dx(3),dfor,df(3),en
c
      call timer_start(26)
      edrs=zero
c
      do 1 i=1,ndrs
      do 2 j=1,6
      xdrs(i,j)=zero
    2 continue
      do 3 j=1,nums
      if(isgan(j).eq.idrs(i,1)) then
      xdrs(i,1)=xs(j,1)
      xdrs(i,2)=xs(j,2)
      xdrs(i,3)=xs(j,3)
      endif
      if(isgan(j).eq.idrs(i,2)) then
      xdrs(i,4)=xs(j,1)
      xdrs(i,5)=xs(j,2)
      xdrs(i,6)=xs(j,3)
      endif
    3 continue
    1 continue
      call timer_stop(26)
c
      call timer_start(27)
      call ga_dgop(mcf_59,xdrs,6*ndrs,'+')
      call timer_stop(27)
c
      call timer_start(28)
      do 4 i=1,ndrs
      dx(1)=xdrs(i,1)-xdrs(i,4)
      dx(2)=xdrs(i,2)-xdrs(i,5)
      dx(3)=xdrs(i,3)-xdrs(i,6)
      d=sqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
      if(d.lt.tiny) then
      di=one
      else
      di=one/d
      endif
      dn=d-rdrs(i,2)
      k=0
      if(d.lt.rdrs(i,1)) then
      k=1
      dd=rdrs(i,1)-rdrs(i,2)
      en=rdrs(i,4)*(dn-half*dd)*dd
      dfor=rdrs(i,4)*di*dd
      elseif(d.gt.rdrs(i,3)) then
      k=4
      dd=rdrs(i,3)-rdrs(i,2)
      en=rdrs(i,5)*(dn-half*dd)*dd
      dfor=rdrs(i,5)*di*dd
      elseif(d.lt.rdrs(i,2)) then
      k=2
      en=half*rdrs(i,4)*dn*dn
      dfor=rdrs(i,4)*di*dn
      else
      k=3
      en=half*rdrs(i,5)*dn*dn
      dfor=rdrs(i,5)*di*dn
      endif
      edrs=edrs+en
      df(1)=dfor*dx(1)
      df(2)=dfor*dx(2)
      df(3)=dfor*dx(3)
c      if(me.eq.0) then
c      write(*,2222) k,d,(rdrs(i,j),j=1,3),en,dfor,df
c 2222 format(i3,4f12.6,e12.5,4f15.3)
c      endif
      do 5 j=1,nums
      if(isgan(j).eq.idrs(i,1)) then
      fs(j,1)=fs(j,1)-df(1)
      fs(j,2)=fs(j,2)-df(2)
      fs(j,3)=fs(j,3)-df(3)
      endif
      if(isgan(j).eq.idrs(i,2)) then
      fs(j,1)=fs(j,1)+df(1)
      fs(j,2)=fs(j,2)+df(2)
      fs(j,3)=fs(j,3)+df(3)
      endif
    5 continue
    4 continue
      call timer_stop(28)
c
      return
      end
      subroutine cf_mlist(nums,xs,fs,isgan,isfrc,isdyn,isct1,
     + msb,ibnd,msh,iang,msd,idih,mso,iimp,ixmul,imul,xmul,fmul,ess)
c
      implicit none
c
#include "cf_common.fh"
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
      call timer_start(23)
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
      call timer_stop(23)
c
      call timer_start(24)
      call ga_dgop(mcf_46,xmul,3*mmult,'+')
      call ga_igop(mcf_47,imul,4*mmult,'+')
      call timer_stop(24)
c
      call timer_start(25)
      if(me.eq.0) then
      if(nmult(1).gt.0) then
      ix=1
      call cf_fsb(nmult(1),ixmul(ix),mbt(2),mbp(2),
     + int_mb(i_ibnd(2)),dbl_mb(i_bnd(2)),dbl_mb(i_rbnd(2)),
     + nmul,mmult,imul(1,1),imul(1,2),imul(1,3),imul(1,4),dbl_mb(i_chg),
     + xmul,fmul,ess,.true.,.true.,me.eq.0,dbl_mb(i_sti))
      endif
c
      if(nmult(2).gt.0) then
      ix=1+nmult(1)
      call cf_fsh(nmult(2),ixmul(ix),mht(2),mhp(2),
     + int_mb(i_iang(2)),dbl_mb(i_ang(2)),dbl_mb(i_rang(2)),
     + nmul,mmult,imul(1,1),imul(1,2),imul(1,3),imul(1,4),dbl_mb(i_chg),
     + xmul,fmul,ess,.true.,.true.,me.eq.0,dbl_mb(i_sti))
      endif
c
      if(nmult(3).gt.0) then
      ix=1+nmult(1)+nmult(2)
      call cf_fsd(nmult(3),ixmul(ix),mdt(2),mdp(2),
     + int_mb(i_idih(2)),dbl_mb(i_dih(2)),dbl_mb(i_rdih(2)),
     + nmul,mmult,imul(1,1),imul(1,2),imul(1,3),
     + xmul,fmul,ess,.true.,.true.,me.eq.0,dbl_mb(i_sti))
      endif
c
      if(nmult(4).gt.0) then
      ix=1+nmult(1)+nmult(2)+nmult(3)
      call cf_fso(nmult(4),ixmul(ix),mit(2),mip(2),
     + int_mb(i_iimp(2)),dbl_mb(i_imp(2)),dbl_mb(i_rimp(2)),
     + nmul,mmult,imul(1,1),imul(1,2),imul(1,3),
     + xmul,fmul,ess,.true.,.true.,me.eq.0,dbl_mb(i_sti))
      endif
      endif
      call timer_stop(25)
c
      call timer_start(24)
      call ga_brdcst(mcf_74,fmul,3*mmult*ma_sizeof(mt_dbl,1,mt_byte),0)
      call timer_stop(24)
c
      call timer_start(25)
      do 3 i=1,nmul
      do 4 j=1,nums
      if(isgan(j).eq.imul(i,1)) then
      fs(j,1)=fs(j,1)+fmul(i,1)
      fs(j,2)=fs(j,2)+fmul(i,2)
      fs(j,3)=fs(j,3)+fmul(i,3)
      endif
    4 continue
    3 continue
      call timer_stop(25)
c
      return
      end
      subroutine cf_gauss(tgauss,frgaus,numw,nums,vw,vs,iwdt,isdt,isat)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      real*8 tgauss,frgaus
      integer numw,nums
      real*8 vw(mwm,3,mwa),vs(msa,3)
      integer iwdt(mwm),isdt(msa),isat(msa)
c
      call cf_gauss1(tgauss,frgaus,
     + numw,nums,vw,vs,iwdt,isdt,int_mb(i_iwa),
     + isat,dbl_mb(i_mas))
c
      return
      end
      subroutine cf_gauss1(tgauss,frgaus,
     + numw,nums,vw,vs,iwdt,isdt,iwat,isat,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "util.fh"
c
      real*8 tgauss,frgaus
      integer numw,nums
      real*8 vw(mwm,3,mwa),vs(msa,3),wgt(mat,6)
      integer iwdt(mwm),isdt(msa),iwat(mwa),isat(msa)
c
      real*8 facdof,rran,tdev
      integer iwa,iwm,ix,isa,iran
c
      tdev=zero
c
      if(numw.gt.0) then
      facdof=rgas*tgauss
      do 1 iwa=1,nwa
      if(facdof.gt.zero) tdev=sqrt(facdof/wgt(iwat(iwa),iset))
      do 2 iwm=1,numw
      if(iwdt(iwm).ne.0) then
      vw(iwm,1,iwa)=zero
      vw(iwm,2,iwa)=zero
      vw(iwm,3,iwa)=zero
      else
      do 3 ix=1,3
      rran=zero
      do 4 iran=1,12
      rran=rran+util_random(0)
    4 continue
      if(facdof.gt.zero) then
      vw(iwm,ix,iwa)=(one-frgaus)*vw(iwm,ix,iwa)+frgaus*(rran-six)*tdev
      else
      vw(iwm,ix,iwa)=zero
      endif
    3 continue
      endif
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      facdof=rgas*tgauss
      do 5 isa=1,nums
      if(isdt(isa).ne.0) then
      vs(isa,1)=zero
      vs(isa,2)=zero
      vs(isa,3)=zero
      else
      do 6 ix=1,3
      rran=zero
      do 7 iran=1,12
      rran=rran+util_random(0)
    7 continue
      if(facdof.gt.zero) then
      vs(isa,ix)=(one-frgaus)*vs(isa,ix)+
     + frgaus*(rran-six)*sqrt(facdof/wgt(isat(isa),iset))
      else
      vs(isa,ix)=zero
      endif
    6 continue
      endif
    5 continue
      endif
c
      return
      end
      subroutine cf_vscale(tgauss,numw,nums,vw,vs,iwdt,isdt,isat)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      real*8 tgauss
      integer numw,nums
      real*8 vw(mwm,3,mwa),vs(msa,3)
      integer iwdt(mwm),isdt(msa),isat(msa)
c
      call cf_vscal(tgauss,
     + numw,nums,vw,vs,iwdt,isdt,int_mb(i_iwa),
     + isat,dbl_mb(i_mas))
c
      return
      end
      subroutine cf_vscal(tgauss,
     + numw,nums,vw,vs,iwdt,isdt,iwat,isat,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "util.fh"
#include "msgids.fh"
c
      real*8 tgauss
      integer numw,nums
      real*8 vw(mwm,3,mwa),vs(msa,3),wgt(mat,6)
      integer iwdt(mwm),isdt(msa),iwat(mwa),isat(msa)
c
      real*8 ekwcg,ekscg,tempgw,tempgs,gaufcw,gaufcs,ekcg(2)
      integer iwm,iwa,ix,isa
c
      ekwcg=zero
      ekscg=zero
c
      if(numw.gt.0) then
      do 1 iwm=1,numw
      if(iwdt(iwm).ne.0) then
      do 2 iwa=1,nwa
      do 3 ix=1,3
      vw(iwm,ix,iwa)=zero
    3 continue
    2 continue
      else
      do 4 iwa=1,nwa
      ekwcg=ekwcg+wgt(iwat(iwa),iset)*(vw(iwm,1,iwa)*vw(iwm,1,iwa)
     + +vw(iwm,2,iwa)*vw(iwm,2,iwa)+vw(iwm,3,iwa)*vw(iwm,3,iwa))
    4 continue
      endif
    1 continue
      endif
c
      if(nums.gt.0) then
      do 5 isa=1,nums
      if(isdt(isa).ne.0) then
      do 6 ix=1,3
      vs(isa,ix)=zero
    6 continue
      else
      ekscg=ekscg+wgt(isat(isa),iset)*(vs(isa,1)*vs(isa,1)+
     + vs(isa,2)*vs(isa,2)+vs(isa,3)*vs(isa,3))
      endif
    5 continue
      endif
c
c     globally sum kinetic energies in ekcg
c
      ekcg(1)=half*ekwcg
      ekcg(2)=half*ekscg
      call ga_dgop(mcf_48,ekcg,2,'+')
      ekwcg=ekcg(1)
      ekscg=ekcg(2)
c
c      if(ivreas.ge.2) then
c      tempgw=fcwtmp*ekwcg
c      tempgs=fcstmp*ekscg
c      else
      tempgw=factmp*(ekwcg+ekscg)
      tempgs=tempgw
c      endif
c
      if(numw.gt.0.and.tempgw.gt.tiny) then
      gaufcw=sqrt(tgauss/tempgw)
      do 7 iwa=1,nwa
      do 8 ix=1,3
      do 9 iwm=1,numw
      vw(iwm,ix,iwa)=gaufcw*vw(iwm,ix,iwa)
    9 continue
    8 continue
    7 continue
      endif
c
      if(nums.gt.0.and.tempgs.gt.tiny) then
      gaufcs=sqrt(tgauss/tempgs)
      do 10 ix=1,3
      do 11 isa=1,nums
      vs(isa,ix)=gaufcs*vs(isa,ix)
   11 continue
   10 continue
      endif
c
      return
      end
      subroutine cf_mcti(fdata)
c
      implicit none
c
#include "cf_common.fh"
c
      real*8 fdata(26)
c
      integer i
c
      do 1 i=1,24
      fdata(i)=deriv(i,1)
    1 continue
      fdata(25)=ep2(1)
      fdata(26)=ep3(1)
c
      return
      end
      subroutine cf_mcti_kin(isat,vs,nums)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer nums
      integer isat(msa)
      real*8 vs(msa,3)
c
      call cf_mcti_w(isat,vs,nums,dbl_mb(i_mas))
c
      return
      end
      subroutine cf_mcti_w(isat,vs,nums,wgt)
c
      implicit none
c
#include "cf_common.fh"
c
      integer nums
      integer isat(msa)
      real*8 vs(msa,3),wgt(mat,6)
c
      integer i
c
      do 1 i=1,nums
      deriv(13,1)=deriv(13,1)+wgt(isat(i),4)*
     + (vs(i,1)*vs(i,1)+vs(i,2)*vs(i,2)+vs(i,3)*vs(i,3))
    1 continue
c
      return
      end
      subroutine cf_center(xw,numw,isfrc,xs,nums,idscb,nscb,icentr)
c
      implicit none
c
#include "cf_common.fh"
#include "msgids.fh"
c
      integer numw,nums,nscb,icentr
      integer isfrc(msa),idscb(nscb)
      real*8 xw(mwm,3,mwa),xs(msa,3)
c
      real*8 xb(4)
      integer i,isa,iwm,iwa,j,n,i0,i9
c
      if(nsa.le.0) return
c
      i0=1
      i9=3
      if(icentr.eq.1) i0=3
      if(icentr.eq.2) i9=2
c
      do 1 i=1,3
      xb(i)=zero
    1 continue
c
      n=0
      do 2 isa=1,nums
      do 3 j=1,nscb
      if(idscb(j).eq.isfrc(isa)) then
      n=n+1
      do 4 i=i0,i9
      xb(i)=xb(i)+xs(isa,i)
    4 continue
      endif
    3 continue
    2 continue
c
      xb(4)=dble(n)
c
      call ga_dgop(mcf_49,xb,4,'+')
c
      do 5 i=i0,i9
      xb(i)=xb(i)/xb(4)
    5 continue
c
      do 6 i=i0,i9
      do 7 isa=1,nums
      xs(isa,i)=xs(isa,i)-xb(i)
    7 continue
    6 continue
c
      do 8 i=i0,i9
      do 9 iwa=1,nwa
      do 10 iwm=1,numw
      xw(iwm,i,iwa)=xw(iwm,i,iwa)-xb(i)
   10 continue
    9 continue
    8 continue
c
      return
      end
      subroutine cf_slow(vw,numw,vs,isat,nums)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw,nums
      integer isat(msa)
      real*8 vw(mwm,3,mwa),vs(msa,3)
c
      call cf_slowt(vw,numw,vs,nums,int_mb(i_iwa),isat,dbl_mb(i_mas))
c
      return
      end
      subroutine cf_slowt(vw,numw,vs,nums,iwat,isat,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "msgids.fh"
c
      integer numw,nums
      real*8 vw(mwm,3,mwa),vs(msa,3),wgt(mat,6)
      integer iwat(mwa),isat(msa)
c
      integer i,j,k
      real*8 vb(3)
c
      do 1 i=1,3
      vb(i)=zero
    1 continue
c
      do 2 k=1,3
      do 3 j=1,nwa
      do 4 i=1,numw
      vb(k)=vb(k)+wgt(iwat(j),iset)*vw(i,k,j)
    4 continue
    3 continue
      do 5 i=1,nums
      vb(k)=vb(k)+wgt(isat(j),iset)*vs(i,k)
    5 continue
      vb(k)=vb(k)/wbox
    2 continue
c
      call ga_dgop(mcf_50,vb,3,'+')
c
      do 6 k=1,3
      do 7 j=1,nwa
      do 8 i=1,numw
      vw(i,k,j)=vw(i,k,j)-vb(k)
    8 continue
    7 continue
      do 9 i=1,nums
      vs(i,k)=vs(i,k)-vb(k)
    9 continue
    6 continue
c
      return
      end
      subroutine cf_print_energy(lfnout)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfnout
c
      if(numb(2).gt.0) call cf_pr_bonds(lfnout,numb(2),mbt(2),mbp(2),
     + int_mb(i_ibnd(2)),dbl_mb(i_bnd(2)),dbl_mb(i_rbnd(2)),
     + byte_mb(i_snam))
c
      if(numh(2).gt.0) call cf_pr_angles(lfnout,numh(2),mht(2),mhp(2),
     + int_mb(i_iang(2)),dbl_mb(i_ang(2)),dbl_mb(i_rang(2)))
c
      if(numd(2).gt.0) call cf_pr_dihes(lfnout,numd(2),mdt(2),mdp(2),
     + int_mb(i_idih(2)),dbl_mb(i_dih(2)),dbl_mb(i_rdih(2)))
c
      if(numi(2).gt.0) call cf_pr_imprs(lfnout,numi(2),mit(2),mip(2),
     + int_mb(i_iimp(2)),dbl_mb(i_imp(2)),dbl_mb(i_rimp(2)))
c
      return
      end
      subroutine cf_pr_bonds(lfnout,nsb,msb,msp,ibnd,bnd,rbnd,snam)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfnout,nsb,msb,msp
      integer ibnd(msb,3)
      real*8 bnd(msb,msp,6),rbnd(msb,2)
      character*16 snam(nsatot)
c
      integer i,j
c
      if(np.gt.1) call ga_dgop(mcf_51,rbnd,2*msb,'+')
      if(me.ne.0) return
c
      write(lfnout,1000)
 1000 format(/,' Solute bonds',/,
     + '      i      j      c     b_0        C          b          U',/)
c
      do 1 i=1,nsb
      if(ibnd(i,3).eq.0) then
      write(lfnout,1001) (ibnd(i,j),j=1,2),
     + (snam(ibnd(i,j))(12:16)//':'//snam(ibnd(i,j))(1:11),j=1,2),
     + (bnd(i,j,iset),j=1,2),(rbnd(i,j),j=1,2)
 1001 format(2i7,2a17,7x,f10.5,e12.5,f10.5,f12.6)
      else
      write(lfnout,1002) (ibnd(i,j),j=1,2),
     + (snam(ibnd(i,j))(12:16)//':'//snam(ibnd(i,j))(1:11),j=1,2),
     + (bnd(i,j,iset),j=1,2),rbnd(i,1)
 1002 format(2i7,2a17,' cnstr ',f10.5,e12.5,f10.5,f12.6)
      endif
    1 continue
c
      return
      end
      subroutine cf_pr_angles(lfnout,nsh,msh,msp,iang,ang,rang)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfnout,nsh,msh,msp
      integer iang(msh,4)
      real*8 ang(msh,msp,6),rang(msh,2)
c
      integer i,j
c
      if(np.gt.1) call ga_dgop(mcf_52,rang,2*msh,'+')
      if(me.ne.0) return
c
      write(lfnout,1000)
 1000 format(/,' Solute angles',/,
     + '      i      j      k      c    th_0   ',
     + '     C         th          U',/)
c
      do 1 i=1,nsh
      write(lfnout,1001) (iang(i,j),j=1,4),(ang(i,j,iset),j=1,2),
     + (rang(i,j),j=1,2)
 1001 format(4i7,f10.5,e12.5,f10.5,f12.6)
    1 continue
c
      return
      end
      subroutine cf_pr_dihes(lfnout,nsd,msd,msp,idih,dih,rdih)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfnout,nsd,msd,msp
      integer idih(msd,5)
      real*8 dih(msd,msp,6),rdih(msd,2)
c
      integer i,j
c
      if(np.gt.1) call ga_dgop(mcf_53,rdih,2*msd,'+')
      if(me.ne.0) return
c
      write(lfnout,1000)
 1000 format(/,' Solute torsions',/,
     + '      i      j      k      l      c  mul   phi_0   ',
     + '     C        phi          U',/)
c
      do 1 i=1,nsd
      write(lfnout,1001) (idih(i,j),j=1,5),(dih(i,j,iset),j=1,3),
     + (rdih(i,j),j=1,2)
 1001 format(5i7,f5.1,f10.5,e12.5,f10.5,f12.6)
    1 continue
c
      return
      end
      subroutine cf_pr_imprs(lfnout,nso,mso,msp,ioop,oop,roop)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfnout,nso,mso,msp
      integer ioop(mso,5)
      real*8 oop(mso,msp,6),roop(mso,2)
c
      integer i,j
c
      if(np.gt.1) call ga_dgop(mcf_54,roop,2*mso,'+')
      if(me.ne.0) return
c
      write(lfnout,1000)
 1000 format(/,' Solute out of plane torsions',/,
     + '      i      j      k      l      c  mul theta_0   ',
     + '     C      theta          U',/)
c
      do 1 i=1,nso
      write(lfnout,1001) (ioop(i,j),j=1,5),(oop(i,j,iset),j=1,3),
     + (roop(i,j),j=1,2)
 1001 format(5i7,f5.1,f10.5,e12.5,f10.5,f12.6)
    1 continue
c
      return
      end
      subroutine cf_fix(lfnout,lfn,numfix,iwgmn,iwdt,numw,
     + isgan,isatt,isdt,issgm,nums,xwm,xs)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfnout,lfn,numfix,numw,nums
      integer iwdt(mwm),iwgmn(mwm),isgan(msa),isatt(msa),isdt(msa)
      integer issgm(msa)
      real*8 xwm(mwm,3),xs(msa,3)
c
      integer i_temp,l_temp
c
      if(.not.ma_push_get(mt_int,nsatot,'itemp',l_temp,i_temp))
     + call errquit('Unable to allocate temp array in cf_fix',0) 
c
      call cf_fix_sw(lfnout,lfn,numfix,iwgmn,iwdt,numw,
     + isgan,isatt,isdt,issgm,nums,byte_mb(i_snam),xwm,xs,
     + int_mb(i_ibnd(2)),mbt(2),int_mb(i_temp))
c
      if(.not.ma_pop_stack(l_temp)) 
     + call errquit('Unable to deallocate temp array in cf_fix',0) 
      return
      end
      subroutine cf_fix_sw(lfnout,lfn,numfix,iwgmn,iwdt,numw,
     + isgan,isatt,isdt,issgm,nums,snam,xwm,xs,idb,nsb,itemp)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer cf_atonum,cf_numcon
      external cf_atonum,cf_numcon
c
      integer lfnout,lfn,numfix,numw,nums,nsb
      integer iwdt(mwm),iwgmn(mwm),isgan(msa),isatt(msa),isdt(msa)
      integer issgm(msa),idb(nsb,4),itemp(nsatot)
      character*16 snam(nsatot)
      real*8 xwm(mwm,3),xs(msa,3)
c
      integer i,j,ndof(3),ifix,ifr,ito,l,len,number,itarget,ntotal(2)
      integer ib,jb
      real*8 fact(3),xtarget(3),rfix,dist2,rfix2
      character*80 string,target,residu
      character*16 tag
c
      if(me.eq.0) then
      write(lfnout,1000) numfix
 1000 format(/,' POSITIONAL CONSTRAINTS',/,
     + ' Number of fix directives is ',i5,/)
      endif
c
      number=0
c
      do 1 l=1,numfix
c
      if(me.eq.0) then
      read(lfn,2000) string
 2000 format(a80)
c
      target=string(41:80)
      len=index(target,':')
      residu=' '
      if(len.gt.0) then
      residu=target(1:len-1)
      target=target(len+1:80)
      target(13-len:11)=residu
      endif
      do 33 i=1,4
      if(target(i:i).eq.'_') target(i:i)=' '
   33 continue
c
      itarget=0
      do 34 i=1,nsatot
      if(snam(i)(6:16).eq.target(1:11)) itarget=i
   34 continue
c
      endif
c
      call ga_brdcst(mrg_d47,itarget,ma_sizeof(mt_int,1,mt_byte),0)
      call util_char_ga_brdcst(mrg_d47,string,0)
c
      if(itarget.gt.0) then
      do 35 i=1,nums
      if(isgan(i).eq.itarget) then
      itarget=i
      xtarget(1)=xs(i,1)
      xtarget(2)=xs(i,2)
      xtarget(3)=xs(i,3)
      goto 36
      endif
   35 continue
      xtarget(1)=zero
      xtarget(2)=zero
      xtarget(3)=zero
      itarget=0
   36 continue
      call ga_dgop(mrg_d47,xtarget,3,'+')
      endif
c
      if(string(1:3).eq.'fix') then
c
      read(string(6:40),'(i3,2i10,f12.6)') ifix,ifr,ito,rfix
      target=string(41:80)
      len=index(target,':')
      residu=' '
      if(len.gt.0) then
      residu=target(1:len-1)
      target=target(len+1:80)
      endif
      len=4
      do 21 i=1,4
      if(target(i:i).eq.'_') target(i:i)=' '
   21 continue
c
c     solvent
c
      if(string(5:5).eq.'w') then
      if(numw.gt.0) then
      if(ifr.eq.0.or.ito.eq.0) then
      if(ifix.eq.0) then
      do 2 i=1,numw
      iwdt(i)=iand(iwdt(i),not(mfixed))
    2 continue
      endif
      if(ifix.eq.1) then
      do 3 i=1,numw
      iwdt(i)=ior(iwdt(i),lfixed)
    3 continue
      endif
      else
      if(ifr.gt.0) then
      do 4 i=1,numw
      if(iwgmn(i).ge.ifr.and.iwgmn(i).le.ito) then
      if(ifix.eq.0) iwdt(i)=iand(iwdt(i),not(mfixed))
      if(ifix.eq.1) iwdt(i)=ior(iwdt(i),lfixed)
      endif
    4 continue
      else
      rfix2=rfix*rfix
      do 41 i=1,numw
      dist2=(xwm(i,1)-xtarget(1))**2+(xwm(i,2)-xtarget(2))**2+
     + (xwm(i,3)-xtarget(3))**2
      if(ito.eq.-1.and.dist2.le.rfix2) then
      if(ifix.eq.0) iwdt(i)=iand(iwdt(i),not(mfixed))
      if(ifix.eq.1) iwdt(i)=ior(iwdt(i),lfixed)
      elseif(ito.eq.-2.and.dist2.ge.rfix2) then
      if(ifix.eq.0) iwdt(i)=iand(iwdt(i),not(mfixed))
      if(ifix.eq.1) iwdt(i)=ior(iwdt(i),lfixed)
      endif
   41 continue
      endif
      endif
      endif
      endif
c
c     solute
c
      if(string(5:5).eq.'s') then
      if(nums.gt.0) then
      if(ifr.eq.0.or.ito.eq.0) then
      if(ifix.eq.0) then
      do 5 i=1,nums
      tag=snam(isgan(i))
      do 51 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   51 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      isdt(i)=iand(isdt(i),not(mfixed))
      number=number+1
      endif
    5 continue
      endif
      if(ifix.eq.1) then
      do 6 i=1,nums
      tag=snam(isgan(i))
      do 61 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   61 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      isdt(i)=ior(isdt(i),lfixed)
      number=number+1
      endif
    6 continue
      endif
      else
      if(ifr.gt.0) then
      do 7 i=1,nums
      if(issgm(i).ge.ifr.and.issgm(i).le.ito) then
      tag=snam(isgan(i))
      do 71 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   71 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      number=number+1
      if(ifix.eq.0) isdt(i)=ior(isdt(i),not(mfixed))
      if(ifix.eq.1) isdt(i)=ior(isdt(i),lfixed)
      endif
      endif
    7 continue
      else
      rfix2=rfix*rfix
      do 72 i=1,nums
      dist2=(xs(i,1)-xtarget(1))**2+(xs(i,2)-xtarget(2))**2+
     + (xs(i,3)-xtarget(3))**2
      if(ito.eq.-1.and.dist2.le.rfix2) then
      if(ifix.eq.0) isdt(i)=iand(isdt(i),not(mfixed))
      if(ifix.eq.1) isdt(i)=ior(isdt(i),lfixed)
      elseif(ito.eq.-2.and.dist2.ge.rfix2) then
      if(ifix.eq.0) isdt(i)=iand(isdt(i),not(mfixed))
      if(ifix.eq.1) isdt(i)=ior(isdt(i),lfixed)
      endif
   72 continue
      endif
      endif
      endif
      endif
c
c     heavy solute atoms only
c
      if(string(5:5).eq.'x') then
      if(nums.gt.0) then
      if(ifr.eq.0.or.ito.eq.0) then
      if(ifix.eq.0) then
      do 8 i=1,nums
      tag=snam(isgan(i))
      do 81 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   81 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      if(cf_atonum(isatt(i),int_mb(i_num)).ne.1)
     + isdt(i)=iand(isdt(i),not(mfixed))
      number=number+1
      endif
    8 continue
      endif
      if(ifix.eq.1) then
      do 9 i=1,nums
      tag=snam(isgan(i))
      do 91 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   91 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      if(cf_atonum(isatt(i),int_mb(i_num)).ne.1)
     + isdt(i)=ior(isdt(i),lfixed)
      number=number+1
      endif
    9 continue
      endif
      else
      do 10 i=1,nums
      if(issgm(i).ge.ifr.and.issgm(i).le.ito) then
      tag=snam(isgan(i))
      do 101 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
  101 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      if(cf_atonum(isatt(i),int_mb(i_num)).ne.1) then
      if(ifix.eq.0) isdt(i)=ior(isdt(i),not(mfixed))
      if(ifix.eq.1) isdt(i)=ior(isdt(i),lfixed)
      endif
      number=number+1
      endif
      endif
   10 continue
      endif
      endif
      endif
c
c     solute atoms constraints
c
      do 76 i=1,nsatot
      itemp(i)=0
   76 continue
c
      do 62 i=1,nums
      itemp(isgan(i))=i
   62 continue
c
      do 64 j=1,4
      do 63 i=1,nsb
      ib=itemp(idb(i,1))
      jb=itemp(idb(i,2))
      if(ib.gt.0.and.jb.gt.0.and.idb(i,3).gt.0) then
      if(iand(isdt(ib),mfixed).eq.lfixed) isdt(jb)=ior(isdt(jb),lfixed)
      if(iand(isdt(jb),mfixed).eq.lfixed) isdt(ib)=ior(isdt(ib),lfixed)
      endif
   63 continue
   64 continue
c
c     make permanent
c
      if(string(5:5).eq.'p') then
      do 77 i=1,nums
      j=iand(isdt(i),mfixed)/4
      isdt(i)=lshift(rshift(isdt(i),2),2)
      isdt(i)=ior(isdt(i),j)
   77 continue
      do 78 i=1,numw
      j=iand(iwdt(i),mfixed)/4
      iwdt(i)=lshift(rshift(iwdt(i),2),2)
      iwdt(i)=ior(iwdt(i),j)
   78 continue
      endif
c
      endif
c
    1 continue
c
      ntotal(1)=0
      ntotal(2)=0
      do 201 i=1,nums
      if(iand(isdt(i),mfixed).eq.lfixed) ntotal(1)=ntotal(1)+1
  201 continue
      do 202 i=1,numw
      if(iand(iwdt(i),mfixed).eq.lfixed) ntotal(2)=ntotal(2)+1
  202 continue
c
      if(np.gt.1) then
      call ga_igop(mrg_d48,ntotal,2,'+')
      endif
c
      if(me.eq.0) then
      write(lfnout,'(a,i8,a,i8,a)')
     +  ' Solute fix applied to  ',ntotal(1),' out of ',nsa,' atoms'
      write(lfnout,'(a,i8,a,i8,a,/)')
     +  ' Solvent fix applied to ',ntotal(2),' out of ',nwm,' molecules'
      endif
c
c      if(numw.gt.0) then
c      if(iand(ifix,1).eq.1) then
c      do 1 i=1,numw
c      iwdt(i)=ior(iwdt(i),lfixed)
c    1 continue
c      endif
c      if(iand(iunfix,1).eq.1) then
c      do 2 i=1,numw
c      iwdt(i)=iand(iwdt(i),not(lfixed))
c    2 continue
c      endif
c      if(iand(ifix,2).eq.2) then
c      do 3 i=1,numw
c      iwdt(i)=ior(iwdt(i),lfixed)
c    3 continue
c      endif
c      if(iand(iunfix,2).eq.2) then
c      do 4 i=1,numw
c      iwdt(i)=iand(iwdt(i),not(lfixed))
c    4 continue
c      endif
c      endif
c
c      if(nums.gt.0) then
c      if(iand(ifix,1).eq.1) then
c      do 5 i=1,nums
c      isdt(i)=ior(isdt(i),lfixed)
c    5 continue
c      endif
c      if(iand(iunfix,1).eq.1) then
c      do 6 i=1,nums
c      isdt(i)=iand(isdt(i),not(lfixed))
c    6 continue
c      endif
c      if(iand(ifix,4).eq.4) then
c      do 7 i=1,nums
c      isdt(i)=ior(isdt(i),lfixed)
c    7 continue
c      endif
c      if(iand(iunfix,4).eq.4) then
c      do 8 i=1,nums
c      isdt(i)=iand(isdt(i),not(lfixed))
c    8 continue
c      endif
c      if(iand(ifix,8).eq.8) then
c      do 9 i=1,nums
c      if(cf_atonum(isatt(i),int_mb(i_num)).ne.1)
c     + isdt(i)=ior(isdt(i),lfixed)
c    9 continue
c      endif
c      if(iand(iunfix,8).eq.8) then
c      do 10 i=1,nums
c      if(cf_atonum(isatt(i),int_mb(i_num)).ne.1)
c     + isdt(i)=iand(isdt(i),not(lfixed))
c   10 continue
c      endif      
c      endif
c
      ndof(1)=0
      ndof(2)=0
      ndof(3)=0
c
      do 11 i=1,numw
      if(iand(iwdt(i),lfixed).ne.lfixed) ndof(1)=ndof(1)+1
   11 continue
      do 12 i=1,nums
      if(iand(isdt(i),lfixed).ne.lfixed) then
      ndof(2)=ndof(2)+1
      ndof(3)=ndof(3)+cf_numcon(int_mb(i_ibnd(2)),mbt(2),isgan(i))
      endif
   12 continue
c
      call ga_igop(mrg_d46,ndof,3,'+')
      nsc=ndof(3)/2
c
      factmw=zero
      factms=zero
      factmp=zero
      if(ndof(1)*(3*nwa-nwc)-3*islow.gt.0)
     + factmw=two/(rgas*dble(ndof(1)*(3*nwa-nwc)-3*islow))
      if(3*ndof(2)-nsc-3*islow.gt.0)
     + factms=two/(rgas*dble(3*ndof(2)-nsc-3*islow))
      if(ndof(1)*(3*nwa-nwc)+3*ndof(2)-nsc-3*islow.gt.0)
     + factmp=two/(rgas*dble(ndof(1)*(3*nwa-nwc)+3*ndof(2)-nsc-3*islow))
c
      fact(1)=factmw
      fact(2)=factms
      fact(3)=factmp
      call ga_brdcst(mrg_d47,fact,3*ma_sizeof(mt_int,1,mt_byte),0)
      factmw=fact(1)
      factms=fact(2)
      factmp=fact(3)
c
      return
      end
      integer function cf_numcon(idb,ndb,ia)
c
      implicit none
c
      integer ia,ndb
      integer idb(ndb,4)
c
      integer i,number
c
      number=0
c
      do 1 i=1,ndb
      if(idb(i,3).eq.1.and.(idb(i,1).eq.ia.or.idb(i,2).eq.ia))
     + number=number+1
    1 continue
c
      cf_numcon=number
c
      return
      end
      real*8 function cf_atmass(i,isat)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      real*8 cf_am
      external cf_am
c
      integer i,isat(msa)
c
      cf_atmass=cf_am(isat(i),dbl_mb(i_mas))
c
      return
      end
      real*8 function cf_am(ia,wgt)
c
      implicit none
c
#include "cf_common.fh"
c
      integer ia
      real*8 wgt(mat,6)
c
      cf_am=wgt(ia,iset)
c
      return
      end
      integer function cf_atonum(ia,num)
c
      implicit none
c
#include "cf_common.fh"
c
      integer ia,num(mat,3)
c
      cf_atonum=num(ia,iset)
c
      return
      end
      integer function cf_atotyp(ityp,ia)
c
      implicit none
c
#include "cf_common.fh"
c
      integer ia
      integer ityp(mat,3)
c
      cf_atotyp=ityp(ia,iset)
c
      return
      end
      integer function cf_element(i)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer cf_atonum
      external cf_atonum
c
      integer i
c
      cf_element=cf_atonum(i,int_mb(i_num))
c
      return
      end
      integer function cf_quantuma(i)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer cf_atotyp
      external cf_atotyp
c
      integer i
c
      cf_quantuma=cf_atotyp(int_mb(i_typ),i)
c
      return
      end
      subroutine cf_links(ltype,nums,isgan,isdt,lst)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer nums,ltype
      integer isgan(msa),isdt(msa),lst(nsa)
c
      call cf_cntlnk(ltype,mbt(2),int_mb(i_ibnd(2)),nums,isgan,isdt,lst)
      return
      end
      subroutine cf_cntlnk(ltype,msb,ibnd,nums,isgan,isdt,lst)
c
      implicit none
c
#include "cf_common.fh"
c
      integer msb,nums,ltype
      integer ibnd(msb,3),isgan(msa),isdt(msa),lst(nsa)
c
      integer i,j
c
      do 1 j=1,nums
      if(iand(isdt(j),ltype).ne.0) then
      do 2 i=1,msb
      if(ibnd(i,1).eq.isgan(j)) lst(ibnd(i,2))=lst(ibnd(i,2))+1
      if(ibnd(i,2).eq.isgan(j)) lst(ibnd(i,1))=lst(ibnd(i,1))+1
    2 continue
      endif
    1 continue
c
      return
      end
      real*8 function cf_charge(iq)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      real*8 cf_chg
      external cf_chg
c
      integer iq
c
      cf_charge=cf_chg(iq,dbl_mb(i_chg))
      return
      end
      real*8 function cf_chg(iq,chg)
c
      implicit none
c
#include "cf_common.fh"
c
      integer iq
      real*8 chg(mqt,mqp,6)
c
      cf_chg=chg(iq,1,iset)/qfac
c
      return
      end
      real*8 function cf_wcharge(iq)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      real*8 cf_chg
      external cf_chg
c
      integer iq
c
      cf_wcharge=cf_chg(int_mb(i_iwq+iq-1),dbl_mb(i_chg))
      return
      end
      subroutine cf_num2nam(i,a)
c
      implicit none
c
      integer i
      character*16 a
c
      character*2 symbol(103)
c
      data symbol/
     + 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     + 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     + 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     + 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     + 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     + 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     + 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     + 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     + 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     + 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     + 'Md','No','Lr'/
c
      a='     Bq         '
      if(i.gt.0.and.i.lt.104) a(6:7)=symbol(i)
      if(a(7:7).eq.' ') a(6:7)=' '//a(6:6)
c
      return
      end
      subroutine cf_num2tag(i,a)
c
      implicit none
c
      integer i
      character*16 a
c
      character*2 symbol(103)
c
      data symbol/
     + 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     + 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     + 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     + 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     + 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     + 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     + 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     + 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     + 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     + 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     + 'Md','No','Lr'/
c
      a='Bq              '
      if(i.gt.0.and.i.lt.104) a(1:2)=symbol(i)
c
      return
      end
      subroutine cf_trjhdr(iunit)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer iunit
c
      call cf_wthdr(iunit,byte_mb(i_snam),
     + int_mb(i_ibnd(1)),mbt(1),int_mb(i_ibnd(2)),mbt(2))
c
      return
      end
      subroutine cf_wthdr(iunit,snam,idbw,nwb,idbs,nsb)
c
      implicit none
c
#include "cf_common.fh"
c
      integer iunit,nwb,nsb
      character*16 snam(nsatot)
      integer idbw(nwb,4),idbs(nsb,4)
c
      integer i
c
      rewind(iunit)
c
      write(iunit,1000)
 1000 format('header')
      write(iunit,1001) nwa,nsatot,numb(2),numb(1),nwmtot
 1001 format(5i10)
      if(nwa.gt.0) then
      do 1 i=1,nwa
      write(iunit,1004) wnam(i)(1:5),wnam(i)(11:15)
 1004 format(2a5)
    1 continue
      endif
      if(nsatot.gt.0) then
      do 2 i=1,nsatot
      write(iunit,1005) snam(i),i,1
 1005 format(a16,i10,i5)
    2 continue
      endif
      if(numb(1).gt.0) then
      do 3 i=1,numb(1)
      write(iunit,1006) idbw(i,1),idbw(i,2)
 1006 format(2i8)
    3 continue
      endif
      if(numb(2).gt.0) then
      do 4 i=1,numb(2)
      write(iunit,1006) idbs(i,1),idbs(i,2)
    4 continue
      endif
c
      return
      end
      subroutine cf_dera_init()
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      call cf_dera_i(dbl_mb(i_sti))
c
      return
      end
      subroutine cf_dera_i(dera)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      real*8 dera(5,nsatot)
c
      integer i
c
      do 1 i=1,nsatot
      dera(1,i)=zero
      dera(2,i)=zero
      dera(3,i)=zero
      dera(4,i)=zero
      dera(5,i)=zero
    1 continue
c
      return
      end
      subroutine cf_print_distar(lfnout)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfnout
c
      call cf_print_dstr(lfnout,byte_mb(i_snam),
     + int_mb(i_idrs),dbl_mb(i_rdrs))
c
      return
      end
      subroutine cf_print_dstr(lfnout,snam,idrs,rdrs)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfnout
      character*16 snam(nsatot)
      integer idrs(ndrs,2)
      real*8 rdrs(ndrs,6)
c
      integer i
c
      if(me.eq.0.and.ndrs.gt.0) then
      write(lfnout,1001)
 1001 format(' Distance restraints',/)
      do 1 i=1,ndrs
      write(lfnout,1002) snam(idrs(i,1)),snam(idrs(i,2)),
     + rdrs(i,1),rdrs(i,2),rdrs(i,3),rdrs(i,4),rdrs(i,5)
 1002 format(1x,2a16,3f12.6,2e12.5)
    1 continue
      endif
c
      return
      end
      subroutine cf_print_dera(lfnout,number)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfnout,number
c
      call cf_print_der(lfnout,number,byte_mb(i_snam),dbl_mb(i_sti))
c
      return
      end
      subroutine cf_print_der(lfnout,number,snam,dera)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfnout,number
      real*8 dera(5,nsatot)
      character*16 snam(nsatot)
c
      integer i,j
c
      call ga_dgop(mcf_68,dera,5*nsatot,'+')
c
      if(me.eq.0) then
c
      write(lfnout,1000) number
 1000 format(//,' Solute MCTI derivatives from last',i8,' steps',//,
     + '     Atom',13x,'Solvent     Solvent      Solute      Solute',
     + 6x,'Solute',/,
     + 21x,'v.d.Waals   electrost   v.d.Waals   electrost',6x,
     + 'bonded',/)
c
      do 1 i=1,nsatot
      write(lfnout,1001) snam(i),(dera(j,i)/dble(number),j=1,5)
 1001 format(1x,a16,5f12.3)
    1 continue
c
      endif
c
      return
      end
      subroutine cf_wrtmro(lfnmro)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfnmro
c
      call cf_wtmro(lfnmro,dbl_mb(i_sti))
c
      return
      end
      subroutine cf_wtmro(lfnmro,dera)
c
      implicit none
c
#include "cf_common.fh"
#include "msgids.fh"
c
      integer lfnmro
      real*8 dera(5,nsatot)
c
      integer i,j
c
      call ga_dgop(mcf_68,dera,5*nsatot,'+')
c
      if(me.eq.0) then
      write(lfnmro) ((dera(i,j),i=1,5),j=1,nsatot)
      else
      do 1 i=1,nsatot
      do 2 j=1,5
      dera(j,i)=zero
    2 continue
    1 continue
      endif
c
      return
      end
      subroutine cf_rdmri(lfnmri)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfnmri
c
      if(npgdec.eq.0) call errquit('rdmri error',0)
c
      call cf_rmri(lfnmri,dbl_mb(i_sti))
c
      return
      end
      subroutine cf_rmri(lfnmri,dera)
c
      implicit none
c
#include "cf_common.fh"
c
      integer lfnmri
      real*8 dera(5,nsatot)
c
      integer i,j
c
      if(me.eq.0) then
      read(lfnmri) ((dera(i,j),i=1,5),j=1,nsatot)
      else
      do 1 i=1,nsatot
      do 2 j=1,5
      dera(j,i)=zero
    2 continue
    1 continue
      endif
c
      return
      end
      subroutine cf_wrtgib(lfngib)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfngib
c
      call cf_wtgib(lfngib,dbl_mb(i_sti))
c
      return
      end
      subroutine cf_wtgib(lfngib,dera)
c
      implicit none
c
#include "cf_common.fh"
#include "msgids.fh"
c
      integer lfngib
      real*8 dera(5,nsatot)
c
      integer i,j
c
      call ga_dgop(mcf_68,dera,5*nsatot,'+')
c
      if(me.eq.0) then
      write(lfngib,1000) ((dera(i,j),i=1,5),j=1,nsatot)
 1000 format(4e20.10)
      else
      do 1 i=1,nsatot
      do 2 j=1,5
      dera(j,i)=zero
    2 continue
    1 continue
      endif
c
      return
      end
      subroutine cf_rdgib(lfngib,dec,dlambd)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfngib
      real*8 dec(5,nsatot),dlambd
c
      if(npgdec.eq.0) call errquit('rdgib error',0)
c
      call cf_rgib(lfngib,dbl_mb(i_sti),dec,dlambd)
c
      return
      end
      subroutine cf_rgib(lfngib,dera,dec,dlambd)
c
      implicit none
c
#include "cf_common.fh"
c
      integer lfngib
      real*8 dera(5,nsatot),dec(5,nsatot),dlambd
c
      integer i,j
c
      if(me.eq.0) then
      read(lfngib,1000) ((dera(i,j),i=1,5),j=1,nsatot)
 1000 format(4e20.10)
      do 1 i=1,nsatot
      do 2 j=1,5
      dec(j,i)=dec(j,i)+dlambd*dera(j,i)
    2 continue
    1 continue
      else
      do 3 i=1,nsatot
      do 4 j=1,5
      dera(j,i)=zero
    4 continue
    3 continue
      endif
c
      return
      end
      subroutine cf_print_deco(lfnout,dec)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfnout
      real*8 dec(5,nsatot)
c
      call cf_print_d(lfnout,byte_mb(i_snam),dec)
c
      return
      end
      subroutine cf_print_d(lfnout,snam,dec)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfnout
      real*8 dec(5,nsatot)
      character*16 snam(nsatot)
c
      integer i,j
c
      write(lfnout,1001)
 1001 format(//,' SOLUTE DECOMPOSITION THERMODYNAMIC INTEGRATION',//,
     + '     Atom',13x,'Solvent     Solvent      Solute      Solute',
     + 6x,'Solute',/,
     + 21x,'v.d.Waals   electrost   v.d.Waals   electrost',6x,
     + 'bonded',/)
      do 1 i=1,nsatot
      write(lfnout,1009) snam(i),(dec(j,i),j=1,5)
 1009 format(1x,a16,5f12.3)
    1 continue
c
      return
      end
      subroutine cf_add_esa(esa)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 esa(nsa,2)
c
      if(numb(2).gt.0) call cf_add_ebonds(numb(2),mbt(2),mbp(2),
     + int_mb(i_ibnd(2)),dbl_mb(i_bnd(2)),dbl_mb(i_rbnd(2)),esa)
c
      if(numh(2).gt.0) call cf_add_eangles(numh(2),mht(2),mhp(2),
     + int_mb(i_iang(2)),dbl_mb(i_ang(2)),dbl_mb(i_rang(2)),esa)
c
      if(numd(2).gt.0) call cf_add_edihes(numd(2),mdt(2),mdp(2),
     + int_mb(i_idih(2)),dbl_mb(i_dih(2)),dbl_mb(i_rdih(2)),esa)
c
      if(numi(2).gt.0) call cf_add_eimprs(numi(2),mit(2),mip(2),
     + int_mb(i_iimp(2)),dbl_mb(i_imp(2)),dbl_mb(i_rimp(2)),esa)
c
      return
      end
      subroutine cf_add_ebonds(nsb,msb,msp,ibnd,bnd,rbnd,esa)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nsb,msb,msp
      integer ibnd(msb,3)
      real*8 bnd(msb,msp,6),rbnd(msb,2),esa(nsa,2)
c
      integer i
c
      do 1 i=1,nsb
      if(ibnd(i,3).eq.0) then
      esa(ibnd(i,1),2)=esa(ibnd(i,1),2)+half*rbnd(i,2)
      esa(ibnd(i,2),2)=esa(ibnd(i,2),2)+half*rbnd(i,2)
      endif
    1 continue
c
      return
      end
      subroutine cf_add_eangles(nsh,msh,msp,iang,ang,rang,esa)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nsh,msh,msp
      integer iang(msh,4)
      real*8 ang(msh,msp,6),rang(msh,2),esa(nsa,2)
c
      integer i
c
      do 1 i=1,nsh
      esa(iang(i,1),2)=esa(iang(i,1),2)+third*rang(i,2)
      esa(iang(i,2),2)=esa(iang(i,2),2)+third*rang(i,2)
      esa(iang(i,3),2)=esa(iang(i,3),2)+third*rang(i,2)
    1 continue
c
      return
      end
      subroutine cf_add_edihes(nsd,msd,msp,idih,dih,rdih,esa)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nsd,msd,msp
      integer idih(msd,5)
      real*8 dih(msd,msp,6),rdih(msd,2),esa(nsa,2)
c
      integer i
c
      do 1 i=1,nsd
      esa(idih(i,1),2)=esa(idih(i,1),2)+fourth*rdih(i,2)
      esa(idih(i,2),2)=esa(idih(i,2),2)+fourth*rdih(i,2)
      esa(idih(i,3),2)=esa(idih(i,3),2)+fourth*rdih(i,2)
      esa(idih(i,4),2)=esa(idih(i,4),2)+fourth*rdih(i,2)
    1 continue
c
      return
      end
      subroutine cf_add_eimprs(nso,mso,msp,ioop,oop,roop,esa)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nso,mso,msp
      integer ioop(mso,5)
      real*8 oop(mso,msp,6),roop(mso,2),esa(nsa,2)
c
      integer i
c
      do 1 i=1,nso
      esa(ioop(i,1),2)=esa(ioop(i,1),2)+fourth*roop(i,2)
      esa(ioop(i,2),2)=esa(ioop(i,2),2)+fourth*roop(i,2)
      esa(ioop(i,3),2)=esa(ioop(i,3),2)+fourth*roop(i,2)
      esa(ioop(i,4),2)=esa(ioop(i,4),2)+fourth*roop(i,2)
    1 continue
c
      return
      end
