      subroutine argos_cafe_final(xw,xwm,yw,vw,vwt,numw,
     + xs,xsm,ys,vs,vst,isat,ismol,isdyn,isfrc,ishop,zs,esk,nums,
     + bx,vl,prs,tmp,tmpw,tmps)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      real*8 pscalm
      real*8 xw(mwm,3,mwa),xwm(mwm,3),yw(mwm,3,mwa)
      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa)
      real*8 xs(msa,3),xsm(msm,3),ys(msa,3),vs(msa,3),vst(msa,3)
      real*8 zs(msf,3,3,2),bx(3),prs,tmp,tmpw,tmps,esk(msf),vl(3,3)
      integer isat(msa),ismol(msa),isdyn(msa),isfrc(msa),ishop(msa)
      integer numw,nums
c
      integer i,j,k
      real*8 prssc3,rdata(4),xwcd,xscd,eskin,ekinc
      real*8 pscx(3)
      real*8 pbuff(4,3,3)
      integer i_vsm,l_vsm
c
      call argos_cafe_vw(xw,yw,vw,vwt,numw,int_mb(i_iwa),dbl_mb(i_mas))
      if(.not.ma_push_get(mt_dbl,3*msm,'vsm',l_vsm,i_vsm))
     + call md_abort('Failed to allocate scratch array vsm',0)
      call argos_cafe_vs(xs,ys,vs,vst,esk,nums,isat,ismol,isfrc,ishop,
     + dbl_mb(i_wsm),dbl_mb(i_mas),dbl_mb(i_vsm))
      if(.not.ma_pop_stack(l_vsm))
     + call md_abort('Failed to deallocate scratch array vsm',0)
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
      if(ipopt.eq.3) then
      pscx(1)=one
      pscx(2)=one
      endif
      if(ipopt.eq.12) pscx(3)=one
      if(ipopt.eq.123) then
      pscx(1)=half*(pscx(1)+pscx(2))
      pscx(2)=pscx(1)
      endif
      pscale=pscx(1)*pscx(2)*pscx(3)
c
      if(np.gt.1) then
      rdata(1)=pscx(1)
      rdata(2)=pscx(2)
      rdata(3)=pscx(3)
      rdata(4)=pscale
      call ga_brdcst(mcf_39,rdata,ma_sizeof(mt_dbl,4,mt_byte),0)
c      write(*,'(i5,4f12.6)') me,rdata
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
      call argos_cafe_volume()
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
c $Id$
