      subroutine cf_alpha
c $Id$
      implicit none
c
#include "cf_common.fh"
#include "cf_funcs_dec.fh"
c
      real*8 a,alow,ahigh,b
      integer i,n
c
#include "cf_funcs_sfn.fh"
c
      a=half
      i=0
    1 continue
      a=two*a
      i=i+1
      b=(erfc(a*rlong))/rlong
      if(b.ge.pmetol) goto 1
c
      n=i+60
      alow=zero
      ahigh=a
      do 2 i=1,n
      a=half*(alow+ahigh)
      b=(erfc(a*rlong))/rlong
      if(b.ge.pmetol) then
      alow=a
      else
      ahigh=a
      endif
    2 continue
      ealpha=a
c
      return
      end
      subroutine cf_chgrid(lset,ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + grid,mgz,lgrid,knode,kfrom,lsldo,node,mnx,mxx,mny,mxy,mnz,mxz)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      logical lset,lgrid
      integer ipset,nwmloc,nsaloc,node,mgz
      real*8 xw(mwm,3,mwa),xs(msa,3)
      integer isq(msa),isgm(msa)
      integer knode(*),kfrom(*)
      logical lsldo(*)
      real*8 grid(ngrx,ngry,mgz)
      integer mnx,mxx,mny,mxy,mnz,mxz
c
      if(lset.or..not.lscr) then
      if(.not.lscr) call cf_salloc
      call cf_pmeset(xw,dbl_mb(i_pmew),nwmloc,xs,dbl_mb(i_pmes),
     + nsaloc,dbl_mb(i_theta),dbl_mb(i_dtheta))
      lset=.false.
      endif
c
      call cf_grid(ipset,xw,dbl_mb(i_pmew),int_mb(i_iwq),nwmloc,
     + xs,dbl_mb(i_pmes),isq,isgm,nsaloc,
     + dbl_mb(i_chg),dbl_mb(i_theta),dbl_mb(i_dtheta),
     + grid,mgz,lgrid,knode,kfrom,lsldo,node,mnx,mxx,mny,mxy,mnz,mxz,
     + int_mb(i_lseq))
c
      return
      end
      subroutine cf_grid(ipset,xw,yw,iwq,nwmloc,xs,ys,isq,isgm,
     + nsaloc,chg,theta,dtheta,
     + grid,mgz,lgrid,knode,kfrom,lsldo,node,
     + minx,maxx,miny,maxy,minz,maxz,lseq)
c
      implicit none
c
#include "cf_common.fh"
c
      logical lgrid
      integer ipset,nwmloc,nsaloc,node,mgz
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 yw(mwm,3,mwa),ys(msa,3)
      integer iwq(mwa),isq(msa),isgm(msa),lseq(mseq)
      real*8 chg(mqt,mqp,mset),theta(morder,3,*),dtheta(morder,3,*)
      integer knode(np),kfrom(np)
      logical lsldo(np,*)
      real*8 grid(ngrx,ngry,mgz)
c
      integer i,i0,j,j0,k,k0,ith1,ith2,ith3,isa,iwm,iwa,indx
      integer minx,miny,minz,maxx,maxy,maxz
      real*8 prod,charge
c
      do 1 k=1,mgz
      do 21 j=1,ngry
      do 31 i=1,ngrx
      grid(i,j,k)=zero
   31 continue
   21 continue
    1 continue
      lgrid=.false.
      indx=0
c
      minx=ngrx
      miny=ngry
      minz=mgz
      maxx=1
      maxy=1
      maxz=1
c
      do 2 iwa=1,mwa
      do 3 iwm=1,nwmloc
      indx=indx+1
      k0=int(yw(iwm,3,iwa))-morder
      do 4 ith3=1,morder
      k0=k0+1
      k=k0+1+(ngz-isign(ngz,k0))/2
      if(k.gt.ngz) k=k-ngz
      if(k.lt.1) k=k+ngz
      if(knode(k).eq.node) then
      lgrid=.true.
      k=k-kfrom(node+1)+1
      j0=int(yw(iwm,2,iwa))-morder
      if(k.lt.minz) minz=k
      if(k.gt.maxz) maxz=k
      lsldo(node+1,k)=.true.
      do 5 ith2=1,morder
      j0=j0+1
      j=j0+1+(ngy-isign(ngy,j0))/2
      miny=min(miny,j)
      maxy=max(maxy,j)
c      j=j0+1+(ngy-isign(ngy,j0))/2
c      if(j.gt.ngy) j=j-ngy
c      if(j.lt.1) j=j+ngy
      prod=theta(ith2,2,indx)*theta(ith3,3,indx)*chg(iwq(iwa),1,ipset)
      i0=int(yw(iwm,1,iwa))-morder
c      if(j.lt.miny) miny=j
c      if(j.gt.maxy) maxy=j
      do 6 ith1=1,morder
      i0=i0+1
      i=i0+1+(ngx-isign(ngx,i0))/2
      minx=min(minx,i)
      maxx=max(maxx,i)
c      i=i0+1+(ngx-isign(ngx,i0))/2
c      if(i.gt.ngx) i=i-ngx
c      if(i.lt.1) i=i+ngx
c      if(i.lt.minx) minx=i
c      if(i.gt.maxx) maxx=i
      grid(i,j,k)=grid(i,j,k)+theta(ith1,1,indx)*prod
    6 continue
    5 continue
      endif
    4 continue
    3 continue
    2 continue
c
      if(isolvo.eq.0) then
      do 7 isa=1,nsaloc
      if(nfhop.eq.0) then
      charge=chg(isq(isa),1,ipset)
      else
      charge=chg(isq(isa),1,lseq(isgm(isa)))
      endif
c      write(*,'(2i5,f12.5)') isa,isgm(isa),charge
      indx=indx+1
      k0=int(ys(isa,3))-morder
      do 8 ith3=1,morder
      k0=k0+1
      k=k0+1+(ngz-isign(ngz,k0))/2
      if(knode(k).eq.node) then
      lgrid=.true.
      k=k-kfrom(node+1)+1
      j0=int(ys(isa,2))-morder
      if(k.lt.minz) minz=k
      if(k.gt.maxz) maxz=k
      lsldo(node+1,k)=.true.
      do 9 ith2=1,morder
      j0=j0+1
      j=j0+1+(ngy-isign(ngy,j0))/2
      miny=min(miny,j)
      maxy=max(maxy,j)
c      j=j0+1+(ngy-isign(ngy,j0))/2
c      if(j.gt.ngy) j=j-ngy
c      if(j.lt.1) j=j+ngy
      prod=theta(ith2,2,indx)*theta(ith3,3,indx)*charge
      i0=int(ys(isa,1))-morder
c      if(j.lt.miny) miny=j
c      if(j.gt.maxy) maxy=j
      do 10 ith1=1,morder
      i0=i0+1
      i=i0+1+(ngx-isign(ngx,i0))/2
      minx=min(minx,i)
      maxx=max(maxx,i)
c      i=i0+1+(ngx-isign(ngx,i0))/2
c      if(i.gt.ngx) i=i-ngx
c      if(i.lt.1) i=i+ngx
c      if(i.lt.minx) minx=i
c      if(i.gt.maxx) maxx=i
      grid(i,j,k)=grid(i,j,k)+theta(ith1,1,indx)*prod
   10 continue
    9 continue
      endif
    8 continue
    7 continue
      endif
c
      return
      end
      subroutine cf_pmeset(xw,yw,nwmloc,xs,ys,nsaloc,theta,dtheta)
c
      implicit none
c
#include "cf_common.fh"
c
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 yw(mwm,3,mwa),ys(msa,3)
      real*8 theta(morder,3,*),dtheta(morder,3,*)
      integer nwmloc,nsaloc
c
      integer i,j,iwa,ix,iwm,isa,indx,ng(3)
      real*8 w,rinv
c
      if(nbxtyp.ne.1) then
      do 1 i=1,3
      do 2 j=1,3
      recip(i,j)=zero
    2 continue
      recip(i,i)=one/box(i)
    1 continue
      else
      do 3 i=1,3
      do 4 j=1,3
      recip(i,j)=vlati(i,j)
    4 continue
    3 continue
      endif
c
c     calculate scaled fractional coordinates
c
      ng(1)=ngx
      ng(2)=ngy
      ng(3)=ngz
      do 5 iwa=1,mwa
      do 6 ix=1,3
      do 7 iwm=1,nwmloc
      w=xw(iwm,1,iwa)*recip(1,ix)+xw(iwm,2,iwa)*recip(2,ix)+
     + xw(iwm,3,iwa)*recip(3,ix)
      yw(iwm,ix,iwa)=dble(ng(ix))*(w-anint(w)+half)
    7 continue
    6 continue
    5 continue
c
      do 8 ix=1,3
      do 9 isa=1,nsaloc
      w=xs(isa,1)*recip(1,ix)+xs(isa,2)*recip(2,ix)+
     + xs(isa,3)*recip(3,ix)
      ys(isa,ix)=dble(ng(ix))*(w-anint(w)+half)
    9 continue
    8 continue
c
c     spline coefficients
c
      indx=0
      do 10 iwa=1,mwa
      do 11 iwm=1,nwmloc
      indx=indx+1
      do 12 ix=1,3
      w=yw(iwm,ix,iwa)-int(yw(iwm,ix,iwa))
      theta(morder,ix,indx)=zero
      theta(2,ix,indx)=w
      theta(1,ix,indx)=one-w
      do 13 i=3,morder-1
      rinv=one/dble(i-1)
      theta(i,ix,indx)=rinv*w*theta(i-1,ix,indx)
      do 14 j=1,i-2
      theta(i-j,ix,indx)=rinv*((w+dble(j))*theta(i-j-1,ix,indx)+
     + (dble(i-j)-w)*theta(i-j,ix,indx))
   14 continue
      theta(1,ix,indx)=rinv*(one-w)*theta(1,ix,indx)
   13 continue
      dtheta(1,ix,indx)=-theta(1,ix,indx)
      do 15 j=2,morder
      dtheta(j,ix,indx)=theta(j-1,ix,indx)-theta(j,ix,indx)
   15 continue
      rinv=one/dble(morder-1)
      theta(morder,ix,indx)=rinv*w*theta(morder-1,ix,indx)
      do 16 j=1,morder-2
      theta(morder-j,ix,indx)=rinv*((w+dble(j))*
     + theta(morder-j-1,ix,indx)+
     + (dble(morder-j)-w)*theta(morder-j,ix,indx))
   16 continue
      theta(1,ix,indx)=rinv*(one-w)*theta(1,ix,indx)
   12 continue
   11 continue
   10 continue
c      
      do 17 isa=1,nsaloc
      indx=indx+1
      do 18 ix=1,3
      w=ys(isa,ix)-int(ys(isa,ix))
      theta(morder,ix,indx)=zero
      theta(2,ix,indx)=w
      theta(1,ix,indx)=one-w
      do 19 i=3,morder-1
      rinv=one/dble(i-1)
      theta(i,ix,indx)=rinv*w*theta(i-1,ix,indx)
      do 20 j=1,i-2
      theta(i-j,ix,indx)=rinv*((w+dble(j))*theta(i-j-1,ix,indx)+
     + (dble(i-j)-w)*theta(i-j,ix,indx))
   20 continue
      theta(1,ix,indx)=rinv*(one-w)*theta(1,ix,indx)
   19 continue
      dtheta(1,ix,indx)=-theta(1,ix,indx)
      do 21 j=2,morder
      dtheta(j,ix,indx)=theta(j-1,ix,indx)-theta(j,ix,indx)
   21 continue
      rinv=one/dble(morder-1)
      theta(morder,ix,indx)=rinv*w*theta(morder-1,ix,indx)
      do 22 j=1,morder-2
      theta(morder-j,ix,indx)=rinv*((w+dble(j))*
     + theta(morder-j-1,ix,indx)+
     + (dble(morder-j)-w)*theta(morder-j,ix,indx))
   22 continue
      theta(1,ix,indx)=rinv*(one-w)*theta(1,ix,indx)
   18 continue
   17 continue
c
      return
      end
      subroutine cf_epme(knode,kfrom,kto,bmod,grid,mgz,energy)
c
      implicit none
c
#include "cf_common.fh"
#include "msgids.fh"
c
      integer mgz
      integer knode(ngz),kfrom(np),kto(np)
      real*8 bmod(ngmax,3),grid(2,ngx,ngy,mgz),energy
c
      real*8 fac,rh1,rh2,rh3,rsq,efac,vfac,ss
      integer i,j,k,kg,m1,m2,m3
      integer ngxy,nfx,nfy,nfz
c
      fac=(pi/ealpha)**2
c
      ngxy=ngx*ngy
      nfx=ngx/2
      if(2*nfx.lt.ngx) nfx=nfx+1
      nfy=ngy/2
      if(2*nfy.lt.ngy) nfy=nfy+1
      nfz=ngz/2
      if(2*nfz.lt.ngz) nfz=nfz+1
c
      epme=zero
c
      if(kfrom(me+1).gt.0) then
      do 2 k=kfrom(me+1),kto(me+1)
      kg=k+1-kfrom(me+1)
c
      do 3 j=1,ngy
      do 4 i=1,ngx
      if(i+j+k.gt.3) then
c
      m1=i-1
      if(i.gt.nfx) m1=m1-ngx
      m2=j-1
      if(j.gt.nfy) m2=m2-ngy
      m3=k-1
      if(k.gt.nfz) m3=m3-ngz
c
      rh1=recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
      rh2=recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
      rh3=recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
      rsq=rh1*rh1+rh2*rh2+rh3*rh3
c
      efac=exp(-fac*rsq)/(pi*volume*bmod(i,1)*bmod(j,2)*bmod(k,3)*rsq)
      vfac=two*(fac*rsq+one)/rsq
      ss=grid(1,i,j,kg)**2+grid(2,i,j,kg)**2
c
      epme=epme+efac*ss
      vpme(1)=vpme(1)+efac*ss*(vfac*rh1*rh1-one)
      vpme(2)=vpme(2)+efac*ss*(vfac*rh1*rh2)
      vpme(3)=vpme(3)+efac*ss*(vfac*rh1*rh3)
      vpme(4)=vpme(4)+efac*ss*(vfac*rh2*rh2-one)
      vpme(5)=vpme(5)+efac*ss*(vfac*rh2*rh3)
      vpme(6)=vpme(6)+efac*ss*(vfac*rh3*rh3-one)
c
      grid(1,i,j,kg)=efac*grid(1,i,j,kg)
      grid(2,i,j,kg)=efac*grid(2,i,j,kg)
c
      endif
    4 continue
    3 continue
c
    2 continue
      endif
c
      energy=0.5d0*epme
c
      return
      end
      subroutine cf_fpme(fw,nwmloc,fs,isq,isgm,nsaloc,grid,mgz,
     + knode,kfrom,lsldo,node)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer nwmloc,nsaloc,node,mgz
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer isq(msa),isgm(msa)
      integer knode(*),kfrom(*)
      logical lsldo(*)
      real*8 grid(ngrx,ngry,mgz)
c
      call cf_forpme(dbl_mb(i_pmew),fw,int_mb(i_iwq),nwmloc,
     + dbl_mb(i_pmes),fs,isq,isgm,nsaloc,
     + dbl_mb(i_chg),dbl_mb(i_theta),dbl_mb(i_dtheta),
     + grid,mgz,knode,kfrom,lsldo,node,int_mb(i_lseq))
c
      return
      end
      subroutine cf_forpme(yw,fw,iwq,nwmloc,ys,fs,isq,isgm,nsaloc,
     + chg,theta,dtheta,
     + grid,mgz,knode,kfrom,lsldo,node,lseq)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer nwmloc,nsaloc,node,mgz
      real*8 yw(mwm,3,mwa),fw(mwm,3,mwa),ys(msa,3),fs(msa,3)
      integer iwq(mwa),isq(msa),isgm(msa),lseq(mseq)
      real*8 chg(mqt,mqp,mset),theta(morder,3,*),dtheta(morder,3,*)
      integer knode(*),kfrom(*)
      logical lsldo(*)
      real*8 grid(ngrx,ngry,mgz)
c
      integer i,i0,j,j0,k,k0,ith1,ith2,ith3,iwm,iwa,isa,indx
      real*8 term,fx,fy,fz,charge
c
      indx=0
c
      do 2 iwa=1,mwa
      do 3 iwm=1,nwmloc
      indx=indx+1
      fx=zero
      fy=zero
      fz=zero
      k0=int(yw(iwm,3,iwa))-morder
      do 4 ith3=1,morder
      k0=k0+1
      k=k0+1+(ngz-isign(ngz,k0))/2
      if(k.gt.ngz) k=k-ngz
      if(k.lt.1) k=k+ngz
      if(knode(k).eq.node) then
      k=k-kfrom(node+1)+1
      j0=int(yw(iwm,2,iwa))-morder
      do 5 ith2=1,morder
      j0=j0+1
      j=j0+1+(ngy-isign(ngy,j0))/2
      i0=int(yw(iwm,1,iwa))-morder
      do 6 ith1=1,morder
      i0=i0+1
      i=i0+1+(ngx-isign(ngx,i0))/2
      term=chg(iwq(iwa),1,iset)*grid(i,j,k)
      fx=fx-ngx*term*dtheta(ith1,1,indx)*theta(ith2,2,indx)*
     + theta(ith3,3,indx)
      fy=fy-ngy*term*dtheta(ith2,2,indx)*theta(ith3,3,indx)*
     + theta(ith1,1,indx)
      fz=fz-ngz*term*dtheta(ith3,3,indx)*theta(ith1,1,indx)*
     + theta(ith2,2,indx)
    6 continue
    5 continue
      endif
    4 continue
      fw(iwm,1,iwa)=fw(iwm,1,iwa)+recip(1,1)*fx+recip(1,2)*fy+
     + recip(1,3)*fz
      fw(iwm,2,iwa)=fw(iwm,2,iwa)+recip(2,1)*fx+recip(2,2)*fy+
     + recip(2,3)*fz
      fw(iwm,3,iwa)=fw(iwm,3,iwa)+recip(3,1)*fx+recip(3,2)*fy+
     + recip(3,3)*fz
    3 continue
    2 continue
c
      do 7 isa=1,nsaloc
      if(nfhop.eq.0) then
      charge=chg(isq(isa),1,iset)
      else
      charge=chg(isq(isa),1,lseq(isgm(isa)))
      endif
      indx=indx+1
      fx=zero
      fy=zero
      fz=zero
      k0=int(ys(isa,3))-morder
      do 8 ith3=1,morder
      k0=k0+1
      k=k0+1+(ngz-isign(ngz,k0))/2
      if(knode(k).eq.node) then
      k=k-kfrom(node+1)+1
      j0=int(ys(isa,2))-morder
      do 9 ith2=1,morder
      j0=j0+1
      j=j0+1+(ngy-isign(ngy,j0))/2
      i0=int(ys(isa,1))-morder
      do 10 ith1=1,morder
      i0=i0+1
      i=i0+1+(ngx-isign(ngx,i0))/2
      term=charge*grid(i,j,k)
      fx=fx-ngx*term*dtheta(ith1,1,indx)*theta(ith2,2,indx)*
     + theta(ith3,3,indx)
      fy=fy-ngy*term*dtheta(ith2,2,indx)*theta(ith3,3,indx)*
     + theta(ith1,1,indx)
      fz=fz-ngz*term*dtheta(ith3,3,indx)*theta(ith1,1,indx)*
     + theta(ith2,2,indx)
   10 continue
    9 continue
      endif
    8 continue
      fs(isa,1)=fs(isa,1)+recip(1,1)*fx+recip(1,2)*fy+recip(1,3)*fz
      fs(isa,2)=fs(isa,2)+recip(2,1)*fx+recip(2,2)*fy+recip(2,3)*fz
      fs(isa,3)=fs(isa,3)+recip(3,1)*fx+recip(3,2)*fy+recip(3,3)*fz
    7 continue
c
      return
      end
      subroutine cf_pmecor(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,ecorr)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer ipset,nwmloc,nsaloc
      real*8 xw(mwm,3,mwa),xs(msa,3),ecorr
      integer isq(msa),isgm(msa)
c
      call cf_pmec(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + ecorr,int_mb(i_iwq),dbl_mb(i_chg),int_mb(i_lseq))
c
      return
      end
      subroutine cf_pmec(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + ecorr,iwq,chg,lseq)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer ipset,nwmloc,nsaloc
      real*8 xw(mwm,3,mwa),xs(msa,3),ecorr
      integer isq(msa),iwq(mwa),isgm(msa),lseq(mseq)
      real*8 chg(mqt,mqp,mset)
c
      integer i,icset
c
      ecorr=zero
c
      if(nwmloc.gt.0) then
      do 1 i=1,nwa
      ecorr=ecorr-chg(iwq(i),1,ipset)*chg(iwq(i),1,ipset)
    1 continue
      ecorr=ecorr*dble(nwmloc)
      endif
c
      if(nfhop.eq.0) then
      if(nsaloc.gt.0) then
      do 2 i=1,nsaloc
      ecorr=ecorr-chg(isq(i),1,ipset)*chg(isq(i),1,ipset)
    2 continue
      endif
      else
      do 3 i=1,nsaloc
      icset=lseq(isgm(i))
      ecorr=ecorr-chg(isq(i),1,icset)*chg(isq(i),1,icset)
    3 continue
      endif
c
      ecorr=ealpha*ecorr/sqrt(pi)
c
      return
      end
