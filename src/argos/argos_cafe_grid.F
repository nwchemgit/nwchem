      subroutine argos_cafe_grid(ipset,xw,yw,iwq,nwmloc,xs,ys,isq,isgm,
     + nsaloc,chg,theta,dtheta,
     + grid,mgz,lgrid,knode,kfrom,lsldo,node,
     + minx,maxx,miny,maxy,minz,maxz,lseq)
c
      implicit none
c
#include "argos_cafe_common.fh"
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
      j0=int(yw(iwm,2,iwa))-1
      if(k.lt.minz) minz=k
      if(k.gt.maxz) maxz=k
      lsldo(node+1,k)=.true.
      miny=min(miny,j0+2)
      maxy=max(maxy,j0+1+morder)
      do 5 ith2=1,morder
      j0=j0+1
      j=j0+1
c      j=j0+1+(ngy-isign(ngy,j0))/2
c      if(j.gt.ngy) j=j-ngy
c      if(j.lt.1) j=j+ngy
      prod=theta(ith2,2,indx)*theta(ith3,3,indx)*chg(iwq(iwa),1,ipset)
      i0=int(yw(iwm,1,iwa))-1
c      if(j.lt.miny) miny=j
c      if(j.gt.maxy) maxy=j
      minx=min(minx,i0+2)
      maxx=max(maxx,i0+1+morder)
      do 6 ith1=1,morder
      i0=i0+1
      i=i0+1
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
      j0=int(ys(isa,2))-1
      if(k.lt.minz) minz=k
      if(k.gt.maxz) maxz=k
      lsldo(node+1,k)=.true.
      miny=min(miny,j0+2)
      maxy=max(maxy,j0+1+morder)
      do 9 ith2=1,morder
      j0=j0+1
      j=j0+1
c      j=j0+1+(ngy-isign(ngy,j0))/2
c      if(j.gt.ngy) j=j-ngy
c      if(j.lt.1) j=j+ngy
      prod=theta(ith2,2,indx)*theta(ith3,3,indx)*charge
      i0=int(ys(isa,1))-1
c      if(j.lt.miny) miny=j
c      if(j.gt.maxy) maxy=j
      minx=min(minx,i0+2)
      maxx=max(maxx,i0+1+morder)
      do 10 ith1=1,morder
      i0=i0+1
      i=i0+1
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
c $Id$
