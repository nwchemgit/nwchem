      subroutine argos_cafe_forpme(yw,fw,iwq,nwmloc,ys,fs,
     + isq,isgm,nsaloc,chg,theta,dtheta,
     + grid,mgz,knode,kfrom,lsldo,node,lseq)
c
      implicit none
c
#include "argos_cafe_common.fh"
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
      j0=int(yw(iwm,2,iwa))-1
      do 5 ith2=1,morder
      j0=j0+1
      j=j0+1
c      j=j0+1+(ngy-isign(ngy,j0))/2
c      if(j.gt.ngy) j=j-ngy
c      if(j.lt.1) j=j+ngy
      i0=int(yw(iwm,1,iwa))-1
      do 6 ith1=1,morder
      i0=i0+1
      i=i0+1
c      i=i0+1+(ngx-isign(ngx,i0))/2
c      if(i.gt.ngx) i=i-ngx
c      if(i.lt.1) i=i+ngx
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
      j0=int(ys(isa,2))-1
      do 9 ith2=1,morder
      j0=j0+1
      j=j0+1
c      j=j0+1+(ngy-isign(ngy,j0))/2
      i0=int(ys(isa,1))-1
      do 10 ith1=1,morder
      i0=i0+1
      i=i0+1
c      i=i0+1+(ngx-isign(ngx,i0))/2
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
c $Id$
