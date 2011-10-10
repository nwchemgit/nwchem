      subroutine argos_cafe_epme(knode,kfrom,kto,bmod,grid,mgz,energy)
c
      implicit none
c
#include "argos_cafe_common.fh"
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
    2 continue
      endif
c
      energy=0.5d0*epme
c
      return
      end
c $Id$
