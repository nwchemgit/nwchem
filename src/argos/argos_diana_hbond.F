      subroutine argos_diana_hbond(xs,xw,isel,ihbnd,idone)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
c
      real*8 argos_diana_angle
      external argos_diana_angle
c
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer ihbnd(msa),isel(msa),idone(mwm)
c
      integer i,j,k,l,n
      real*8 d,a
c
      write(lfnhba,1000) time
 1000 format('frame',/,f12.6,/,'solute-solute')
      do 1 i=1,nsa-1
      do 2 j=i+1,nsa
      if(ihbnd(i)*ihbnd(j).lt.0.and.(isel(i).gt.0.or.isel(j).gt.0)) then
      if(ihbnd(i).lt.0.and.ihbnd(j).eq.i) goto 2
      if(ihbnd(j).lt.0.and.ihbnd(i).eq.j) goto 2
      d=sqrt((xs(i,1)-xs(j,1))**2+(xs(i,2)-xs(j,2))**2+
     + (xs(i,3)-xs(j,3))**2)
      if(d.gt.rhbmin.and.d.lt.rhbmax) then
      if(ihbnd(i).gt.0) then
      a=argos_diana_angle(xs,msa,xw,mwm,mwa,ihbnd(i),0,i,0,j,0)
      else
      a=argos_diana_angle(xs,msa,xw,mwm,mwa,ihbnd(j),0,j,0,i,0)
      endif
      if(a.gt.hbdmin.and.a.lt.hbdmax) then
      write(lfnhba,1001) i,j,d,a
 1001 format(2i7,2f12.6)
      endif
      endif
      endif
    2 continue
    1 continue
c
      if(lxw.and.numwhb.gt.0) then
      write(lfnhba,1002)
 1002 format('solute-solvent')
c
      do 3 i=1,nwm
      idone(i)=0
    3 continue
c
      do 4 i=1,nsa
      do 5 j=1,nwm
      do 6 k=1,nwa
c
      if(ihbnd(i)*ihbndw(k).lt.0.and.isel(i).gt.0) then
      d=sqrt((xs(i,1)-xw(j,k,1))**2+(xs(i,2)-xw(j,k,2))**2+
     + (xs(i,3)-xw(j,k,3))**2)
      if(d.gt.rhbmin.and.d.lt.rhbmax) then
      if(ihbnd(i).gt.0) then
      a=argos_diana_angle(xs,msa,xw,mwm,mwa,ihbnd(i),0,i,0,j,k)
      else
      a=argos_diana_angle(xs,msa,xw,mwm,mwa,j,ihbndw(k),j,k,i,0)
      endif
      if(a.gt.hbdmin.and.a.lt.hbdmax) then
      write(lfnhba,1003) i,j,k,d,a
 1003 format(3i7,2f12.6)
      idone(j)=1
      endif
      endif
      endif
c
    6 continue
    5 continue
    4 continue
c
      if(numwhb.gt.1) then
      write(lfnhba,1004)
 1004 format('solvent-solvent')
      do 7 n=1,numwhb-1
      do 8 i=1,nwm-1
      do 9 j=i+1,nwm
      if((idone(i).eq.n.and.idone(j).eq.0).or.
     + (idone(i).eq.0.and.idone(j).eq.n)) then
      do 10 k=1,nwa
      do 11 l=1,nwa
      if(ihbndw(k)*ihbndw(l).lt.0) then
      d=sqrt((xw(i,k,1)-xw(j,l,1))**2+(xw(i,k,2)-xw(j,l,2))**2+
     + (xw(i,k,3)-xw(j,l,3))**2)
      if(d.gt.rhbmin.and.d.lt.rhbmax) then
      if(ihbndw(i).gt.0) then
      a=argos_diana_angle(xs,msa,xw,mwm,mwa,i,ihbndw(i),i,k,j,l)
      else
      a=argos_diana_angle(xs,msa,xw,mwm,mwa,j,ihbndw(l),j,l,i,k)
      endif
      if(a.gt.hbdmin.and.a.lt.hbdmax) then
      write(lfnhba,1005) i,k,j,l,d,a,n
 1005 format(4i7,2f12.6,i5)
      if(idone(i).eq.0) idone(i)=n+1
      if(idone(j).eq.0) idone(j)=n+1
      endif
      endif
      endif
   11 continue
   10 continue
      endif
    9 continue
    8 continue
    7 continue
      endif
c
      endif
c
      return
      end
c $Id$
