      subroutine argos_diana_evalord(xtrj)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "argos_diana_common.fh"
c
      real*8 xtrj(ntrj,3,nfrdim)
c
      integer i,j,k,ndordr,ndxo,inum
      real*8 s,d
      logical lopen
c
      do 4 i=1,nord
      rord(i)=0.0d0
    4 continue
c
      do 1 i=1,nord
      ndxo=iord(i,1,3)
      ndordr=iord(i,3,3)
      do 2 k=1,nfrme
      do 3 j=2,ndordr+1
      d=(xtrj(iord(i,j,1),1,k)-xtrj(iord(i,1,1),1,k))**2+
     +(xtrj(iord(i,j,1),2,k)-xtrj(iord(i,1,1),2,k))**2+
     +(xtrj(iord(i,j,1),3,k)-xtrj(iord(i,1,1),3,k))**2
      s=3.0d0*(xtrj(iord(i,j,1),3,k)-xtrj(iord(i,1,1),3,k))**2/d-1.0d0
      rord(ndxo)=rord(ndxo)+s
    3 continue
    2 continue
    1 continue
c
      call ga_dgop(mag_d09,rord,mord,'+')
c
      if(me.eq.0) then
      inquire(unit=lfnord,opened=lopen)
      if(.not.lopen)
     + open(unit=lfnord, file=filord(1:index(filord,' ')-1),
     + form='formatted',status='unknown')
      k=0
      j=0
      do 5 i=1,nord
      if(iord(i,1,3).ne.k) then
      k=iord(i,1,3)
      j=j+1
      endif
    5 continue
c
      write(lfnord,1000) j,nord,nfrtot
 1000 format(2i5,i10)
c
      k=0
      j=0
      do 6 i=1,nord
      if(iord(i,1,3).ne.k) then
      k=iord(i,1,3)
      j=0
      endif
      j=j+1
      write(lfnord,1001) k,j,cord(i)
 1001 format(2i5,a)
    6 continue
c
      k=0
      do 7 i=1,nord
      if(iord(i,1,3).ne.k) then
      inum=2*iord(i,3,3)*iord(i,2,3)*nfrtot
      k=iord(i,1,3)
      write(lfnord,1002) k,iord(i,3,3),iord(i,2,3),rord(k)/dble(inum)
 1002 format(3i5,f12.5)
      endif
    7 continue
c
      close(unit=lfnord,status='keep')
c
      endif
c
      return
      end
c $Id$
