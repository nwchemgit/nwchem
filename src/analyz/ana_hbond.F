      subroutine ana_hbond_init(snam,ibnd,ihbnd)
c
c $Id: ana_hbond.F
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
c
      character*16 snam(msa)
      integer ibnd(msb,2),ihbnd(msa)
c
      integer i,j,k
c
      if(np.gt.1)
     + call errquit('hydrogen bond analysis on single processor only',0)
c
      do 1 i=1,nsa
      ihbnd(i)=0
    1 continue
c
      do 2 j=1,2
      do 3 i=1,msb
      if(snam(ibnd(i,j))(6:7).eq.' O'.or.
     + snam(ibnd(i,j))(6:7).eq.' N'.or.
     + snam(ibnd(i,j))(6:7).eq.' S') ihbnd(ibnd(i,j))=-1
    3 continue
    2 continue
c
      do 4 j=1,2
      k=2
      if(j.eq.2) k=1
      do 5 i=1,msb
      if(ihbnd(ibnd(i,k)).eq.-1) then
      if(snam(ibnd(i,j))(6:7).eq.' H'.or.
     + snam(ibnd(i,j))(6:7).eq.'2H'.or.
     + snam(ibnd(i,j))(6:7).eq.'3H'.or.
     + snam(ibnd(i,j))(6:7).eq.'4H'.or.
     + snam(ibnd(i,j))(6:7).eq.'5H'.or.
     + snam(ibnd(i,j))(6:7).eq.'6H'.or.
     + snam(ibnd(i,j))(6:7).eq.'7H'.or.
     + snam(ibnd(i,j))(6:7).eq.'8H'.or.
     + snam(ibnd(i,j))(6:7).eq.'9H') ihbnd(ibnd(i,j))=ibnd(i,k)
      endif
    5 continue
    4 continue
c
      do 6 i=1,nwa
      ihbndw(i)=0
    6 continue
c
      do 7 j=1,2
      do 8 i=1,mwb
      if(wnam(ibndw(i,j))(6:7).eq.' O'.or.
     + wnam(ibndw(i,j))(6:7).eq.' N'.or.
     + wnam(ibndw(i,j))(6:7).eq.' S') ihbndw(ibnd(i,j))=-1
    8 continue
    7 continue
c
      do 9 j=1,2
      k=2
      if(j.eq.2) k=1
      do 10 i=1,mwb
      if(ihbndw(ibndw(i,k)).eq.-1) then
      if(wnam(ibndw(i,j))(6:7).eq.' H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'2H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'3H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'4H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'5H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'6H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'7H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'8H'.or.
     + wnam(ibndw(i,j))(6:7).eq.'9H') ihbndw(ibndw(i,j))=ibndw(i,k)
      endif
   10 continue
    9 continue
c
      return
      end
      subroutine ana_hbond_header(snam,isel,ihbnd)
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
c
      character*16 snam(msa)
      integer isel(msa),ihbnd(msa)
c
      integer i
c
      write(lfnhba,1000) nsa
 1000 format('hydrogen bond analysis',/,i7)
      do 1 i=1,nsa
      write(lfnhba,1001) i,snam(i),isel(i),ihbnd(i)
 1001 format(i7,3x,a,2i7)
    1 continue
c
      return
      end
      subroutine ana_hbond(xs,xw,isel,ihbnd)
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
c
      real*8 ana_angle
      external ana_angle
c
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer ihbnd(msa),isel(msa)
c
      integer i,j,k
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
      a=ana_angle(xs,msa,xw,mwm,mwa,ihbnd(i),0,i,0,j,0)
      else
      a=ana_angle(xs,msa,xw,mwm,mwa,ihbnd(j),0,j,0,i,0)
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
      if(lxw) then
      write(lfnhba,1002)
 1002 format('solute-solvent')
c
      do 3 i=1,nsa
      do 4 j=1,nwm
      do 5 k=1,nwa
c
      if(ihbnd(i)*ihbndw(k).lt.0.and.isel(i).gt.0) then
      d=sqrt((xs(i,1)-xw(j,k,1))**2+(xs(i,2)-xw(j,k,2))**2+
     + (xs(i,3)-xw(j,k,3))**2)
      if(d.gt.rhbmin.and.d.lt.rhbmax) then
      if(ihbnd(i).gt.0) then
      a=ana_angle(xs,msa,xw,mwm,mwa,ihbnd(i),0,i,0,j,k)
      else
      a=ana_angle(xs,msa,xw,mwm,mwa,j,ihbndw(k),j,k,i,0)
      endif
      if(a.gt.hbdmin.and.a.lt.hbdmax) then
      write(lfnhba,1003) i,j,k,d,a
 1003 format(3i7,2f12.6)
      endif
      endif
      endif
c
    5 continue
    4 continue
    3 continue
c
      endif
c
      return
      end
