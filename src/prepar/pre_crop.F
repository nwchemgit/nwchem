      logical function pre_crop(xs,msa,nsa,lfnout)
c
      implicit none
c
      integer msa,nsa,lfnout
      real*8 xs(3,msa)
c
      integer i,k
      real*8 xmax(3),xmin(3)
c
      xmax(1)=xs(1,1)
      xmax(2)=xs(2,1)
      xmax(3)=xs(3,1)
      xmin(1)=xs(1,1)
      xmin(2)=xs(2,1)
      xmin(3)=xs(3,1)
c
      do 1 k=1,3
      do 2 i=1,nsa
      xmax(k)=max(xmax(k),xs(k,i))
      xmin(k)=min(xmin(k),xs(k,i))
    2 continue
    1 continue
c
      write(lfnout,1000) xmin,xmax
 1000 format(/,' Minimum coordinates ',3f12.6,/,
     + ' Maximum coordinates ',3f12.6,/)
c
      pre_crop=.true.
c
      return
      end

c $Id$
