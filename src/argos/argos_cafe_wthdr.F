      subroutine argos_cafe_wthdr(iunit,snam,idbw,nwb,idbs,nsb)
c
      implicit none
c
#include "argos_cafe_common.fh"
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
      write(iunit,1001) nwa,nsatot,numb(2),numb(1),nwmtot,nprec
 1001 format(6i10)
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
c $Id$
