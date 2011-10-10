      subroutine argos_diana_ramhdr(sgmnam,iram)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      integer iram(msgm,7)
c
      integer i,j,number
c
      number=0
      do 1 i=1,nsa
      if(i.eq.nsa) goto 2
      if(sgmnam(i)(11:16).ne.sgmnam(i+1)(11:16)) goto 2
      goto 1
    2 continue
      read(sgmnam(i)(11:16),'(i6)') j
      if(iram(j,1).gt.0.or.iram(j,5).gt.0) number=number+1
    1 continue
c
      write(lfnram,1000) number
 1000 format(i5)
c
      do 3 i=1,nsa
      if(i.eq.nsa) goto 4
      if(sgmnam(i)(11:16).ne.sgmnam(i+1)(11:16)) goto 4
      goto 3
    4 continue
      read(sgmnam(i)(11:16),'(i6)') j
      if(iram(j,1).gt.0.or.iram(j,5).gt.0)
     + write(lfnram,1002) sgmnam(i)(1:5),sgmnam(i)(11:16)
 1002 format(a5,a6)
    3 continue
c
      return
      end
c $Id$
