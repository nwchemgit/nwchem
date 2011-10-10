      subroutine argos_wrtime
c
      implicit none
c
#include "argos.fh"
#include "global.fh"
#include "msgids.fh"
c
      external timer_wall,timer_wall_total
      real*8 timer_wall,timer_wall_total
c
      integer i,j
      real*8 tim(56,1024)
c
      do 1 i=1,np
      do 2 j=1,56
      tim(j,i)=zero
    2 continue
    1 continue
c
      do 3 i=1,55
      tim(i,me+1)=timer_wall_total(i)
      tim(56,me+1)=tim(56,me+1)+tim(i,me+1)
    3 continue
c
      call ga_dgop(mrg_d04,tim,56*np,'+')
c
      if(me.eq.0) then
      write(lfntim,1000) stime
 1000 format('timings',/,f12.6)
      do 1002 j=1,np
      write(lfntim,1001) (tim(i,j),i=1,56)
 1001 format(10f7.3)
 1002 continue
      endif
c
      return
      end
c $Id$
