      subroutine argos_cafe_hopinit(lsthop,timhop)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "util.fh"
c
      integer lsthop(2,*)
      real*8 timhop(*)
c
      integer i,n
c
      n=nhops
      nhops=0
      do 1 i=1,n
      if(stime-timhop(i).le.thop) then
      nhops=nhops+1
      lsthop(1,nhops)=lsthop(1,i)
      lsthop(2,nhops)=lsthop(2,i)
      timhop(nhops)=timhop(i)
      if(me.eq.0) then
      if(util_print('qhop',print_debug)) then
      write(lfnhop,'(a,3i5,f12.6)')
     + ' Possible backhop ',nhops,lsthop(1,i),lsthop(2,i),
     + stime-timhop(i)
      endif
      endif
      endif
    1 continue
c
      return
      end
c $Id$
