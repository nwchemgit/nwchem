      subroutine crossing(c,c_tst,v_tst,side,side_prev,nxing,forward,backward,forside,backside)
      implicit none

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
      logical forward,backward
      double precision dot,c(natom),c_tst(natom),v_tst(natom)
      integer natom,i,nxing,side,side_prev
      integer forside,backside

C ...c is the current coordinate, c_tst is teh coordinate of the TST
C ...v_tst is the imaginary mode vector 
C ... you need to know on which side were you at the previous(side_prev)point, so that you can determine if 
C ... you have crossed accross the dividing surface at this point (side)
  
      side_prev=side
      dot =0.0
      do i=1,3*natom
         dot=dot+(c(i)-c_tst(i))*v_tst(i)
      enddo
C
      if abs(dot).lt.small then
         side=0
      else
         if(dot.lt.0.0) then
            side=-1
         else
            side=1
         endif
      endif
C
C     ----- Side ? Xing ? ... -----
C
      if ((side.ne.side_prev).and.(side_prev.ne.0)) then
         nxing=nxing+1
      endif

      if (forward) then
        forside=side
      endif
      if (backward) then
        backside=side
      endif

      end
