      logical function xc_chkdispauto()
C
      implicit none
c
c     see if any functional has automatic dispersion included
c
#include "cdft.fh"
c     note in cdft.fh, they define cfac(numfunc) and xfac(numfunc)
c 
      logical lhasautodisp
      lhasautodisp = .false.
c
c     beck97-d, ssb-d
      lhasautodisp = xccomb(37).or.xccomb(46).or.xccomb(60).or.
     &               xccomb(61).or.xccomb(62).or.xccomb(63)
c
      if (lhasautodisp) then
         xc_chkdispauto = .true.
      else
         xc_chkdispauto = .false.
      endif
c
      return
      end
c $Id$
