      logical function xc_chkdispauto()
C
      implicit none
c
c     see if any functional has automatic dispersion included
c
#include "cdft.fh"
c     note in cdft.fh, they define cfac(numfunc) and xfac(numfunc)
c 
      double precision tot,eps
      parameter(eps=1.d-10)
c
c     beck97-d, ssb-d
      tot = xccomb(37) + xccomb(46)   
c
      if (abs(tot).gt.eps) then
         xc_chkdispauto = .true.
      else
         xc_chkdispauto = .false.
      endif
c
      return
      end
