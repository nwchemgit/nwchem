      logical function xc_chktau()
C
C $Id: xc_chktau.F,v 1.2 2003-05-27 22:35:13 pollack Exp $ 
C
      implicit none
 
c
c     see if any metaGGA Functional using TAU is not null
c
 
 
#include "cdft.fh"
C  note in cdft.fh,they define cfac(numfunc) and xfac(numfunc)
 
 
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(18) + cfac(25)
      if (abs(tot).gt.eps)then
         xc_chktau = .true.
      else
         xc_chktau = .false.
      endif
      return
      end



