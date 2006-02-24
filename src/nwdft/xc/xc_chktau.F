      logical function xc_chktau()
C
C $Id: xc_chktau.F,v 1.4 2006-02-24 20:46:44 edo Exp $ 
C
      implicit none
 
c
c     see if any metaGGA Functional using TAU is not null
c
 
 
#include "cdft.fh"
C  note in cdft.fh,they define cfac(numfunc) and xfac(numfunc)
 
 
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(18) + cfac(25)+xfac(21)+cfac(27) +
     +     xfac(28) + xfac(29) +
     +     cfac(31) + cfac(32) + cfac(33) + cfac(34) + cfac(35)
      if (abs(tot).gt.eps)then
         xc_chktau = .true.
      else
         xc_chktau = .false.
      endif
      return
      end



