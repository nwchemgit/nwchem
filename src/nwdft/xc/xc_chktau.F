      logical function xc_chktau()
C
C $Id: xc_chktau.F,v 1.8 2009-01-27 01:01:02 niri Exp $ 
C
      implicit none
 
c
c     see if any metaGGA Functional using TAU is not null
c
 
 
#include "cdft.fh"
C  note in cdft.fh,they define cfac(numfunc) and xfac(numfunc)
 
 
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(18) + cfac(25) + xfac(21) + cfac(27) 
     +       + xfac(28) + xfac(29) + xfac(33) + xfac(34) + xfac(35)
     +       + xfac(36) + cfac(29) + cfac(31) + cfac(32) + cfac(33) 
     +       + cfac(34) + cfac(35) + cfac(30) + cfac(38)
     +       + cfac(39) + cfac(40)

      if (abs(tot).gt.eps)then
         xc_chktau = .true.
      else
         xc_chktau = .false.
      endif
      return
      end
      logical function xc_getxcfac(rtdb)
      implicit none
      integer rtdb
c     load cfac and xfac into cdft.fh
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      xc_getxcfac=rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc,
     &xfac)
      xc_getxcfac=xc_getxcfac.and.
     A     rtdb_get(rtdb, 'dft:cfac', mt_dbl, numfunc,
     &     cfac)
      return
      end



