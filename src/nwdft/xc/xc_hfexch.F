      double precision function xc_hfexch()
      implicit none
c
c     This function returns the amount of Hartree-Fock exchange in a
c     "functional". If no functional is used (i.e. not a DFT
c     calculation) then the fraction of Hartree-Fock exchange is 1.0.
c
c
#include "cdft.fh"
c
      logical xc_gotxc
      external xc_gotxc
c
      xc_hfexch = 1.0d0
      if (xc_gotxc()) then
        xc_hfexch = xfac(1)
      endif
c
      return
      end
c $Id$
