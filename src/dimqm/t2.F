      function t2r(r, dist, invdist, a1, a2, pol1, pol2) result (tr)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
c
c     Input variables
      integer, intent(in) ::  a1, a2
      double precision, intent(in) :: dist, invdist
      double precision, intent(in) :: pol1, pol2
      double precision, intent(in) :: r(3)
c
c     Local variables
      double precision dist3
      double precision invdist3
      double precision invdist5
      double precision temp1, temp2
      double precision Rp1, Rp2
      double precision Rpp, invRpp
      double precision aa, ab, ag
      double precision bb, bg
      double precision gg
c
c     Returns
      double precision tr(3,3)
      double precision util_erf
      external util_erf
c 
c     Compute powers and inverses once
      dist3 = dist * dist * dist
      invdist3 = invdist  * invdist * invdist
      invdist5 = invdist3 * invdist * invdist
c
c     Calculate R values
      Rp1   = (SQRTTWOINVPI * (pol1 * THIRD))**THIRD
      Rp2   = (SQRTTWOINVPI * (pol2 * THIRD))**THIRD
      Rpp   = SQRT(Rp1**2 + Rp2**2)
      invRpp = ONE / Rpp
c
c     Calculate screening values
      temp1 = util_erf(dist * invRpp) - TWOINVSQRTPI
     $         * (dist * invRpp) * EXP(-(dist * invRpp)**2)
      temp2 = FOUR / (SQRTPI * Rpp**3) * EXP(-(dist * invRpp)**2)
c     Include these prefactors into the inverse distance terms
      invdist5 = invdist5 * (temp1 * THREE - temp2 * dist3)
      invdist3 = invdist3 * temp1
c     Predetermine all possible values allowed by symmetry
      aa = r(1)**2 * invdist5 - invdist3
      bb = r(2)**2 * invdist5 - invdist3
      gg = r(3)**2 * invdist5 - invdist3
c
      ab = r(1) * r(2) * invdist5
      ag = r(1) * r(3) * invdist5
      bg = r(2) * r(3) * invdist5
c     Fill the temp array
      tr(1,1) = aa
      tr(2,1) = ab
      tr(3,1) = ag
      tr(1,2) = ab
      tr(2,2) = bb
      tr(3,2) = bg
      tr(1,3) = ag
      tr(2,3) = bg
      tr(3,3) = gg
      end function t2r

      function t2c(r, dist, invdist, a1, a2, pol1, pol2) result (tc)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
c
c     Input variables
      integer, intent(in) ::  a1, a2
      double precision, intent(in) :: dist, invdist
      double complex, intent(in) :: pol1, pol2
      double precision, intent(in) :: r(3)
c
c     Local variables
      double precision dist3
      double precision invdist3
      double precision invdist5
      double precision temp1, temp2
      double precision Rp1, Rp2
      double precision Rpp, invRpp
      double precision aa, ab, ag
      double precision bb, bg
      double precision gg
      double precision util_erf
      external util_erf
c
c     Returns
      double complex tc(3,3)
c 
c     Compute powers and inverses once
      dist3 = dist * dist * dist
      invdist3 = invdist  * invdist * invdist
      invdist5 = invdist3 * invdist * invdist
c
c     Calculate R values
      Rp1   = (SQRTTWOINVPI * (pol1 * THIRD))**THIRD
      Rp2   = (SQRTTWOINVPI * (pol2 * THIRD))**THIRD
      Rpp   = SQRT(Rp1**2 + Rp2**2)
      invRpp = ONE / Rpp
c
c     Calculate screening values
      temp1 = util_erf(dist * invRpp) - TWOINVSQRTPI
     $         * (dist * invRpp) * EXP(-(dist * invRpp)**2)
      temp2 = FOUR / (SQRTPI * Rpp**3) * EXP(-(dist * invRpp)**2)
c     Include these prefactors into the inverse distance terms
      invdist5 = invdist5 * (temp1 * THREE - temp2 * dist3)
      invdist3 = invdist3 * temp1
c     Predetermine all possible values allowed by symmetry
      aa = r(1)**2 * invdist5 - invdist3
      bb = r(2)**2 * invdist5 - invdist3
      gg = r(3)**2 * invdist5 - invdist3
c
      ab = r(1) * r(2) * invdist5
      ag = r(1) * r(3) * invdist5
      bg = r(2) * r(3) * invdist5
c     Fill the temp array
      tc(1,1) = aa
      tc(2,1) = ab
      tc(3,1) = ag
      tc(1,2) = ab
      tc(2,2) = bb
      tc(3,2) = bg
      tc(1,3) = ag
      tc(2,3) = bg
      tc(3,3) = gg
      end function t2c

