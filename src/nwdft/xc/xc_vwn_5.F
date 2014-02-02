c     Correlation energy and potential as determined via the local 
c     density functional of Vosko, Wilk and Nusair (Can. J. of Phys. 58, 
c     1200 (1980)), functional V.
*
* $Id$
*
#define FNL_VWN_5
c xc_vwn_a contains the restricted calculation stuff.
#include "xc_vwn_a.fh"
c
c     Compute the function deltaEc(rs,zeta) function and its derivatives
c     wrt rs and zeta for the spin-unrestricted case - the rest has the
c     same form for all VWN functionals and is handled in the header files.
c
            dec_rsz = devwn_rsz
            d1dec_rsz(1) = d1devwn_rsz(1)
            d1dec_rsz(2) = d1devwn_rsz(2)
#ifdef SECOND_DERIV
            d2dec_rsz(1) = d2devwn_rsz(1)
            d2dec_rsz(2) = d2devwn_rsz(2)
            d2dec_rsz(3) = d2devwn_rsz(3)
#endif
#ifdef THIRD_DERIV
            d3dec_rsz(1) = d3devwn_rsz(1)
            d3dec_rsz(2) = d3devwn_rsz(2)
            d3dec_rsz(3) = d3devwn_rsz(3)
            d3dec_rsz(4) = d3devwn_rsz(4)
#endif
c
c xc_vwn_a contains the unrestricted calculation stuff.
#include "xc_vwn_b.fh"
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_vwn_5.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_vwn_5.F"
#endif
