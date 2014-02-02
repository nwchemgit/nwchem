c     Correlation energy and potential as determined via the local 
c     density functional of Vosko, Wilk and Nusair (Can. J. of Phys. 58, 
c     1200 (1980)), functional I, the Monte Carlo parameterization.
*
* $Id$
*
#define FNL_VWN_1
#include "xc_vwn_a.fh"
c
c     Compute the function deltaEc(rs,zeta) function and its derivatives
c     wrt rs and zeta for the spin-unrestricted case - the rest has the
c     same form for all VWN functionals and is handled in the header files.
c
            dec_rsz = dec_rs1*fz
            d1dec_rsz(1) = d1dec_rs1*fz
            d1dec_rsz(2) = dec_rs1*d1fz
#ifdef SECOND_DERIV
            d2dec_rsz(1) = d2dec_rs1*fz
            d2dec_rsz(2) = d1dec_rs1*d1fz
            d2dec_rsz(3) = dec_rs1*d2fz
#endif
#ifdef THIRD_DERIV
            d3dec_rsz(1) = d3dec_rs1*fz
            d3dec_rsz(2) = d2dec_rs1*d1fz
            d3dec_rsz(3) = d1dec_rs1*d2fz
            d3dec_rsz(4) = dec_rs1*d3fz
#endif
c
#include "xc_vwn_b.fh"
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_vwn_1.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_vwn_1.F"
#endif
