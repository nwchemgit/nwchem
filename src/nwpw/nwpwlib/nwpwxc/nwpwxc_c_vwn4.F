#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_vwn4.F
C> The Vosko, Wilk, Nusair correlation functional #4
C>
C> @}
#endif
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The homogenous electron gas correlation functional of
C> Vosko, Wilk and Nusair #4 of 1980
C>
C> The correlation functional #4 of Vosko, Wilk and Nusair parametrized
C> against Ceperley-Alder Monte Carlo data [1].
C>
C> [1] S.H. Vosko, L. Wilk, M. Nusair,
C>     <i>"Accurate spin-dependent electron liquid correlation energies
C>     for local spin density calculations: a critical analysis"</i>,
C>     Can. J. Phys. <b>58</b>, 1200-1211 (1980), DOI:
C>     <a href="https://doi.org/10.1139/p80-159">10.1139/p80-159</a>.
C>
*
* $Id$
*
#define FNL_VWN_4
#include "nwpwxc_c_vwn_a.fh"
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
c
#include "nwpwxc_c_vwn_b.fh"
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwpwxc_c_vwn4.F"
#endif
C>
C> @}
