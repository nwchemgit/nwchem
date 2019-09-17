#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_vwn1_rpa.F
C> The Vosko, Wilk, Nusair correlation functional #1 with RPA parametrization
C>
C> @}
#endif
#endif
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The homogenous electron gas correlation functional of
C> Vosko, Wilk and Nusair #1 of 1980 with RPA parametrization
C>
C> The correlation functional #1 of Vosko, Wilk and Nusair parametrized
C> against RPA data [1].
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
#define FNL_VWN_1_RPA
#include "nwxc_c_vwn_a.fh"
c
c     Compute the function deltaEc(rs,zeta) function and its derivatives
c     wrt rs and zeta for the spin-unrestricted case - the rest has the
c     same form for all VWN functionals and is handled in the header files.
c
            dec_rsz = dec_rs1*fz
c           d1dec_rsz(1) = d1dec_rs1*fz
c           d1dec_rsz(2) = dec_rs1*d1fz
#ifdef SECOND_DERIV
c           d2dec_rsz(1) = d2dec_rs1*fz
c           d2dec_rsz(2) = d1dec_rs1*d1fz
c           d2dec_rsz(3) = dec_rs1*d2fz
#endif
#ifdef THIRD_DERIV
c           d3dec_rsz(1) = d3dec_rs1*fz
c           d3dec_rsz(2) = d2dec_rs1*d1fz
c           d3dec_rsz(3) = d1dec_rs1*d2fz
c           d3dec_rsz(4) = dec_rs1*d3fz
#endif
c
#include "nwxc_c_vwn_b.fh"
c
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_vwn1_rpa.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_vwn1_rpa.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_vwn1_rpa.F"
#endif
#undef NWAD_PRINT
C>
C> @}
