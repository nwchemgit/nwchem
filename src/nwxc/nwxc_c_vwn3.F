#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_vwn3.F
C> The Vosko, Wilk, Nusair correlation functional #3
C>
C> @}
#endif
#endif
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The homogenous electron gas correlation functional of
C> Vosko, Wilk and Nusair #3 of 1980
C>
C> The correlation functional #3 of Vosko, Wilk and Nusair parametrized
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
#define FNL_VWN_3
#include "nwxc_c_vwn_a.fh"
c
c     Compute the function deltaEc(rs,zeta) function and its derivatives
c     wrt rs and zeta for the spin-unrestricted case - the rest has the
c     same form for all VWN functionals and is handled in the header files.
c
            dec_rsz = (dec_rs1/derpa_rs1)*devwn_rsz
c           d1dec_rsz(1) = (d1dec_rs1/derpa_rs1)*devwn_rsz
c    &           - (dec_rs1*d1derpa_rs1/derpa_rs1**2)*devwn_rsz
c    &           + (dec_rs1/derpa_rs1)*d1devwn_rsz(1)
c           d1dec_rsz(2) = (dec_rs1/derpa_rs1)*d1devwn_rsz(2)
#ifdef SECOND_DERIV
c           d2dec_rsz(1) = (dec_rs1/derpa_rs1)*d2devwn_rsz(1)
c    &           + (d2dec_rs1/derpa_rs1)*devwn_rsz
c    &           - (dec_rs1*d2derpa_rs1/derpa_rs1**2)*devwn_rsz
c    &           + 2.d0*((d1dec_rs1/derpa_rs1)*d1devwn_rsz(1)
c    &           - (d1dec_rs1*d1derpa_rs1/derpa_rs1**2)*devwn_rsz
c    &           - (dec_rs1*d1derpa_rs1/derpa_rs1**2)*d1devwn_rsz(1)
c    &           + (dec_rs1*d1derpa_rs1**2/derpa_rs1**3)*devwn_rsz)
c           d2dec_rsz(2) = (dec_rs1/derpa_rs1)*d2devwn_rsz(2)
c    &           + ((d1dec_rs1/derpa_rs1)
c    &           - (dec_rs1*d1derpa_rs1/derpa_rs1**2))*d1devwn_rsz(2)
c           d2dec_rsz(3) = (dec_rs1/derpa_rs1)*d2devwn_rsz(3)
#endif
#ifdef THIRD_DERIV
c           d3dec_rsz(1) = (dec_rs1/derpa_rs1)*d3devwn_rsz(1)
c    &           + 3.0d0*(derpa_rs1*d1dec_rs1
c    &              - dec_rs1*d1derpa_rs1)*d2devwn_rsz(1)/derpa_rs1**2
c    &           + 3.0d0*(derpa_rs1**2*d2dec_rs1 
c    &             - derpa_rs1*(2.0d0*d1dec_rs1*d1derpa_rs1
c    &              +dec_rs1*d2derpa_rs1) 
c    &              +2.0d0*dec_rs1*d1derpa_rs1**2)
c    &              *d1devwn_rsz(1)/derpa_rs1**3
c    &           + (d3dec_rs1*derpa_rs1**3 
c    &              + 6.0d0*derpa_rs1*d1derpa_rs1
c    &               *(d1dec_rs1*d1derpa_rs1+dec_rs1*d2derpa_rs1)
c    &              - derpa_rs1**2*(3.0d0*d2dec_rs1*d1derpa_rs1
c    &               + 3.0d0*d1dec_rs1*d2derpa_rs1 
c    &               + dec_rs1*d3derpa_rs1)
c    &              - 6.0d0*dec_rs1*d1derpa_rs1**3)
c    &              *devwn_rsz/derpa_rs1**4
c           d3dec_rsz(2) = (dec_rs1/derpa_rs1)*d3devwn_rsz(2)
c    &           + 2.0d0*(derpa_rs1*d1dec_rs1
c    &              - dec_rs1*d1derpa_rs1)*d2devwn_rsz(2)/derpa_rs1**2
c    &           + (derpa_rs1**2*d2dec_rs1 
c    &             - derpa_rs1*(2.0d0*d1dec_rs1*d1derpa_rs1
c    &              +dec_rs1*d2derpa_rs1) 
c    &              +2.0d0*dec_rs1*d1derpa_rs1**2)
c    &              *d1devwn_rsz(2)/derpa_rs1**3
c           d3dec_rsz(3) = (dec_rs1/derpa_rs1)*d3devwn_rsz(3)
c    &           + (derpa_rs1*d1dec_rs1
c    &              - dec_rs1*d1derpa_rs1)*d2devwn_rsz(3)/derpa_rs1**2
c           d3dec_rsz(4) = (dec_rs1/derpa_rs1)*d3devwn_rsz(4)
#endif
c
#include "nwxc_c_vwn_b.fh"
c
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_vwn3.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_vwn3.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_vwn3.F"
#endif
#undef NWAD_PRINT
C>
C> @}
