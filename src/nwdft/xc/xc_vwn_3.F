c     Correlation energy and potential as determined via the local 
c     density functional of Vosko, Wilk and Nusair (Can. J. of Phys. 58, 
c     1200 (1980)), functional III.
*
* $Id$
*
#define FNL_VWN_3
#include "xc_vwn_a.fh"
c
c     Compute the function deltaEc(rs,zeta) function and its derivatives
c     wrt rs and zeta for the spin-unrestricted case - the rest has the
c     same form for all VWN functionals and is handled in the header files.
c
            dec_rsz = (dec_rs1/derpa_rs1)*devwn_rsz
            d1dec_rsz(1) = (d1dec_rs1/derpa_rs1)*devwn_rsz
     &           - (dec_rs1*d1derpa_rs1/derpa_rs1**2)*devwn_rsz
     &           + (dec_rs1/derpa_rs1)*d1devwn_rsz(1)
            d1dec_rsz(2) = (dec_rs1/derpa_rs1)*d1devwn_rsz(2)
#ifdef SECOND_DERIV
            d2dec_rsz(1) = (dec_rs1/derpa_rs1)*d2devwn_rsz(1)
     &           + (d2dec_rs1/derpa_rs1)*devwn_rsz
     &           - (dec_rs1*d2derpa_rs1/derpa_rs1**2)*devwn_rsz
     &           + 2.d0*((d1dec_rs1/derpa_rs1)*d1devwn_rsz(1)
     &           - (d1dec_rs1*d1derpa_rs1/derpa_rs1**2)*devwn_rsz
     &           - (dec_rs1*d1derpa_rs1/derpa_rs1**2)*d1devwn_rsz(1)
     &           + (dec_rs1*d1derpa_rs1**2/derpa_rs1**3)*devwn_rsz)
            d2dec_rsz(2) = (dec_rs1/derpa_rs1)*d2devwn_rsz(2)
     &           + ((d1dec_rs1/derpa_rs1)
     &           - (dec_rs1*d1derpa_rs1/derpa_rs1**2))*d1devwn_rsz(2)
            d2dec_rsz(3) = (dec_rs1/derpa_rs1)*d2devwn_rsz(3)
#endif
#ifdef THIRD_DERIV
            d3dec_rsz(1) = (dec_rs1/derpa_rs1)*d3devwn_rsz(1)
     &           + 3.0d0*(derpa_rs1*d1dec_rs1
     &              - dec_rs1*d1derpa_rs1)*d2devwn_rsz(1)/derpa_rs1**2
     &           + 3.0d0*(derpa_rs1**2*d2dec_rs1 
     &             - derpa_rs1*(2.0d0*d1dec_rs1*d1derpa_rs1
     &              +dec_rs1*d2derpa_rs1) 
     &              +2.0d0*dec_rs1*d1derpa_rs1**2)
     &              *d1devwn_rsz(1)/derpa_rs1**3
     &           + (d3dec_rs1*derpa_rs1**3 
     &              + 6.0d0*derpa_rs1*d1derpa_rs1
     &               *(d1dec_rs1*d1derpa_rs1+dec_rs1*d2derpa_rs1)
     &              - derpa_rs1**2*(3.0d0*d2dec_rs1*d1derpa_rs1
     &               + 3.0d0*d1dec_rs1*d2derpa_rs1 
     &               + dec_rs1*d3derpa_rs1)
     &              - 6.0d0*dec_rs1*d1derpa_rs1**3)
     &              *devwn_rsz/derpa_rs1**4
            d3dec_rsz(2) = (dec_rs1/derpa_rs1)*d3devwn_rsz(2)
     &           + 2.0d0*(derpa_rs1*d1dec_rs1
     &              - dec_rs1*d1derpa_rs1)*d2devwn_rsz(2)/derpa_rs1**2
     &           + (derpa_rs1**2*d2dec_rs1 
     &             - derpa_rs1*(2.0d0*d1dec_rs1*d1derpa_rs1
     &              +dec_rs1*d2derpa_rs1) 
     &              +2.0d0*dec_rs1*d1derpa_rs1**2)
     &              *d1devwn_rsz(2)/derpa_rs1**3
            d3dec_rsz(3) = (dec_rs1/derpa_rs1)*d3devwn_rsz(3)
     &           + (derpa_rs1*d1dec_rs1
     &              - dec_rs1*d1derpa_rs1)*d2devwn_rsz(3)/derpa_rs1**2
            d3dec_rsz(4) = (dec_rs1/derpa_rs1)*d3devwn_rsz(4)
#endif
c
#include "xc_vwn_b.fh"
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_vwn_3.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_vwn_3.F"
#endif
