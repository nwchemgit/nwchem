C$Id$
************************************************************************
*                                                                      *
      subroutine ecp_c2snorm (l,X,ldX,nc,inorm)
*                                                                      *
*   Normalize solid spherical harmonic transformation coefficients     *
*                                                                      *
*   l (inp) - angular momentum of s.s.h.                               *
*   X (i/o) - solid spherical harmonic transformation coefficients     *
*   nc (inp) - number of cartesians, (l+1)*(l+2)/2                     *
*   inorm (inp) - flag for type of normalization. 0 = none,            *
*                 1 = spherical tensor, 2 = spherical harmonic         *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "errquit.fh"
      integer i,inorm,l,m,ldX,nc
      double precision wa,wb,wc,zero,half
      double precision X(ldX,-l:l)
      parameter (zero = 0.0D0, half = 0.5D0)
*
      if (inorm .eq. 0) return
      if (inorm .eq. 1) then
        wa = 1.0d00
      else if (inorm .eq. 2) then
        wb = 2*l+1
        wc = acos(zero)
        wa = half*sqrt(wb/wc)
        wb = wa*sqrt(half)
        do i = 1,nc
          X(i,0) = X(i,0)*wb
        end do
      else
        call errquit (' Illegal parameter inorm in ecp_c2snorm',99,
     &       UNKNOWN_ERR)
      end if
      do m = 1,l
        wb = l-m+1
        wc = l+m
        wa = wa/sqrt(wb*wc)
        do i = 1,nc
          X(i,m) = X(i,m)*wa
          X(i,-m) = X(i,-m)*wa
        end do
      end do
*
      return
      end
