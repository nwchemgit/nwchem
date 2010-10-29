C$Id$
************************************************************************
*                                                                      *
      subroutine ecp_s2cnorm (l,X,ldX,nc,inorm)
*                                                                      *
*   Normalize transformation from solid spherical harmonics to         *
*   cartesian monomials, by multiplying by inverse of overlap.         *
*                                                                      *
*   l (inp) - angular momentum of s.s.h.                               *
*   X (i/o) -  transformation coefficients                             *
*   nc (inp) - number of cartesians, (l+1)*(l+2)/2                     *
*   inorm (inp) - flag for type of normalization. 0 = none,            *
*                 1 = spherical tensor, 2 = spherical harmonic         *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "errquit.fh"
      integer i,inorm,k,l,ldX,m,n,nc
      double precision wa,wb,wc,zero
      double precision X(ldX,nc)
      parameter (zero = 0.0D0)
*
      i = inorm-1
      if (abs(i) .gt. 1) 
     &    call errquit (' Illegal parameter inorm in ecp_s2cnorm',99,
     &       BASIS_ERR)
      n = 2*l+1
      if (inorm .eq. 2) then
        wa = acos(zero)*8
        do m = 1,n
          do i = 1,nc
            X(m,i) = X(m,i)*wa
          end do
        end do
      else
        wb = n
        k = l+1
        do i = 1,nc
          X(k,i) = X(k,i)*wb
        end do
        wa = 2*n
        do m = 1,l
          if (inorm .eq. 0) then
            wb = l-m+1
            wc = l+m
            wa = wa/(wb*wc)
          end if
          do i = 1,nc
            X(k-m,i) = X(k-m,i)*wa
            X(k+m,i) = X(k+m,i)*wa
          end do
        end do
      end if
*
      return
      end
