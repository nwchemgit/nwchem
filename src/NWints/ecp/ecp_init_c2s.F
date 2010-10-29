C$Id$
************************************************************************
*                                                                      *
      subroutine ecp_init_c2s (l_max,tco,lstart,ltr,inorm,itrans,
     &    DryRun,mem_max)
*                                                                      *
*   Initialize transformation from cartesian monomials x^i y^j z^k to  *
*   real solid spherical tensor functions r^n G_lm(theta,phi).         *
*                                                                      *
*   l_max (inp) - maximum angular momentum required for transformation *
*   tco (out) - array of cartesian to spherical transformations        *
*   lstart (out) - integer array of pointers to start of block for     *
*                  each l value in Xc2s and Xs2c.                      *
*   ltr (inp) - length of transformation array                         *
*   inorm (inp) - 0 for unnormalized sphericals, 1 for spherical       *
*                 tensors, 2 for spherical harmonics                   *
*   itrans (inp) - 0 for c2s only, 1 for both c2s and s2c              *
*   DryRun (inp) - logical for dry run, i.e. evaluate memory required  *
*   mem_max (out) - memory required by this routine                    *
*                                                                      *
*   Notes:                                                             *
*                                                                      *
*   1. The array tco and the array lstart are in fact the same array   *
*   passed twice so that it can be used as both integer and real.      *
*   2. The array tco is used as scratch as well as for final storage   *
*   of the transformation array.                                       *
*                                                                      *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "errquit.fh"
      integer i,inorm,itrans,i_df,i_free,i_m,i_o,i_p,i_scr,l,l_max,ltr,
     &    mem_max,ncm,nco,ncp,lstart(0:1,0:l_max)
      double precision zero,one,wa
      double precision tco(ltr)
      logical DryRun
      parameter (zero = 0.0d00, one = 1.0d00)
*
      if ((itrans .gt. 1) .or. (itrans .lt. 0))
     &    call errquit (' Invalid parameter itrans in ecp_init_c2s',99,
     &       BASIS_ERR)
*
*      Set up pointers to array tco
*
      if (DryRun) mem_max = 2*(l_max+1)
      i_o = 2*l_max+3
      nco = 0
      do l = 0,l_max
        nco = nco+l+1
        do i = 0,itrans
          if (.not.DryRun) lstart(i,l) = i_o
          i_o = i_o+nco*nco
        end do
      end do
      i_free = i_o

      i_df  = i_o               ! take care of compiler warnings
      i_scr = i_o+l_max+2

      if (itrans .eq. 1) then
        i_df = i_o
        i_scr = i_o+l_max+2
        i_free = i_scr+nco*nco
      end if
      if (DryRun) then
        mem_max = i_free-1
        return
      else
        if (ltr .lt. i_free-1) call errquit (
     &    ' Insufficient memory in ecp_init_c2s',99, MEM_ERR)
      end if
*
*   Generate unnormalized transformation coefficients by recursion
*
      do i = lstart(0,0),i_free-1
        tco(i) = zero
      end do
      i_p = lstart(0,0)
      i_o = i_p
      nco = 1
      ncp = 1
      tco(i_p) = one
      do l = 0,l_max-1
        i_m = i_o
        i_o = i_p
        i_p = lstart(0,l+1)
        ncm = nco
        nco = ncp
        ncp = nco+l+2
        call ecp_c2srecur(l,tco(i_p),tco(i_o),tco(i_m),ncp,nco,ncm)
      end do
*
*   Normalize coefficients (if required)
*
      nco = 1
      do l = 0,l_max
        i_o = lstart(0,l)
        call ecp_c2snorm(l,tco(i_o),nco,nco,inorm)
        nco = nco+l+2
      end do
*
*   Generate coefficients for lower angular momentum for given n.
*
      ncm = 1
      ncp = 6
      do l = 2,l_max
        i_p = lstart(0,l)+ncp*(2*l+1)
        i_m = lstart(0,l-2)
        call ecp_c2slower(l,tco(i_p),tco(i_m),ncp,ncm)
        ncm = ncm+l
        ncp = ncp+l+2
      end do
      if (itrans .eq. 0) return
*
*   Generate inverse transformation
*
      i_o = 1
      nco = 1
      tco(i_df) = 1
      do l = 0,l_max
        i_o = lstart(0,l)
        wa = 2*l+1
        tco(i_df+l+1) = tco(i_df+l)*wa
        i_p = lstart(1,l)
        call ecp_s2ccoef (l,tco(i_o),tco(i_p),tco(i_scr),nco,tco(i_df))
        i_m = i_p
        do i = l,0,-2
          call ecp_s2cnorm(i,tco(i_m),nco,nco,inorm)
          i_m = i_m+2*i+1
        end do
        nco = nco+l+2
      end do
*
      return
      end
