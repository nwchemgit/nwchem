C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_cart_xpd (n,n_n,n_all,X,Y,Z,xn,yn,zn,tmp,coef,
     &    sphcart,csco,lcsco)
*                                                                      *
*   Expand cartesian function on one centre about second centre.       *
*                                                                      *
*   n (inp) - angular momentum of original cartesian                   *
*   n_n (inp) - number of components of original cartesian             *
*               = (n+1)*(n+2)/2                                        *
*   n_all (inp) - number of components of cartesians on new centre     *
*           = sum_m^n (m+1)*(m+2)/2 = (n+1)*(n+2)*(n+3)/6              *
*   X,Y,Z (inp) - relative cartesian coordinates of new centre         *
*   xn,yn,zn (scr) - work arrays to store powers of X,Y,Z              *
*   tmp (scr) - work array for spherical transformation                *
*   coef (out) - array of coefficients of the cartesians in terms of   *
*                sphericals on the new centre                          *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "errquit.fh"
      integer n,n_n,n_all
      integer lcsco,sphcart
      integer a,b,c,i,j,k,ii,jj,kk,m,ind_m,ind_n,ind_s,l_comp
      integer binomial(231)
      double precision X,Y,Z,xn(0:n),yn(0:n),zn(0:n),
     &    tmp(n_n,n_n),coef(n_n,n_all)
      double precision ifac,jfac,kfac,zero
      double precision csco(lcsco)
      parameter (zero = 0.0d00)
      data binomial/1, 1,1, 1,2,1, 1,3,3,1, 1,4,6,4,1, 1,5,10,10,5,1,
     &    1,6,15,20,15,6,1, 1,7,21,35,35,21,7,1, 1,8,28,56,70,56,28,8,1,
     &    1,9,36,84,126,126,84,36,9,1, 1,10,45,120,210,252,210,120,45,
     &    10,1, 1,11,55,165,330,462,462,330,165,55,11,1, 1,12,66,220,
     &    495,792,924,792,495,220,66,12,1, 1,13,78,286,715,1287,1716,
     &    1716,1287,715,286,78,13,1, 1,14,91,364,1001,2002,3003,3432,
     &    3003,2002,1001,364,91,14,1, 1,15,105,455,1365,3003,5005,6435,
     &    6435,5005,3003,1365,455,105,15,1, 1,16,120,560,1820,4368,8008,
     &    11440,12870,11440,8008,4368,1820,560,120,16,1, 1,17,136,680,
     &    2380,6188,12376,19448,24310,24310,19448,12376,6188,2380,680,
     &    136,17,1, 1,18,153,816,3060,8568,18564,31824,43758,48620,
     &    43758,31824,18564,8568,3060,816,153,18,1, 1,19,171,969,3876, 
     &    11628,27132,50388,75582,92378,92378,75582,50388,27132,11628,
     &    3876,969,171,19,1, 1,20,190,1140,4845,15504,38760,77520,
     &    125970,167960,184756,167960,125970,77520,38760,15504,4845,
     &    1140,190,20,1/
*
      if (n .gt. 20) call errquit (
     &    'Cannot handle angular momentum > 20 in ecp_cart_xpd',99,
     &       UNKNOWN_ERR)
*
*   Set up monomials in X, Y and Z
*
      xn(0) = 1.0d00
      yn(0) = 1.0d00
      zn(0) = 1.0d00
      do i = 1,n
        xn(i) = xn(i-1)*X
        yn(i) = yn(i-1)*Y
        zn(i) = zn(i-1)*Z
      end do
*
*     Loop over angular momenta on new centre
*
      ind_s = 1
      do m = 0,n
*
*   Loop over components of original cartesian
*
        ind_n = 0
        do i = n,0,-1
          ii = i*(i+1)/2+1
          do j = n-i,0,-1
            jj = j*(j+1)/2+1
            k = n-i-j
            kk = k*(k+1)/2+1
            ind_n = ind_n+1
*
*     Loop over components on new centre
*
            ind_m = 0
            do a = m,0,-1
              if (i .ge. a) then
                ifac = binomial(ii+a)
                ifac = ifac*xn(i-a)
                do b = m-a,0,-1
                  c = m-a-b
                  ind_m = ind_m+1
                  tmp(ind_n,ind_m) = 0d0
                  if (j .ge. b.and.k.ge.c) then
                    jfac = binomial(jj+b)
                    jfac = jfac*yn(j-b)
                    kfac = binomial(kk+c)
                    kfac = kfac*zn(k-c)
                    tmp(ind_n,ind_m) = ifac*jfac*kfac
                  end if
                end do
              else
                 call dcopy(m-a+1,zero,0,tmp(ind_n,ind_m+1),n_n)
                ind_m = ind_m+m-a+1
              end if
            end do
*
          end do
        end do
*
*   Now transform to spherical harmonics. This must be done with the 
*   inverse transformation because what we are really doing is inserting 
*   the resolution of the identity into the cartesian expansion; the 
*   cartesians on the new centre get the "real" transformation, the 
*   coefficients get the transpose of the inverse transformation.
*
        if (sphcart .eq. 0) then
          call dcopy (n_n*ind_m,tmp,1,coef(1,ind_s),1)
        else
          call ecp_cstrans (m,ind_m,ind_n,m,mod(m,2),l_comp,tmp,n_n,
     &        coef(1,ind_s),n_n,csco,lcsco,csco,1,1,-1)
        end if
        ind_s = ind_s+ind_m
*
      end do
*
      return
      end
