C$Id$
************************************************************************
*                                                                      *
      subroutine ecp_s2ccoef (n,Xc2s,Xs2c,tmp,nc,df)
*                                                                      *
*   Set up inverse transformation, from sphericals to cartesians       *
*                                                                      *
*   n (inp) - angular momentum of s.s.h.                               *
*   Xc2s (inp) - coefficents of cartesian to spherical transformation  *
*   Xs2c (out) - coefficents of spherical to cartesian transformation  *
*   tmp (scr) - array for storing cartesian overlap matrix             *
*   nc (inp) - number of cartesians, (l+1)*(l+2)/2                     *
*   df (inp) - array of double factorials                              *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer i,ia,ib,ii,ja,jb,jj,ka,kb,kk,n,nc
      double precision d,one,zero
      double precision Xc2s(nc,nc),Xs2c(nc,nc),tmp(nc*nc),df(0:n+1)
      parameter (zero = 0.0D0, one = 1.0D0)
*
      i = 0
      d = df(n+1)
      do ia = n,0,-1
        do ja = n-ia,0,-1
          ka = n-ia-ja
          do ib = n,0,-1
            ii = ia+ib
            do jb = n-ib,0,-1
              jj = ja+jb
              kb = n-ib-jb
              kk = ka+kb
              i = i+1
              if ((mod(ii,2) .eq. 0) .and.
     &            (mod(jj,2) .eq. 0) .and.
     &            (mod(kk,2) .eq. 0)) then
                tmp(i) = df(ii/2)*df(jj/2)*df(kk/2)/d
              else
                tmp(i) = zero
              end if
            end do
          end do
        end do
      end do
      call dgemm ('T','N',nc,nc,nc,one,Xc2s,nc,tmp,nc,zero,Xs2c,nc)
*
      return
      end
