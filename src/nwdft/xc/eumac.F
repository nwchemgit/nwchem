      Subroutine eumac(alpha,r,w,nr,rm)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Dimension r(*),w(*)
c
c Definition of Euler-Maclaurin numerical quadrature points and weights
c for radial integrals.
c
c***************************************************************************
      nr = nr+1

      fmn = dble(rm)/dble(nr)

      do 10 i = 1,nr-1

       qi = dble(i)/dble(nr)
       ri = alpha*(qi/(1.D0 - qi))**rm
       wi = fmn*(ri*ri*ri)/(qi*(1.D0 - qi))

       r(i) = ri
       w(i) = wi

   10 continue

      nr = nr-1

      end
