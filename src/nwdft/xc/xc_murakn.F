      Subroutine xc_murakn(alpha,r,w,nr,rm)

C$Id: xc_murakn.F,v 1.1 1999-01-12 22:24:45 edo Exp $
      implicit none

      double precision r(*) ! grid pts coord [output]
      double precision w(*) ! grid pts weights [output]
      double precision alpha !   [input]
      double precision rm    !   [input]
c
      integer nr,i
      double precision fmn,qi,ri,wi
c
      double precision alpha
c
c Definition of Euler-Maclaurin numerical quadrature points and weights
c for radial integrals.
c Transformation from 0<r<intfy to 0<x<1 according to
c ME Mura and PJ Knowles, J Chem Phys 104, 9848 (1996)
c
c***************************************************************************
      nr = nr+1
      alpha=5.d0
      fmn = rm/dble(nr)

      do 10 i = 1,nr-1

       qi = dble(i)/dble(nr)
       ri = -alpha*log(1.D0 - qi**rm)
       wi = fmn*alpha*(ri*ri)/(1.D0 - qi**rm)*qi**(rm-1.d0)

       r(i) = ri
       w(i) = wi

   10 continue

      nr = nr-1

      end
