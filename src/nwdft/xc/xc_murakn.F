      Subroutine xc_murakn(r,w,nr,range,lreset)

C$Id: xc_murakn.F,v 1.9 2000-02-02 02:37:40 edo Exp $
      implicit none

      double precision r(*) ! grid pts coord [output]
      double precision w(*) ! grid pts weights [output]
c
      double precision rm    
      parameter (rm=3d0)
c
      integer nr,i
      double precision fmn,qi,ri,wi,alpha
      double precision range ! max extent of basis f
c
      logical lreset
      double precision alphad(36)
      data alphad/ 5d0, 5d0,
c Li-N
     *       7d0,  7d0,  5d0,  5d0,  5d0,  5d0,  5d0,  5d0,
C Na-Ar
     *       7d0,  7d0,  5d0,  5d0,  5d0,  5d0,  5d0,  5d0,
C K-Kr
     *       7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,
C Cu-Kr
     *       7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0/
c
c Definition of Euler-Maclaurin numerical quadrature points and weights
c for radial integrals.
c Transformation from 0<r<infty to 0<x<1 according to
c ME Mura and PJ Knowles, J Chem Phys 104, 9848 (1996)
c
c***************************************************************************
      nr = nr+1
      if(lreset) then
         alpha=-range/log(1d0-(nr/(1d0+nr))**3)
      else
c
c     huub recipe
c
         alpha=3.3d0*range
c         alpha=5d0
c         if(znumber.lt.37) alpha=alphad(znumber)
      endif
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
