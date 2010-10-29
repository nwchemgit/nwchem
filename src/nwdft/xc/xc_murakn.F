      Subroutine xc_murakn(r,w,nr,range,lreset)

C$Id$
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
cedo      double precision alphad(36)
cedo      data alphad/ 5d0, 5d0,
cedoc Li-N
cedo     *       7d0,  7d0,  5d0,  5d0,  5d0,  5d0,  5d0,  5d0,
cedoC Na-Ar
cedo     *       7d0,  7d0,  5d0,  5d0,  5d0,  5d0,  5d0,  5d0,
cedoC K-Kr
cedo     *       7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,
cedoC Cu-Kr
cedo     *       7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0,  7d0/
c
c Definition of Euler-Maclaurin numerical quadrature points and weights
c for radial integrals.
c Transformation from 0<r<infty to 0<x<1 according to
c ME Mura and PJ Knowles, J Chem Phys 104, 9848 (1996)
c
c***************************************************************************
      if(lreset) then
         alpha=-range/log(1d0-((1d0+nr)/(2d0+nr))**3)
      else
c
c     huub recipe
c
         alpha=3.3d0*range
c         alpha=5d0
c         if(znumber.lt.37) alpha=alphad(znumber)
      endif
      fmn = rm/(1d0+nr)

      do 10 i = 1,nr

       qi = dble(i)/(nr+1d0)
       ri = -alpha*log(1.D0 - qi**rm)
       wi = fmn*alpha*(ri*ri)/(1.D0 - qi**rm)*qi**(rm-1.d0)

       r(i) = ri
       w(i) = wi

   10 continue

      end
