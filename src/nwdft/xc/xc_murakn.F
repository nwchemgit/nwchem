      Subroutine xc_murakn(znumber,r,w,nr,range,lreset)

C$Id: xc_murakn.F,v 1.5 2000-01-06 00:41:35 edo Exp $
      implicit none

      double precision r(*) ! grid pts coord [output]
      double precision w(*) ! grid pts weights [output]
      integer znumber !   atomic number [input]
c
      double precision rm    
      parameter (rm=3d0)
c
      integer nr,i
      double precision fmn,qi,ri,wi,alpha
      double precision range ! max extent of basis f
c
      logical lreset
c
c Definition of Euler-Maclaurin numerical quadrature points and weights
c for radial integrals.
c Transformation from 0<r<infty to 0<x<1 according to
c ME Mura and PJ Knowles, J Chem Phys 104, 9848 (1996)
c
c***************************************************************************
      nr = nr+1
cedo      if(lreset) then
cedo         alpha=-range/log(1d0-(nr/(1d0+nr))**3)
cedo      else
         alpha=5d0
cedo      endif
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
