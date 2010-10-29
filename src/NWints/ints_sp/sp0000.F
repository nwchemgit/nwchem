      subroutine sp0000(gout)
c $Id$
c     *****  special fast routine for -p- loop for 0000 ****
      implicit none
c     
      double precision a0, b0, c0, abc
      common/tabint/
     +     a0(333),b0(333),c0(333),abc(5001)
c     
#include "auxvar.fh"
#include "miscg.fh"
#include "ginf.fh"
#include "pgeom.fh"
#include "shllfo.fh"
#include "spgeom.fh"
#include "qgeom.fh"
#include "maxc.fh"
c     
      double precision gout(*)
      double precision sixty
      double precision qq, theta
      integer n
      double precision cq
      integer k, l, i
      double precision g0000, csck, gcrcds, ecd
      double precision p, pp, ppp, qqtest, h0000, gdecd, f0
      double precision theta2, theta3, theta4, aqx, aqz,qperp2
      data sixty/60.0d0/
c     
      g0000 = 0.d0
      do 220 k = 1,ngc
         gc = cg(k)
         csck = csc(k)
         gcrcds = gc*rcdsq
         do 210 l = 1,ngd
            gd = dg(l)
            gcd = gc+gd
            ecd = 1.d0/gcd
            gdecd = gd*ecd
            ppp = -gdecd*gcrcds
            if(ppp+sixty)480,500,500
 480        pp = 0.0d0
            go to 520
 500        pp = ecd* dexp(ppp)
 520        qqtest = pp*cmaxc(k)*cmaxd(l)
            if (qqtest .le. error1) go to 100
            ismlq = 0
            go to 120
 100        if (qqtest .le. error2) go to 210
            ismlq = 1
 120        cq = gdecd*rcd
            aqx = acx+sing*cq
            aqz = acz+cosg*cq
            qperp2 = aqx*aqx+acy2
            h0000 = 0.d0
            do 200 i = 1,ngangb
               isml = ismlq+ismlp(i)
               if (isml .ge. 2) go to 200
               auxvar = var(isml+1)
               p = ((aqz-app(i))**2+qperp2)/(ep(i)+ecd)
               if (p .le. auxvar) go to 180
               h0000 = h0000+dp00p(i)*
     $              dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
               go to 200
 180           continue
               qq = p*12.5d0
               n =  idint(qq)
               theta = qq- dble(n)
               theta2 = theta*(theta-1.d0)
               theta3 = theta2*(theta-2.d0)
               theta4 = theta2*(theta+1.d0)
               f0 = a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2)
               h0000 = h0000+ dp00p(i)*f0/dsqrt(gp(i)+gcd)
 200        continue
            g0000 = g0000+h0000*csck*csd(l)*pp
 210     continue
 220  continue
      gout(1) = g0000
      return
      end
