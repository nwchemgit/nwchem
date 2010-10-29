      subroutine sp0001(gout)
c $Id$
c        *****  special fast routine for -p- loop for 0001 ****
      implicit none
      double precision gout(*)
c
      double precision a0, b0, c0, a1, b1, c1, a2
      common/tabint/a0(333),b0(333),c0(334),a1(333),b1(333),c1(334),
     +  a2(4000)
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
      double precision qq, theta, cq, dq
      integer n
      integer k, l, i
      double precision sixty, tenm12, gout1, gout2, gout3, gout4
      double precision q0, q1, ecd, qqq, qqtest
      double precision dq00, dq01, dq10, dq11, v, aqx, aqz
      double precision qperp, qperp2, cosp, sinp, p, q, pqab
      double precision h0000, h0001, h0003, g0001, g0002, g0003
      double precision theta2, theta3, theta4, u, g, t1, t2, t3
c
      data sixty,tenm12/60.0d0,1.0d-12/
c
      gout1 = 0.0d0
      gout2 = 0.0d0
      gout3 = 0.0d0
      gout4 = 0.0d0
      do 940 k = 1,ngc
      gc = cg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = 1.0d0/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 480,500,500
  480 v = 0.0d0
      go to 520
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = 1.0d0
      sinp = 0.0d0
  660 h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      pqab = aqz-app(i)
      g = 1.d0/(ep(i)+ecd)
      p = (pqab*pqab+qperp2)*g
      if (p .le. auxvar) go to 140
      q0 = dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      q1 = 0.5d0*q0/p
      go to 160
  140 q = dp00p(i)/dsqrt(gp(i)+gcd)
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      q0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      q1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*q
  160 u = g*q1
      h0000 = h0000+q0
      h0001 = h0001+u
      h0003 = h0003-u*pqab
  180 continue
      h0001 = h0001*ecd*qperp
      h0003 = h0003*ecd
      p = dq*h0000
      g0001 = h0001*cosp+p*sing
      g0002 = h0001*sinp
      g0003 = h0003+p*cosg
      gout1 = gout1+dq00*h0000
      gout2 = gout2+dq01*g0001
      gout3 = gout3+dq01*g0002
      gout4 = gout4+dq01*g0003
 940  continue
      gout(1) = gout1
      t1 = gout2
      t2 = gout3
      t3 = gout4
      gout(2) = p11*t1+p21*t2+p31*t3
      gout(3) = p12*t1+p22*t2+p32*t3
      gout(4) = p13*t1+p23*t2+p33*t3
      return
      end
