      subroutine sp0011(gout)
c $Id$
c        *****  special fast routine for -p- loop for 0011 *****
      implicit none
      double precision gout(*)
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
      double precision a0, b0, c0, a1, b1, c1, a2, b2, c2, a3, b3, c3
      double precision a4, b4, c4, a5, b5, c5, abc1, abc2, abc3
      double precision abc4, abc5, abc6
c
      common/tabint/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
c
      integer k, l, i
      integer ind, i1, i2, i3
      double precision dzero, done, sixty, tenm12
      double precision qqq, qqtest, qperp, qperp2, q
      double precision dq00, dq01, dq10, dq11
      double precision aqx, aqz, cosp, sinp
      double precision h0000, h0001, h0003, h0011, h0013
      double precision h0022, h0033
      double precision g0000, g0001, g0002, g0003, g0010, g0011, g0012
      double precision g0013, g0020, g0021, g0022, g0023
      double precision g0030, g0031, g0032, g0033
      double precision theta2, theta3, theta4
      double precision v77, v44, v74, v47, v57, v54, v
      double precision r13, r14, r33, r34, t1, t2, t3, u12
      double precision f0, f1, f2, g, gy, ggy, gtx
      double precision ecd, pqab2, pqab, p
c
      data dzero,done/0.0d0,1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
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
  640 cosp = done
      sinp = 0.0d0
 660  h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0011 = 0.d0
      h0013 = 0.d0
      h0033 = 0.d0
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      g = 1.d0/(ep(i)+ecd)
      p = g*(pqab2+qperp2)
      g = g*ecd
      if (p .le. auxvar) go to 140
      f0 = dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = g/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      go to 160
  140 q = dp00p(i)/dsqrt(gp(i)+gcd)
      gy = g*q
      ggy = g*gy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
  160 h0000 = h0000+f0
      h0001 = h0001+f1
      h0003 = h0003-f1*pqab
      h0011 = h0011+f2
      h0013 = h0013-f2*pqab
      h0033 = h0033+f2*pqab2
  180 continue
      h0022 = 0.5d0*ecd*(h0000-h0001)
      h0001 = h0001*qperp
      h0011 = h0011*qperp2+h0022
      h0013 = h0013*qperp
      h0033 = h0033+h0022
      if(sinp)120,100,120
 100  if(cosp)1000,120,920
 120  u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      g0011 = v44*h0011+v47*h0022
      g0012 = v54*h0011+v57*h0022
      g0022 = v74*h0011+v77*h0022
      g0013 = cosp*h0013
      g0023 = sinp*h0013
      g0033 = h0033
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0000 = h0000
      go to 2000
  920 g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
      go to 2000
1000  g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = -h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      g0010 = g0001
      g0020 = g0002
      g0021 = g0012
      g0030 = g0003
      g0031 = g0013
      g0032 = g0023
      if (rcdsq) 220,220,200
 200  g0010 = g0010+r13*g0000
      g0011 = g0011+r13*g0001
      g0012 = g0012+r13*g0002
      g0013 = g0013+r13*g0003
      g0030 = g0030+r33*g0000
      g0031 = g0031+r33*g0001
      g0032 = g0032+r33*g0002
      g0033 = g0033+r33*g0003
      g0001 = g0001+r14*g0000
      g0011 = g0011+r14*g0010
      g0021 = g0021+r14*g0020
      g0031 = g0031+r14*g0030
      g0003 = g0003+r34*g0000
      g0013 = g0013+r34*g0010
      g0023 = g0023+r34*g0020
      g0033 = g0033+r34*g0030
220   gout( 1) = gout( 1)+g0000*dq00
      gout( 2) = gout( 2)+g0001*dq01
      gout( 3) = gout( 3)+g0002*dq01
      gout( 4) = gout( 4)+g0003*dq01
      gout( 5) = gout( 5)+g0010*dq10
      gout( 6) = gout( 6)+g0011*dq11
      gout( 7) = gout( 7)+g0012*dq11
      gout( 8) = gout( 8)+g0013*dq11
      gout( 9) = gout( 9)+g0020*dq10
      gout( 10) = gout( 10)+g0021*dq11
      gout( 11) = gout( 11)+g0022*dq11
      gout( 12) = gout( 12)+g0023*dq11
      gout( 13) = gout( 13)+g0030*dq10
      gout( 14) = gout( 14)+g0031*dq11
      gout( 15) = gout( 15)+g0032*dq11
      gout( 16) = gout( 16)+g0033*dq11
 940  continue
      ind = 0
      do 700 l = 1,4
      ind = ind+1
      i1 = 4+ind
      i2 = 8+ind
      i3 = 12+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
  700 continue
      ind = -3
      do 720 k = 1,4
      ind = ind+4
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i3 ) = p13*t1+p23*t2+p33*t3
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
  720 continue
      return
      end
