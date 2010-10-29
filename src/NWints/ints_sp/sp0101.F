      subroutine sp0101(gout)
c $Id$
c        *****  special fast routine for -p- loop for 0101 *****
      implicit none
      double precision gout(*)
c
c
      double precision a0, b0, c0, a1, b1, c1, a2, b2, c2, a3, b3, c3
      double precision a4, b4, c4, a5, b5, c5, abc1, abc2, abc3
      double precision abc4, abc5, abc6
      common/tabint/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
c
#include "const.fh"
#include "auxvar.fh"
#include "miscg.fh"
#include "ginf.fh"
#include "pgeom.fh"
#include "shllfo.fh"
#include "spgeom.fh"
#include "qgeom.fh"
#include "maxc.fh"
c
      double precision qq, theta, cq, dq, bp
      integer n
c
      integer k, l, i, loopkl, j
      integer ind, i1, i2, i3
      double precision dzero, done, sixty, tenm12
      double precision qqq, qqtest, qperp, qperp2, q
      double precision dq00, dq01, dq10, dq11
      double precision aqx, aqz, cosp, sinp
      double precision h0000, h0001, h0003, h0100, h0101, h0103, h0300
      double precision h0301, h0303
      double precision g0000, g0001, g0002, g0003
      double precision h0102, h0202, h0201, eab, dp00
      double precision g03, g0201, g0202, g0200, g0203
      double precision g0301, g0302, g0300, g0303
      double precision g0101, g0102, g0100, g0103
      double precision theta2, theta3, theta4
      double precision v
      double precision r14, r34, t1, t2, t3, u12
      double precision f0, f1, f2, g, gy, ggy, gtx
      double precision ecd, pqab2, pqab, p
c
      data dzero/0.0d0/,done/1.0d0/
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
  660 h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0100 = 0.d0
      h0101 = 0.d0
      h0103 = 0.d0
      h0300 = 0.d0
      h0301 = 0.d0
      h0303 = 0.d0
c        *****  begin -p- loop                   *****
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      eab = ep(i)
      dp00 = dp00p(i)
      bp = bpp(i)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      g = 1.d0/(ep(i)+ecd)
      p = (pqab2+qperp2)*g
      if (p .le. auxvar) go to 140
      f0 = conp(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = g/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      go to 160
  140 q = conp(i)/dsqrt(gp(i)+gcd)
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
  160 continue
      g03 = -pqab*f1
      h0000 = h0000+f0 *dp00
      h0001 = h0001+f1 *dp00
      h0003 = h0003+g03*dp00
      h0100 = h0100-f1
      h0101 = h0101-f2
      h0103 = h0103+pqab*f2
      h0300 = h0300-g03+bp*f0
      h0301 = h0301+bp*f1
      h0303 = h0303-pqab2*f2+bp*g03
  180 continue
      p = qperp*ecd
      h0001 = h0001*p
      h0003 = h0003*ecd
      h0202 = -0.5d0*ecd*h0100
      h0100 = h0100*qperp
      h0101 = h0101*qperp2*ecd
      h0103 = h0103*p
      h0301 = h0301*p
      h0303 = h0303*ecd
      h0301 = h0301+h0103
      h0101 = h0101+h0202
      h0303 = h0303+h0202
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      g0101 = cosp*h0101
      g0102 = sinp*h0101
      g0201 = u12*h0202
      g0202 = cosp*h0202
      g0301 = cosp*h0301
      g0302 = sinp*h0301
      g0303 = h0303
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0300 = h0300
      g0000 = h0000
      h0101 = g0101
      h0102 = g0102
      h0201 = g0201
      h0202 = g0202
      g0101 = cosp*h0101+u12*h0201
      g0102 = cosp*h0102+u12*h0202
      g0103 = cosp*h0103
      g0201 = sinp*h0101+cosp*h0201
      g0202 = sinp*h0102+cosp*h0202
      g0203 = sinp*h0103
      g0100 = cosp*h0100
      g0200 = sinp*h0100
      go to 2000
  920 g0100 = h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = h0301
      g0302 = dzero
      g0303 = h0303
      g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      go to 2000
 1000 g0100 = -h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = -h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = -h0301
      g0302 = dzero
      g0303 = h0303
      g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
2000  r14 = dq*sing
      r34 = dq*cosg
      if (rcdsq) 720,720,700
  700 g0001 = g0001+r14*g0000
      g0101 = g0101+r14*g0100
      g0201 = g0201+r14*g0200
      g0301 = g0301+r14*g0300
      g0003 = g0003+r34*g0000
      g0103 = g0103+r34*g0100
      g0203 = g0203+r34*g0200
      g0303 = g0303+r34*g0300
  720 gout( 1) = gout( 1)+g0000*dq00
      gout( 2) = gout( 2)+g0001*dq01
      gout( 3) = gout( 3)+g0002*dq01
      gout( 4) = gout( 4)+g0003*dq01
      gout( 17) = gout( 17)+g0100*dq00
      gout( 18) = gout( 18)+g0101*dq01
      gout( 19) = gout( 19)+g0102*dq01
      gout( 20) = gout( 20)+g0103*dq01
      gout( 33) = gout( 33)+g0200*dq00
      gout( 34) = gout( 34)+g0201*dq01
      gout( 35) = gout( 35)+g0202*dq01
      gout( 36) = gout( 36)+g0203*dq01
      gout( 49) = gout( 49)+g0300*dq00
      gout( 50) = gout( 50)+g0301*dq01
      gout( 51) = gout( 51)+g0302*dq01
      gout( 52) = gout( 52)+g0303*dq01
 940  continue
c ***
c ***
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 0101
c
c
c
      ind = 0
      do 1 loopkl = 1,16
      ind = ind+1
      i1 = 16+ind
      i2 = 32+ind
      i3 = 48+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
   1  continue
      ind = -15
      do 2 j = 1,4
      ind = ind+16
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
   2  continue
      return
      end
