      subroutine sp1111(gout)
c $Id$
      implicit none
      double precision g(256),h(256)
      double precision dzero, done
      parameter (dzero=0.0d0)
      integer kq1off(10),kq2off(6),kq3off(6),kq4off(4),kq5off(6)
c
      double precision aa, ba, ca, abc1, ab, bb, cb, abc2
      double precision ac, bc, cc, abc3, ad, bd, cd, abc4
      double precision ae, be, ce, abc5, af, bf, cf, abc6
      common/tabint/
     +  aa(333),ba(333),ca(333),abc1,
     +  ab(333),bb(333),cb(333),abc2,
     +  ac(333),bc(333),cc(333),abc3,
     +  ad(333),bd(333),cd(333),abc4,
     +  ae(333),be(333),ce(333),abc5,
     +  af(333),bf(333),cf(333),abc6
c
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
      double precision qq, theta, cq, dq, bp, ap
      integer n
c
      double precision gout(*)
c
      double precision sixty, tenm12
      integer k, l, i
      integer ind, i1, i2, i3, kq1, loopkl
      integer kkq1, kkkq1, jq1, jkl, loopij, ijk
c
      double precision qqq, qqtest, qperp, qperp2, q
      double precision dq00, dq01, dq10, dq11
      double precision aqx, aqz, cosp, sinp
      double precision theta2, theta3, theta4
      double precision v77, v44, v45, v74, v47, v55, v57, v54, v
      double precision r13, r14, r33, r34, t1, t2, t3, u12
      double precision f0, f1, f2, f3, gy, ggy, gtx
      double precision ecd, pqab2, pqab, p
      double precision p1, p2, p3, p4, p5, p6
      double precision q1, q2, q3, q4, q5, q6
      double precision r1, r2, r3, r4, r5, r6, r7, r8, r9
      double precision eab, dp00, gabcd, gggy, hecd, ecd2, qecd
      double precision qecd2, q2ecd, q2ecd2
      double precision f1pqab, f2pqab, f3pqab, f2pqa2
      double precision t4, t5, t6, t7, t8, t9, t10, t11, t12, t13
      double precision v1, v2, v3, v4, v5, v6
      double precision w1, w2, w3, w4, w5, w6, w7, w8, w9
      double precision s1, s2, s3, s4, s6, s7, s8, s9
      double precision s10, s11, s12, s13, s14
      double precision c1, c2, c3, c4, c5, c6
      double precision a1, a2, a3, a4, a5, a6, a8, a9, a10
      double precision b1, b2, b3, b4, b5, b6, b8, b9, b10
      double precision dp01, dp10, apdp10, edp10
      double precision f1pqa2, f3pqa2, f2pqa3, f4pqa3, f4pqa2, f3pqa3
      double precision f4pqab, t14, f4, apbp, eab2, bpdp01, edp01
      double precision bqz, hecd2, hqecd, hqecd2, q3ecd, q3ecd2
      double precision temp, hxxyy, dq01x, dq10x, dq00x
c
      data done/1.0d0/
      data kq1off/1,17,49,65,81,113,161,193,209,241/
      data kq2off/33,97,129,145,177,225/
      data kq3off/1,49,81,161,193,241/
      data kq4off/17,65,113,209/
      data kq5off/33,97,129,145,177,225/
      data sixty,tenm12/60.0d0,1.0d-12/
c
c ***
c *** this is the fps version of sp1111.
c ***
c *** as much code as possible reduced to loops (>=4)
c *** to avoid ps cache misses and to enhance compiler
c *** optimisation. will probably run like a drain on
c *** the cray-1s]
c ***
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
660   p1 = 0.d0
      p2 = 0.d0
      p3 = 0.d0
      p4 = 0.d0
      p5 = 0.d0
      p6 = 0.d0
      q1 = 0.d0
      q2 = 0.d0
      q3 = 0.d0
      q4 = 0.d0
      q5 = 0.d0
      q6 = 0.d0
      r1 = 0.d0
      r2 = 0.d0
      r3 = 0.d0
      r4 = 0.d0
      r5 = 0.d0
      r6 = 0.d0
      r7 = 0.d0
      r8 = 0.d0
      r9 = 0.d0
      v1 = 0.d0
      v2 = 0.d0
      v3 = 0.d0
      v4 = 0.d0
      v5 = 0.d0
      v6 = 0.d0
      w1 = 0.d0
      w2 = 0.d0
      w3 = 0.d0
      w4 = 0.d0
      w5 = 0.d0
      w6 = 0.d0
      w7 = 0.d0
      w8 = 0.d0
      w9 = 0.d0
      s1 = 0.d0
      s2 = 0.d0
      s3 = 0.d0
      s4 = 0.d0
      s6 = 0.d0
      s7 = 0.d0
      s8 = 0.d0
      s9 = 0.d0
      s10 = 0.d0
      s11 = 0.d0
      s12 = 0.d0
      s13 = 0.d0
      s14 = 0.d0
      t1 = 0.d0
      t2 = 0.d0
      t3 = 0.d0
      t4 = 0.d0
      t5 = 0.d0
      t6 = 0.d0
      t7 = 0.d0
      t8 = 0.d0
      t9 = 0.d0
      t10 = 0.d0
      t11 = 0.d0
      t12 = 0.d0
      t13 = 0.d0
      t14 = 0.d0
      c1 = 0.d0
      c2 = 0.d0
      c3 = 0.d0
      c4 = 0.d0
      c5 = 0.d0
      c6 = 0.d0
      do 180 ind = 1,ngangb
      isml = ismlq+ismlp(ind)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
       eab = ep(ind)
      dp00 = dp00p(ind)
      dp01 = dp01p(ind)
      dp10 = dp10p(ind)
      ap = app(ind)
      bp = bpp(ind)
      pqab = aqz-ap
      pqab2 = pqab*pqab
      gabcd = 1.d0/(eab+ecd)
      p = gabcd*(qperp2+pqab2)
      if (p .le. auxvar) go to 140
      f0 = dsqrt(0.7853981625d0/(p*(gp(ind)+gcd)))*conp(ind)
      gtx = gabcd/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      f3 = 2.5d0*f2*gtx
      f4 = 3.5d0*f3*gtx
      go to 160
  140 q = conp(ind)/dsqrt(gp(ind)+gcd)
      gy = gabcd*q
      ggy = gabcd*gy
      gggy = gabcd*ggy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (aa(n+1)+theta*ba(n+1)-theta3*ca(n+1)+theta4*ca(n+2))*q
      f1 = (ab(n+1)+theta*bb(n+1)-theta3*cb(n+1)+theta4*cb(n+2))*gy
      f2 = (ac(n+1)+theta*bc(n+1)-theta3*cc(n+1)+theta4*cc(n+2))*ggy
      f3 = (ad(n+1)+theta*bd(n+1)-theta3*cd(n+1)+theta4*cd(n+2))*gggy
      f4 = (ae(n+1)+theta*be(n+1)-theta3*ce(n+1)+theta4*ce(n+2))*gggy*
     &     gabcd
  160 apbp = ap*bp
      eab2 = eab*eab
      bpdp01 = bp*dp01
      apdp10 = ap*dp10
      edp01 = eab*dp01
      edp10 = eab*dp10
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f4pqab = f4*pqab
      f1pqa2 = f1*pqab2
      f2pqa2 = f2*pqab2
      f3pqa2 = f3*pqab2
      f4pqa2 = f4*pqab2
      f2pqa3 = f2pqa2*pqab
      f3pqa3 = f3pqa2*pqab
      f4pqa3 = f4pqa2*pqab
      p1 = p1+f0 *dp00
      p2 = p2+f1 *dp00
      p3 = p3+f2 *dp00
      p4 = p4+f1pqab*dp00
      p5 = p5+f2pqab*dp00
      p6 = p6+f2pqa2*dp00
      r1 = r1+f1 *edp01
      r2 = r2+f2 *edp01
      r3 = r3+f3 *edp01
      r4 = r4+f1pqab *edp01
      r5 = r5+f2pqab *edp01
      r6 = r6+f3pqab *edp01
      r7 = r7+f2pqa2 *edp01
      r8 = r8+f3pqa2 *edp01
      r9 = r9+f3pqa3 *edp01
      w1 = w1+f1 *edp10
      w2 = w2+f2 *edp10
      w3 = w3+f3 *edp10
      w4 = w4+f1pqab *edp10
      w5 = w5+f2pqab *edp10
      w6 = w6+f3pqab *edp10
      w7 = w7+f2pqa2 *edp10
      w8 = w8+f3pqa2 *edp10
      w9 = w9+f3pqa3 *edp10
      s1 = s1+f0 *eab
      s2 = s2+f1 *eab
      s3 = s3+f2 *eab
      s4 = s4+f3 *eab
      s6 = s6+f1pqab*eab
      s7 = s7+f2pqab*eab
      s8 = s8+f3pqab*eab
      s9 = s9+f1pqa2*eab
      s10 = s10+f2pqa2*eab
      s11 = s11+f3pqa2*eab
      s12 = s12+f2pqa3*eab
      s13 = s13+f3pqa3*eab
      s14 = s14+f3pqa3*pqab*eab
      t1 = t1+f0 *eab2
      t2 = t2+f1 *eab2
      t3 = t3+f2 *eab2
      t4 = t4+f3 *eab2
      t5 = t5+f4 *eab2
      t6 = t6+f2pqab*eab2
      t7 = t7+f3pqab*eab2
      t8 = t8+f4pqab*eab2
      t9 = t9+f2pqa2*eab2
      t10 = t10+f3pqa2*eab2
      t11 = t11+f4pqa2*eab2
      t12 = t12+f3pqa3*eab2
      t13 = t13+f4pqa3*eab2
      t14 = t14+f4pqa3*pqab*eab2
      if (rabsq .eq. 0.0d0) go to 180
      q1 = q1+f0 *bpdp01
      q2 = q2+f1 *bpdp01
      q3 = q3+f2 *bpdp01
      q4 = q4+f1pqab*bpdp01
      q5 = q5+f2pqab*bpdp01
      q6 = q6+f2pqa2*bpdp01
      v1 = v1+f0 *apdp10
      v2 = v2+f1 *apdp10
      v3 = v3+f2 *apdp10
      v4 = v4+f1pqab*apdp10
      v5 = v5+f2pqab*apdp10
      v6 = v6+f2pqa2*apdp10
      c1 = c1+f0 *apbp
      c2 = c2+f1 *apbp
      c3 = c3+f2 *apbp
      c4 = c4+f1pqab*apbp
      c5 = c5+f2pqab*apbp
      c6 = c6+f2pqa2*apbp
  180 continue
      a1 = aqz*s2-s6
      a2 = aqz*s3-s7
      a3 = aqz*s4-s8
      a4 = aqz*s6-s9
      a5 = aqz*s7-s10
      a6 = aqz*s8-s11
      a8 = aqz*s10-s12
      a9 = aqz*s11-s13
      a10 = aqz*s13-s14
      bqz = aqz-rab
      b1 = bqz*s2-s6
      b2 = bqz*s3-s7
      b3 = bqz*s4-s8
      b4 = bqz*s6-s9
      b5 = bqz*s7-s10
      b6 = bqz*s8-s11
      b8 = bqz*s10-s12
      b9 = bqz*s11-s13
      b10 = bqz*s13-s14
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      hecd2 = 0.5d0*ecd2
      qecd = qperp*ecd
      hqecd = 0.5d0*qecd
      qecd2 = qperp*ecd2
      hqecd2 = 0.5d0*qecd2
      q2ecd = qperp2*ecd
      q3ecd = qperp*q2ecd
      q2ecd2 = qperp2*ecd2
      q3ecd2 = q2ecd2*qperp
      h(  1) = p1
      h(  2) = qecd*p2
      h(  4) = -ecd*p4
      h( 11) = hecd*(p1-ecd*p2)
      h(  6) = h( 11)+q2ecd2*p3
      h(  8) = -qecd2*p5
      h( 16) = h( 11)+ecd2*p6
      h( 17) = -qperp*r1
      h( 49) = r4+q1
      h( 35) = hecd*r1
      h( 18) = h( 35)-q2ecd*r2
      h( 20) = qecd*r5
      h( 50) = h( 20)+qecd*q2
      h( 52) = h( 35)-ecd*r7-ecd*q4
      h( 39) = hqecd2*r2
      h( 44) = -hecd2*r5
      h( 27) = h( 39)-qperp*h( 35)
      h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
      h( 24) = h( 44)+q2ecd2*r6
      h( 56) = h( 39)-qecd2*(r8+q5)
      h( 22) = h( 27)+h( 39)+h( 39)-q3ecd2*r3
      h( 32) = h( 27)-qecd2*r8
      h( 54) = h( 59)+q2ecd2*(r6+q3)
      h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
      h( 65) = -qperp*w1
      h(193) = w4+v1
      h(131) = hecd*w1
      h( 66) = h(131)-q2ecd*w2
      h( 68) = qecd*w5
      h(194) = h( 68)+qecd*v2
      h(196) = h(131)-ecd*w7-ecd*v4
      h(135) = hqecd2*w2
      h(140) = -hecd2*w5
      h( 75) = h(135)-qperp*h(131)
      h(203) = h(140)+hecd*(h(193)-ecd*v2)
      h( 72) = h(140)+q2ecd2*w6
      h(200) = h(135)-qecd2*(w8+v5)
      h( 70) = h( 75)+h(135)+h(135)-q3ecd2*w3
      h( 80) = h( 75)-qecd2*w8
      h(198) = h(203)+q2ecd2*(w6+v3)
      h(208) = h(203)+h(140)+h(140)+ecd2*(w9+v6)
      h(161) = 0.5d0*(s1-t2)
      h( 81) = h(161)+qperp2*t3
      h(113) = -qperp*(t6+b1)
      h(209) = -qperp*(t6+a1)
      h(241) = h(161)+t9+a4+b4+c1
      h(162) = hqecd*(s2-t3)
      h( 82) = h(162)-qecd*t3+q3ecd*t4
      temp = hecd*t6-q2ecd*t7
      h(114) = temp+hecd*b1-q2ecd*b2
      h(210) = temp+hecd*a1-q2ecd*a2
      h(242) = h(162)+qecd*(t10+a5+b5+c2)
      h( 99) = -hqecd*t3
      h(147) = h( 99)
      h(179) = hecd*(t6+b1)
      h(227) = hecd*(t6+a1)
      h(164) = hecd*(t6-s6)
      h( 84) = h(164)-q2ecd*t7
      temp = -hqecd*t3+qecd*t10
      h(116) = temp+qecd*b5
      h(212) = temp+qecd*a5
      h(244) = h(164)+ecd*(t6-t12-a8-b8-c4)+hecd*(a1+b1)
      h(103) = 0.25d0*ecd2*t3-0.5d0*q2ecd2*t4
      h(151) = h(103)
      h(183) = hqecd2*(t7+b2)
      h(231) = hqecd2*(t7+a2)
      h(108) = hqecd2*t7
      h(156) = h(108)
      h(188) = hecd2*(0.5d0*t3-t10-b5)
      h(236) = hecd2*(0.5d0*t3-t10-a5)
      hxxyy = 0.25d0*(ecd*(s1-t2)-ecd2*(s2-t3))
      h(171) = hxxyy+hecd2*t3
      h( 91) = hxxyy+0.5d0*(q2ecd*t3-q2ecd2*t4)
      temp = hqecd*(ecd*t7-t6)
      h(123) = temp+hqecd*(ecd*b2-b1)
      h(219) = temp+hqecd*(ecd*a2-a1)
      h(251) = hxxyy+hecd*(t9+a4+b4+c1)-hecd2*(t10+a5+b5+c2)
      h(166) = hxxyy+0.5d0*q2ecd2*(s3-t4)
      h( 86) = hxxyy+(hecd2+0.5d0*q2ecd)*t3+q2ecd2*(-3.d0*t4+
     +    0.5d0*s3+qperp2*t5)
      h(118) = 1.5d0*qecd2*(t7+b2)-hqecd*(t6+b1)-q3ecd2*(b3+t8)
      h(214) = 1.5d0*qecd2*(t7+a2)-hqecd*(t6+a1)-q3ecd2*(a3+t8)
      h(246) = hxxyy-hecd2*(qperp2*t4+t10+a5+b5)+hecd*(t9+a4+b4+c1-ecd*c
     +     2)+q2ecd2*(t11+0.5d0*s3+a6+b6+c3)
      h(168) = hqecd2*(t7-s7)
      h( 88) = 1.5d0*qecd2*t7-hqecd2*s7-q3ecd2*t8
      temp = hecd2*(0.5d0*t3-t10)+q2ecd2*(t11-0.5d0*t4)
      h(120) = temp-hecd2*b5+q2ecd2*b6
      h(216) = temp-hecd2*a5+q2ecd2*a6
      h(248) = qecd2*(1.5d0*t7-t13-a9-b9-c5)-hqecd2*(s7-a2-b2)
      h(176) = hxxyy+hecd2*(s10-t10)
      h( 96) = hxxyy-hecd2*(qperp2*t4+t10-s10)+0.5d0*q2ecd*t3+q2ecd2*
     +     t11
      h(128) = qecd2*(1.5d0*t7-t13-b9)-hqecd*(t6+b1)+hqecd2*b2
      h(224) = qecd2*(1.5d0*t7-t13-a9)-hqecd*(t6+a1)+hqecd2*a2
      h(256) = hxxyy+hecd2*(-3.d0*(a5+b5)+t3+s10-c2)+ecd2*(-3.d0*t10+
     +     t14+a10+b10+c6)+hecd*(t9+a4+b4+c1)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
 120  u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
      do 103 kq1=22,214,48
          g(kq1  ) = v44*h(kq1) + v47*h(kq1+5)
          g(kq1+1) = v54*h(kq1) + v57*h(kq1+5)
103       g(kq1+5) = v74*h(kq1) + v77*h(kq1+5)
      do 101 kq1=24,216,48
          g(kq1  ) = cosp*h(kq1)
101       g(kq1+4) = sinp*h(kq1)
      do 102 kq1=18,210,48
          g(kq1  ) = cosp*h(kq1)
102       g(kq1+1) = sinp*h(kq1)
      g( 80) = h( 80)
      g( 86) = v44*h( 86)+v47*h( 91)
      g( 87) = v54*h( 86)+v57*h( 91)
      g( 91) = v74*h( 86)+v77*h( 91)
      g( 88) = cosp*h( 88)
      g( 92) = sinp*h( 88)
      g( 96) = h( 96)
      g(102) = v45*h(103)
      g(103) = v55*h(103)
      g(107) = -g(102)
      g(104) = u12*h(108)
      g(108) = cosp*h(108)
      g(112) = dzero
      g(128) = h(128)
      g(134) = v45*h(135)
      g(135) = v55*h(135)
      g(139) = -g(134)
      g(136) = u12*h(140)
      g(140) = cosp*h(140)
      g(144) = dzero
      g(150) = v45*h(151)
      g(151) = v55*h(151)
      g(155) = -g(150)
      g(152) = u12*h(156)
      g(156) = cosp*h(156)
      g(160) = dzero
      g(176) = h(176)
      g(182) = v45*h(183)
      g(183) = v55*h(183)
      g(187) = -g(182)
      g(184) = u12*h(188)
      g(188) = cosp*h(188)
      g(192) = dzero
      g(198) = v44*h(198)+v47*h(203)
      g(199) = v54*h(198)+v57*h(203)
      g(203) = v74*h(198)+v77*h(203)
      g(200) = cosp*h(200)
      g(204) = sinp*h(200)
      g(230) = v45*h(231)
      g(231) = v55*h(231)
      g(235) = -g(230)
      g(232) = u12*h(236)
      g(236) = cosp*h(236)
      g(240) = dzero
      g(246) = v44*h(246)+v47*h(251)
      g(247) = v54*h(246)+v57*h(251)
      g(251) = v74*h(246)+v77*h(251)
      g(248) = cosp*h(248)
      g(252) = sinp*h(248)
      g( 38) = v45*h( 39)
      g( 39) = v55*h( 39)
      g( 43) = -g( 38)
      g( 40) = u12*h( 44)
      g( 44) = cosp*h( 44)
      g( 48) = dzero
      g( 54) = v44*h( 54)+v47*h( 59)
      g( 55) = v54*h( 54)+v57*h( 59)
      g( 59) = v74*h( 54)+v77*h( 59)
      g( 56) = cosp*h( 56)
      g( 60) = sinp*h( 56)
      g(  6) = v44*h(  6)+v47*h( 11)
      g(  7) = v54*h(  6)+v57*h( 11)
      g( 11) = v74*h(  6)+v77*h( 11)
      g(  8) = cosp*h(  8)
      g( 12) = sinp*h(  8)
      g( 68) = h( 68)
      g( 82) = cosp*h( 82)
      g( 83) = sinp*h( 82)
      g( 84) = h( 84)
      g( 98) = u12*h( 99)
      g( 99) = cosp*h( 99)
      g(100) = dzero
      g(116) = h(116)
      g(130) = u12*h(131)
      g(131) = cosp*h(131)
      g(132) = dzero
      g(146) = u12*h(147)
      g(147) = cosp*h(147)
      g(148) = dzero
      g(164) = h(164)
      g(178) = u12*h(179)
      g(179) = cosp*h(179)
      g(180) = dzero
      g(194) = cosp*h(194)
      g(195) = sinp*h(194)
      g(226) = u12*h(227)
      g(227) = cosp*h(227)
      g(228) = dzero
      g(242) = cosp*h(242)
      g(243) = sinp*h(242)
      g( 34) = u12*h( 35)
      g( 35) = cosp*h( 35)
      g( 36) = dzero
      g( 50) = cosp*h( 50)
      g( 51) = sinp*h( 50)
      g(  2) = cosp*h(  2)
      g(  3) = sinp*h(  2)
      g( 65) = h( 65)
      g( 81) = h( 81)
      g( 97) = dzero
      g(113) = h(113)
      g(129) = dzero
      g(145) = dzero
      g(161) = h(161)
      g(177) = dzero
      g(225) = dzero
      g( 33) = dzero
      h( 80) = cosp*g( 80)
      h( 96) = cosp*g( 96)
      h(112) =           u12*g(176)
      h(128) = cosp*g(128)
      h(144) = sinp*g( 80)
      h(160) = sinp*g( 96)
      h(176) =           cosp*g(176)
      h(192) = sinp*g(128)
      do 121 kq1=70,118,16
          h(kq1   ) = cosp*g(kq1   ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1   ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+ 1) + u12*g(kq1+65)
          h(kq1+65) = sinp*g(kq1+ 1) + cosp*g(kq1+65)
          h(kq1+ 2) = cosp*g(kq1+ 2) + u12*g(kq1+66)
          h(kq1+66) = sinp*g(kq1+ 2) + cosp*g(kq1+66)
          h(kq1+ 5) = cosp*g(kq1+ 5) + u12*g(kq1+69)
          h(kq1+69) = sinp*g(kq1+ 5) + cosp*g(kq1+69)
          h(kq1+ 6) = cosp*g(kq1+ 6) + u12*g(kq1+70)
          h(kq1+70) = sinp*g(kq1+ 6) + cosp*g(kq1+70)
121   continue
      h( 68) = cosp*g( 68)
      h( 84) = cosp*g( 84)
      h(100) =           u12*g(164)
      h(116) = cosp*g(116)
      h(132) = sinp*g( 68)
      h(148) = sinp*g( 84)
      h(164) =           cosp*g(164)
      h(180) = sinp*g(116)
      do 122 kq1=66,114,16
          h(kq1   ) = cosp*g(kq1  ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1  ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+1) + u12*g(kq1+65)
122       h(kq1+65) = sinp*g(kq1+1) + cosp*g(kq1+65)
      do 1221 kq1=2,50,16
          h(kq1   ) = g(kq1   )
          h(kq1+ 1) = g(kq1+ 1)
          h(kq1+ 4) = g(kq1+ 4)
          h(kq1+ 5) = g(kq1+ 5)
          h(kq1+ 6) = g(kq1+ 6)
1221      h(kq1+ 9) = g(kq1+ 9)
      do 1222 kq1=12,60,16
          h(kq1    ) = g(kq1    )
          h(kq1+182) = g(kq1+182)
          h(kq1+183) = g(kq1+183)
          h(kq1+186) = g(kq1+186)
          h(kq1+187) = g(kq1+187)
1222      h(kq1+188) = g(kq1+188)
      do 1223 kq1=203,251,16
          h(kq1  ) = g(kq1  )
1223      h(kq1+1) = g(kq1+1)
      h( 65) = cosp*g( 65)
      h( 81) = cosp*g( 81)
      h( 97) =           u12*g(161)
      h(113) = cosp*g(113)
      h(129) = sinp*g( 65)
      h(145) = sinp*g( 81)
      h(161) =           cosp*g(161)
      h(177) = sinp*g(113)
      h( 48) = g( 48)
      h( 36) = g( 36)
      h(228) = g(228)
      h(240) = g(240)
      h(225) = g(225)
      h( 33) = g( 33)
      do 123 kq1=22,214,64
          g(kq1   ) = cosp*h(kq1   ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1   ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+ 1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+ 1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+ 2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+ 2) + cosp*h(kq1+18)
          g(kq1+ 5) = cosp*h(kq1+ 5) + u12* h(kq1+21)
          g(kq1+21) = sinp*h(kq1+ 5) + cosp*h(kq1+21)
          g(kq1+ 6) = cosp*h(kq1+ 6) + u12* h(kq1+22)
          g(kq1+22) = sinp*h(kq1+ 6) + cosp*h(kq1+22)
          g(kq1+10) = cosp*h(kq1+10) + u12* h(kq1+26)
123       g(kq1+26) = sinp*h(kq1+10) + cosp* h(kq1+26)
      do 124 kq1=17,209,64
          g(kq1   ) = cosp*h(kq1  ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1  ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+2) + cosp*h(kq1+18)
          g(kq1+ 3) = cosp*h(kq1+3) + u12* h(kq1+19)
124       g(kq1+19) = sinp*h(kq1+3) + cosp* h(kq1+19)
      do 125 kq1=49,177,64
          kkq1=kq1
          do 126 kkkq1=1,32
              g(kkq1)=h(kkq1)
              kkq1=kkq1+1
126       continue
125   continue
      do 127 kq1=1,16
          g(kq1)=h(kq1)
127       g(kq1+240)=h(kq1+240)
      goto 2000
920   do 921 kkq1=1,10
          kq1=kq1off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) = h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) = h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
921       g(kq1+15) = h(kq1+15)
      do 922 kkq1=1,6
          kq1=kq2off(kkq1)
          g(kq1   ) = dzero
          g(kq1+ 1) = dzero
          g(kq1+ 2) = h(kq1+ 2)
          g(kq1+ 3) = dzero
          g(kq1+ 5) = dzero
          g(kq1+ 6) = h(kq1+ 6)
          g(kq1+ 7) = dzero
          g(kq1+10) = dzero
          g(kq1+11) = h(kq1+11)
922       g(kq1+15) = dzero
      go to 2000
1000  do 1001 kkq1=1,6
          kq1=kq3off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) =-h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) =-h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
1001      g(kq1+15) = h(kq1+15)
      do 1002 kkq1=1,4
          kq1=kq4off(kkq1)
          g(kq1   ) = -h(kq1   )
          g(kq1+ 1) =  h(kq1+ 1)
          g(kq1+ 2) =  dzero
          g(kq1+ 3) = -h(kq1+ 3)
          g(kq1+ 5) = -h(kq1+ 5)
          g(kq1+ 6) =  dzero
          g(kq1+ 7) =  h(kq1+ 7)
          g(kq1+10) = -h(kq1+10)
          g(kq1+11) =  dzero
1002      g(kq1+15) = -h(kq1+15)
      do 1003 kkq1=1,6
          kq1=kq5off(kkq1)
          g(kq1   ) =  dzero
          g(kq1+ 1) =  dzero
          g(kq1+ 2) =  h(kq1+ 2)
          g(kq1+ 3) =  dzero
          g(kq1+ 5) =  dzero
          g(kq1+ 6) = -h(kq1+ 6)
          g(kq1+ 7) =  dzero
          g(kq1+10) =  dzero
          g(kq1+11) =  h(kq1+11)
1003      g(kq1+15) =  dzero
          g(99)=-g(99)
          g(108)=-g(108)
          g(147)=-g(147)
          g(156)=-g(156)
          g(103)=-g(103)
          g(151)=-g(151)
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      do 2001 kq1=2,242,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+12) = g(kq1+ 6)
2001      g(kq1+13) = g(kq1+10)
      if (rcdsq) 1200,1200,1300
1300  do 1301 kq1=1,4
          kkq1=kq1
          do 1302 jq1=1,16
              g(kkq1+4) = r13*g(kkq1) + g(kkq1+4)
              g(kkq1+12)= r33*g(kkq1) + g(kkq1+12)
1302          kkq1=kkq1+16
1301  continue
c ***
      do 1303 kq1=1,253,4
          g(kq1+1) = r14*g(kq1) + g(kq1+1)
1303      g(kq1+3) = r34*g(kq1) + g(kq1+3)
c ***
1200  do 1201 kq1=1,256
1201      gout(kq1) = dq11*g(kq1) + gout(kq1)
c1200  continue
c      call daxpy(256,dq11,g,1,gout,1)
      dq01x=dq01-dq11
      do 1202 kq1=2,242,16
          gout(kq1  ) = dq01x*g(kq1  ) + gout(kq1  )
          gout(kq1+1) = dq01x*g(kq1+1) + gout(kq1+1)
1202      gout(kq1+2) = dq01x*g(kq1+2) + gout(kq1+2)
      dq10x=dq10-dq11
      dq00x=dq00-dq11
      do 1203 kq1=1,241,16
          gout(kq1   ) = dq00x*g(kq1   ) + gout(kq1   )
          gout(kq1+ 4) = dq10x*g(kq1+ 4) + gout(kq1+ 4)
          gout(kq1+ 8) = dq10x*g(kq1+ 8) + gout(kq1+ 8)
1203      gout(kq1+12) = dq10x*g(kq1+12) + gout(kq1+12)
940   continue
c
c     --------------------------
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 1111
c
c
      i1 = 64
      i2 = 128
      i3 = 192
      do 1 jkl = 1,64
      i1 = i1+1
      i2 = i2+1
      i3 = i3+1
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
    1 continue
      ind = -48
      do 2 i = 1,4
      ind = ind+48
      do 2 loopkl = 1,16
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
    2 continue
      ind = -12
      do 3 loopij = 1,16
      ind = ind+12
      do 3 l = 1,4
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
    3 continue
      ind = -3
      do 4 ijk = 1,64
      ind = ind+4
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
    4 continue
c
      return
      end
