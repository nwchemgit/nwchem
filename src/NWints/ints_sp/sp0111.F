      subroutine sp0111(gout)
c $Id$
c        *****  special fast routine for -p- loop for 0111 *****
      implicit none
      double precision gout(*)
c
      double precision g(64),h(64)
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
      integer k, l, i, j
      integer ind, i1, i2, i3, kq1, loopkl, loopjk
      double precision dzero, done, sixty, tenm12
c
      double precision qqq, qqtest, qperp, qperp2, q
      double precision dq00, dq01, dq10, dq11, dq01dd
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
      double precision t4, t5, t9, t13
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
  660 p1 = 0.d0
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
      gabcd = 1.d0/(eab+ecd)
      p = gabcd*(pqab2+qperp2)
      if (p .le. auxvar) go to 140
      f0 = conp(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = gabcd/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      f3 = 2.5d0*f2*gtx
      go to 160
  140 q = conp(i)/dsqrt(gp(i)+gcd)
      gy = gabcd*q
      ggy = gabcd*gy
      gggy = gabcd*ggy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
      f3 = (a3(n+1)+theta*b3(n+1)-theta3*c3(n+1)+theta4*c3(n+2))*gggy
  160 continue
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f2pqa2 = f2*pqab2
      p1 = p1+f0 *dp00
      p2 = p2+f1 *dp00
      p3 = p3+f2 *dp00
      p4 = p4+f1pqab*dp00
      p5 = p5+f2pqab*dp00
      p6 = p6+f2pqa2*dp00
      q1 = q1+f0 *bp
      q2 = q2+f1 *bp
      q3 = q3+f2 *bp
      q4 = q4+f1pqab*bp
      q5 = q5+f2pqab*bp
      q6 = q6+f2pqa2*bp
      r1 = r1+f1
      r2 = r2+f2
      r3 = r3+f3
      r4 = r4+f1pqab
      r5 = r5+f2pqab
      r6 = r6+f3pqab
      r7 = r7+f2pqa2
      r8 = r8+f3*pqab2
      r9 = r9+f3pqab*pqab2
  180 continue
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      qecd = qperp*ecd
      qecd2 = qperp*ecd2
      q2ecd = qperp2*ecd
      q2ecd2 = qperp2*ecd2
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
      h( 39) = 0.5d0*qecd2*r2
      h( 44) = -0.5d0*ecd2*r5
      h( 27) = h( 39)-qperp*h( 35)
      h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
      h( 24) = h( 44)+q2ecd2*r6
      h( 56) = h( 39)-qecd2*(r8+q5)
      h( 22) = h( 27)+h( 39)+h( 39)-q2ecd2*qperp*r3
      h( 32) = h( 27)-qecd2*r8
      h( 54) = h( 59)+q2ecd2*(r6+q3)
      h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
      g( 22) = v44*h( 22)+v47*h( 27)
      g( 23) = v54*h( 22)+v57*h( 27)
      g( 27) = v74*h( 22)+v77*h( 27)
      g( 24) = cosp*h( 24)
      g( 28) = sinp*h( 24)
      g( 38) = v45*h( 39)
      g( 39) = v55*h( 39)
      g( 43) = -g( 38)
      g( 40) = u12*h( 44)
      g( 44) = cosp*h( 44)
      g( 54) = v44*h( 54)+v47*h( 59)
      g( 55) = v54*h( 54)+v57*h( 59)
      g( 59) = v74*h( 54)+v77*h( 59)
      g( 56) = cosp*h( 56)
      g( 60) = sinp*h( 56)
      g( 64) = h( 64)
      g(  6) = v44*h(  6)+v47*h( 11)
      g(  7) = v54*h(  6)+v57*h( 11)
      g( 11) = v74*h(  6)+v77*h( 11)
      g(  8) = cosp*h(  8)
      g( 12) = sinp*h(  8)
      g( 16) = h( 16)
      g( 18) = cosp*h( 18)
      g( 19) = sinp*h( 18)
      g( 20) = h( 20)
      g( 34) = u12*h( 35)
      g( 35) = cosp*h( 35)
      g( 50) = cosp*h( 50)
      g( 51) = sinp*h( 50)
      g( 52) = h( 52)
      g(  2) = cosp*h(  2)
      g(  3) = sinp*h(  2)
      g(  4) = h(  4)
      g( 49) = h( 49)
      g(  1) = h(  1)
      h( 22) = g( 22)
      h( 23) = g( 23)
      h( 24) = g( 24)
      h( 27) = g( 27)
      h( 28) = g( 28)
      h( 38) = g( 38)
      h( 39) = g( 39)
      h( 40) = g( 40)
      h( 43) = g( 43)
      h( 44) = g( 44)
      h( 18) = g( 18)
      h( 19) = g( 19)
      h( 34) = g( 34)
      h( 35) = g( 35)
      g( 22) = cosp*h( 22)+u12*h( 38)
      g( 23) = cosp*h( 23)+u12*h( 39)
      g( 24) = cosp*h( 24)+u12*h( 40)
      g( 27) = cosp*h( 27)+u12*h( 43)
      g( 28) = cosp*h( 28)+u12*h( 44)
      g( 32) = cosp*h( 32)
      g( 38) = sinp*h( 22)+cosp*h( 38)
      g( 39) = sinp*h( 23)+cosp*h( 39)
      g( 40) = sinp*h( 24)+cosp*h( 40)
      g( 43) = sinp*h( 27)+cosp*h( 43)
      g( 44) = sinp*h( 28)+cosp*h( 44)
      g( 48) = sinp*h( 32)
      g( 18) = cosp*h( 18)+u12*h( 34)
      g( 19) = cosp*h( 19)+u12*h( 35)
      g( 20) = cosp*h( 20)
      g( 34) = sinp*h( 18)+cosp*h( 34)
      g( 35) = sinp*h( 19)+cosp*h( 35)
      g( 36) = sinp*h( 20)
      g( 17) = cosp*h( 17)
      g( 33) = sinp*h( 17)
      go to 2000
  920 g( 17) = h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = h( 20)
      g( 22) = h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = h( 27)
      g( 28) = dzero
      g( 32) = h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
      g(  1) = h(  1)
      g(  2) = h(  2)
      g(  3) = dzero
      g(  4) = h(  4)
      g(  6) = h(  6)
      g(  7) = dzero
      g(  8) = h(  8)
      g( 11) = h( 11)
      g( 12) = dzero
      g( 16) = h( 16)
      go to 2000
 1000 g( 17) = -h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = -h( 20)
      g( 22) = -h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = -h( 27)
      g( 28) = dzero
      g( 32) = -h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = -h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = -h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = -h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
      g(  1) = h(  1)
      g(  2) = -h(  2)
      g(  3) = dzero
      g(  4) = h(  4)
      g(  6) = h(  6)
      g(  7) = dzero
      g(  8) = -h(  8)
      g( 11) = h( 11)
      g( 12) = dzero
      g( 16) = h( 16)
2000  r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      do 2001 kq1=2,50,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+12) = g(kq1+ 6)
2001      g(kq1+13) = g(kq1+10)
      if (rcdsq) 720,720,700
700   do 701 kq1=1,49,16
          t1=g(kq1)
          t5=g(kq1+4)
          t13=g(kq1+12)
          t5=t5+r13*t1
          t13=t13+r33*t1
          t2=g(kq1+1)
          g(kq1+5)=g(kq1+5)+r13*t2+r14*t5
          g(kq1+13)=g(kq1+13)+r33*t2+r14*t13
          g(kq1+1)=t2+r14*t1
          t3=g(kq1+2)
          g(kq1+6)=g(kq1+6)+t3*r13
          g(kq1+14)=g(kq1+14)+t3*r33
          t9=g(kq1+8)
          g(kq1+9)=g(kq1+9)+r14*t9
          g(kq1+11)=g(kq1+11)+r34*t9
          t4=g(kq1+3)
          g(kq1+7)=g(kq1+7)+t4*r13+t5*r34
          g(kq1+4)=t5
          g(kq1+15)=g(kq1+15)+t4*r33+t13*r34
          g(kq1+12)=t13
701       g(kq1+3)=t4+t1*r34
720   do 721 kq1=2,62,4
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq11
          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq11
721       gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq11
      dq01dd=dq01-dq11
      do 722 kq1=1,49,16
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq00
          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq01dd
          gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq01dd
722       gout(kq1+3)=gout(kq1+3)+g(kq1+3)*dq01dd
      do 723 kq1=5,53,16
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq10
          gout(kq1+4)=gout(kq1+4)+g(kq1+4)*dq10
723       gout(kq1+8)=gout(kq1+8)+g(kq1+8)*dq10
  940 continue
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
      ind = -12
      do 2 j = 1,4
      ind = ind+12
      do 2 l = 1,4
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
   2  continue
      ind = -3
      do 3 loopjk = 1,16
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
   3  continue
      return
      end
