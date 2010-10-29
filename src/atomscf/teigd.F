      subroutine teigd(pcap,qcap,u,t,dt,dos)
C$Id$
      implicit none
      logical klnemn
      real *8  pcap(*),qcap(*),u(*),t(*),dos(*),dt(*)
c.......................................................................
c
c     two-electron integral routine for s,p,d, and f functions.
c.......................................................................
#include "cguess.fh"
c
      integer i, k, l, im, m, n
      integer j, nstep1, nstep2, kin, kmx, kl, mmx, mmxp, mn
      integer nmx, ntest, min
      logical open
c
      real *8 twopow(14), pifac(14)
      real *8 pi, on6, on15, on70, on35, pot, potn, factkl
      real *8 prfac1, prfac2, prfac3, prfac4
      real *8 zp, zq, zpq, zr, zpr, zqr, zs, zps, zrs, zqs, zpqrs
      real *8 zpqrs2, zprzqs, zpszqr
      real *8 xterm, factmn, term
      real *8 xfac1, xfac2, xfac3, xfac11, xfac21, xfsum
      real *8 x0, y0, pj, qj, y1, y2, y3, x2, x4, y01, y02
      real *8 y21, y22, y4, y5, y6
c
c.......................................................................
c     angular factors for exchange integrals. obtained as sum of squares
c     of slater coefficients c(kappa;l1,m1;l2,m2) divided by
c     2*(2*l1+1)*(2*l2+1) for any kappa, l1, and l2.
c.......................................................................
      real *8 ss0, sp1, pp0, pp2, sd2, pd1, pd3, dd0, dd2
      real *8 sf3, pf2, pf4, df1, df3, df5, ff0, ff2, ff4, ff6
      real *8 a, b
      real *8 f0pol, f2pol, df1pol, df3pol
      data ss0,sp1,pp0,pp2,sd2,pd1,pd3,dd0,dd2,
     x     sf3,pf2,pf4,df1,df3,df5,ff0,ff2,ff4,ff6
     x/ .50000000000d+00, .16666666667d+00, .16666666667d+00,
     x  .66666666667d-01, .10000000000d+00, .66666666667d-01,
     x  .42857142857d-01, .10000000000d+00, .28571428571d-01,
     x  .71428571429d-01, .42857142857d-01,
     x  .31746031746d-01, .42857142857d-01, .19047619048d-01,
     x  .21645021645d-01, .71428571429d-01, .19047619048d-01,
     x  .12987012987d-01, .16650016650d-01/
      f0pol(a,b) = 3*(16*a**6+104*a**5*b+286*a**4*b**2+429*(a*b)
     +             **3+286*a**2*b**4+104*a*b**5+16*b**6)
      f2pol(a,b) = 8*a**4 + 52*a**3*b + 143*(a*b)**2 + 52*a*b**3 +
     +             8*b**4
      df1pol(a,b) = 8*a**4 + 44*a**3*b + 99*(a*b)**2 + 44*a*b**3 +
     +              8*b**4
      df3pol(a,b) = 2*a**2 + 11*a*b + 2*b**2
c.......................................................................
c     two-electron integral routine for lcgo atom scf.
c     restricted to principal quantum numbers 1,2 and 3 for respectively
c     s,p, and d orbitals.
c.......................................................................
      pi = 3.14159265d0
      on6 = 1.0d0/6.0d0
      on15 = 1.0d0/15.d0
      on70 = 1.0d0/70.d0
      on35 = 2*on70
      pifac(1) = sqrt(3.14159265d0)
      twopow(1) = 1
      do i = 2 , 14
         twopow(i) = twopow(i-1)*0.5d0
         pifac(i) = pifac(i-1)*(2*i-1)
      enddo
c.......................................................................
c
c     this part sets up the coefficients lambda,p,q and mu,r,s.
c.......................................................................
      j = 0
      nstep1 = 0
      kmx = 0
      factkl = 1
      kl = 0
      pot = 0.0d0
      potn = 0.0d0
      cin = 0.0d0
      do i = 1 , nsym
         prfac1 = twopow(i+1)*factkl
         kin = kmx + 1
         kmx = kin + nbas(i) - 1
         do k = kin , kmx
            zp = zeta(k)
            do l = kin , k
               kl = kl + 1
               pcap(kl) = 0.0d0
               qcap(kl) = 0.0d0
               zq = zeta(l)
               zpq = zp + zq
               prfac2 = prfac1*u(kl)
               xfac1 = prfac2*zpq**i*pifac(i)
               nstep2 = 0
               mmx = 0
               factmn = 1.0d0
               mn = 0
               do im = 1 , i
                  open = (nosh(i).ne.0 .and. nosh(im).ne.0)
                  prfac3 = prfac2*pifac(im)*twopow(im)*factmn
                  xfac2 = xfac1*factmn*twopow(im+1)
                  min = mmx + 1
                  mmx = min - 1 + nbas(im)
                  mmxp = mmx
                  if (im.eq.i) mmxp = k
                  do m = min , mmxp
                     zr = zeta(m)
                     zpr = zp + zr
                     zqr = zq + zr
                     nmx = l
                     if (m.lt.k) nmx = m
                     do n = min , nmx
                        mn = mn + 1
                        klnemn = (kl.ne.mn)
                        j = j + 1
c.......................................................................
c
c     j is the number label of the matrix elements to be calculated
c     i=lambda+1,k=p,l=q,im=mu+1,m=r,n=s
c.......................................................................
                        zs = zeta(n)
                        zrs = zr + zs
                        zqs = zq + zs
                        zps = zp + zs
                        zpqrs = zpq + zrs
                        zpqrs2 = 2*zpqrs**2
                        zprzqs = zpr*zqs
                        zpszqr = zps*zqr
                        xterm = (1.0d0/sqrt(zpqrs))**(2*(i+im)-3)
                        prfac4 = prfac3*u(mn)*xterm
                        xfac3 = xfac2*u(mn)*xterm
                        xfac11 = xfac3*(zrs/zprzqs)**im
                        xfac21 = xfac3*(zrs/zpszqr)**im
                        xfsum = xfac11 + xfac21
                        ntest = i*(i-1)/2 + im
                        go to (30,40,50,60,70,80,90,100,110,120,140) ,
     +                         ntest
c.......................................................................
c
c     i=1,im=1,(ss)-loop. x0=j0(ss),y0=k0(ss)
c.......................................................................
 30                     x0 = prfac4
                        y0 = xfsum
                        pj = x0 - y0*ss0
                        qj = -ajmn(1)*y0
                        go to 130
c.......................................................................
c
c     i=2,im=1,(sp)-loop. x0=j0(sp),y1=k1(sp)
c.......................................................................
 40                     x0 = prfac4*(3*zpq+2*zrs)
                        y1 = xfsum
                        pj = x0 - y1*sp1
                        qj = -ajmn(2)*y1
                        go to 130
c.......................................................................
c
c     i=2,im=2,(pp)-loop. x0=j0(pp),y0=k0(pp),y2=k2(pp)
c.......................................................................
 50                     x0 = prfac4*(zpqrs2+zpq*zrs)
                        y0 = xfsum*zpqrs2
                        xfsum = xfac11*zpr*zqs + xfac21*zps*zqr
                        y0 = y0 + xfsum
                        y2 = xfsum*5
                        pj = x0 - y0*pp0 - y2*pp2
                        qj = -(ajmn(3)*y0+ajmn(4)*y2)
                        go to 130
c.......................................................................
c
c     i=3,im=1,(sd)-loop. x0=j0(sd),y2=k2(sd)
c.......................................................................
 60                     x0 = prfac4*(15*zpq**2+20*zpq*zrs+8*zrs**2)
                        y2 = xfsum
                        pj = x0 - y2*sd2
                        qj = -ajmn(5)*y2
                        go to 130
c.......................................................................
c
c     i=3,im=2,(pd)-loop. x0=j0(pd),y1=k1(pd),y3=k3(pd),x2=j2(pd)
c.......................................................................
 70                     x0 = prfac4*(10*zpq**3+35*zpq**2*zrs+
     +                       28*zpq*zrs**2+8*zrs**3)
                        y1 = xfsum*zpqrs2
                        xfsum = xfac11*zpr*zqs + xfac21*zps*zqr
                        y1 = y1 + 3*xfsum
                        y3 = 7*xfsum
                        pj = x0 - y1*pd1 - y3*pd3
                        qj = 0.d0
                        if (open) then
                           x2 = prfac4*5*zpq*zrs*(7*zpq+2*zrs)
                           qj = ajmn(21)*x2 - (ajmn(6)*y1+ajmn(7)*y3)
                        end if
                        go to 130
c.......................................................................
c
c     i=3,im=3,(dd)-loop. x0=j0(dd),y0=k0(dd),y2=k2(dd),y4=k4(dd)
c.......................................................................
 80                     zprzqs = zpr*zqs
                        zpszqr = zps*zqr
                        x0 = prfac4*((zpqrs2+zpq*zrs)
     +                       *zpqrs2*2+7*(zpq*zrs)**2)
                        y01 = xfac11*((zpqrs2+zprzqs)
     +                        *zpqrs2*2+7*zprzqs**2)
                        y02 = xfac21*((zpqrs2+zpszqr)
     +                        *zpqrs2*2+7*zpszqr**2)
                        y0 = y01 + y02
                        xfac11 = xfac11*7*zprzqs
                        xfac21 = xfac21*7*zpszqr
                        y21 = xfac11*(zpqrs2+5*zprzqs)
                        y22 = xfac21*(zpqrs2+5*zpszqr)
                        y2 = y21 + y22
                        y4 = (xfac11*zprzqs+xfac21*zpszqr)*9
                        pj = x0 - y0*dd0 - (y2+y4)*dd2
                        qj = -ajmn(8)*y0 - ajmn(9)*y2 - ajmn(10)*y4
                        go to 130
c.......................................................................
c
c     i=4,im=1,(sf)-loop. x0=j0(sf),y3=k3(sf)
c.......................................................................
 90                     x0 = prfac4*(105*zpq**3+210*zpq**2*zrs+
     +                       168*zpq*zrs**2+48*zrs**3)
                        y3 = xfsum
                        pj = x0 - y3*sf3
                        qj = -ajmn(11)*y3
                        go to 130
c.......................................................................
c
c     i=4,im=2,(pf)-loop. x0=j0(pf),y2=k2(pf),y4=k4(pf),x2=j2(pf)
c.......................................................................
 100                    x0 = prfac4*(70*zpq**4+315*zpq**3*zrs+
     +                       378*(zpq*zrs)**2+216*zpq*zrs**3+48*zrs**4)
c.......................................................................
c                       y2 = xfsum*zpqrs2
c.......................................................................
                        y2 = xfac11*(2*zpr**2+9*zpr*zqs+2*zqs**2)
     +                       + xfac21*(2*zps**2+9*zps*zqr+2*zqr**2)
                        xfsum = xfac11*zpr*zqs + xfac21*zps*zqr
c.......................................................................
c                       y2 = y2 + 5*xfsum
c.......................................................................
                        y4 = 9*xfsum
                        pj = x0 - y2*pf2 - y4*pf4
                        qj = 0.0d0
                        if (open) then
                           x2 = prfac4*5*zpq*zrs*(63*zpq**2+26*zpq*zrs+
     +                          8*zrs**2)
                           qj = ajmn(22)*x2 - ajmn(12)*y2 - ajmn(13)*y4
                        end if
                        go to 130
c.......................................................................
c
c     i=4,im=3,(df)-loop. x0=j0(df),y1=k1(df),y3=k3(df),y5=k5(df)
c                         x2=j2(df),x4=j4(df).
c.......................................................................
 110                    x0 = prfac4*(56*zpq**5+308*zpq**4*zrs+
     +                       693*zpq**3*zrs**2+594*zpq**2*zrs**3+
     +                       264*zpq*zrs**4+48*zrs**5)
                        y1 = xfac11*df1pol(zpr,zqs)
     +                       + xfac21*df1pol(zps,zqr)
                        xfac11 = 9*zpr*zqs*xfac11
                        xfac21 = 9*zps*zqr*xfac21
                        y3 = xfac11*df3pol(zpr,zqs)
     +                       + xfac21*df3pol(zps,zqr)
                        y5 = 11*(xfac11*zpr*zqs+xfac21*zps*zqr)
                        pj = x0 - y1*df1 - y3*df3 - y5*df5
                        qj = 0.0d0
                        if (open) then
                           prfac4 = prfac4*7*zpq*zrs
                           x2 = prfac4*(18*zpq**3+99*zpq**2*zrs+
     +                          44*zpq*zrs**2+8*zrs**3)
                           x4 = 9*prfac4*(11*zpq+2*zrs)*zpq*zrs
                           qj = x2*ajmn(23) + x4*ajmn(24) - y1*ajmn(14)
     +                          - y3*ajmn(15) - y5*ajmn(16)
                        end if
                        go to 130
c.......................................................................
c
c     i=4,im=4,(ff)-loop. x0=j0(ff),y0=k0(ff),y2=k2(ff),y4=k4(ff)
c              y6=k6(ff)
c.......................................................................
 120                    x0 = prfac4*f0pol(zpq,zrs)
                        y0 = xfac11*f0pol(zpr,zqs)
     +                       + xfac21*f0pol(zps,zqr)
                        xfac11 = xfac11*zpr*zqs*9
                        xfac21 = xfac21*zps*zqr*9
                        y2 = xfac11*f2pol(zpr,zqs)
     +                       + xfac21*f2pol(zps,zqr)
                        xfac11 = xfac11*zpr*zqs*11
                        xfac21 = xfac21*zps*zqr*11
                        y4 = xfac11*(2*zpr**2+13*zpr*zqs+2*zqs**2)
     +                       + xfac21*(2*zps**2+13*zps*zqr+2*zqr**2)
                        y6 = 13*(xfac11*zpr*zqs+xfac21*zps*zqr)
                        pj = x0 - y0*ff0 - y2*ff2 - y4*ff4 - y6*ff6
                        qj = -y0*ajmn(17) - y2*ajmn(18) - y4*ajmn(19)
     +                       - y6*ajmn(20)
 130                    pcap(kl) = pcap(kl) + pj*dt(mn)
                        if (klnemn) pcap(mn) = pcap(mn) + pj*dt(kl)
                        term = dt(kl)*pj*dt(mn)
                        if (open) then
                           qcap(kl) = qcap(kl) + qj*dos(mn)
                           if (klnemn) qcap(mn) = qcap(mn) + qj*dos(kl)
                           term = term - dos(kl)*qj*dos(mn)
                        end if
                        if (.not.klnemn) term = term*0.5d0
                        pot = pot + term
 140                 continue
                     enddo
                  enddo
                  factmn = factmn/im
               enddo
               potn = potn + u(kl)*dt(kl)
               cin = cin + t(kl)*dt(kl)
            enddo
         enddo
         factkl = factkl/i
      enddo
      pot = pot - zn*potn
      energ = cin + pot
      vir = pot/cin
      return
      end
