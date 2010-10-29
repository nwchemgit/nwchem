C$PRAGMA SUN OPT=2
      Subroutine hf1mke3(Axyz,Bxyz,Cxyz,alpha,G,GT,ABC2I,E,
     &                   NABC,La,Lb,Lc)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Dimension Axyz(3),Bxyz(3),Cxyz(3),alpha(4,NABC),G(3,NABC)

      Dimension GT(NABC,3),ABC2I(3*NABC)

      Dimension E(NABC,3,0:(La+Lb+Lc),0:La,0:Lb,0:Lc)
c
c Define the Hermite linear expansion coefficients for the product of 
c three monomials. These coefficients are scaled by a factor appropriate
c for 3-ctr OIs, defined as
c
c         /    PI     \ 3/2     /   a b   __2 \      /   (a + b) c   __2 \
c   ES = | ----------- |    EXP| - -----  AB   | EXP| - -----------  PC   |
c         \ a + b + c /         \  a + b      /      \   a + b + c       /
c
c This scaling factor is passed as the 4th index of the exponents array
c (i.e., "alpha(4,m)").
c
c     Recursion Formulae:
c
c       0,0,0
c     Ex      = ES
c       0
c
c       0,0,0
c     Ey      = 1.0
c       0
c
c       0,0,0
c     Ez      = 1.0
c       0
c
c       Ia+1,Ib,Ic           Ia,Ib,Ic         Ia,Ib,Ic            Ia,Ib,Ic
c     Ex           = ABC2I Ex         + GAx Ex         + (Ip+1) Ex
c       Ip                   Ip-1             Ip                  Ip+1
c
c       Ia,Ib+1,Ic           Ia,Ib,Ic         Ia,Ib,Ic            Ia,Ib,Ic
c     Ex           = ABC2I Ex         + GBx Ex         + (Ip+1) Ex
c       Ip                   Ip-1             Ip                  Ip+1
c
c       Ia,Ib,Ic+1           Ia,Ib,Ic         Ia,Ib,Ic            Ia,Ib,Ic
c     Ex           = ABC2I Ex         + GCx Ex         + (Ip+1) Ex
c       Ip                   Ip-1             Ip                  Ip+1
c
c     Indices for E(m,k,Ip,Ia,Ib,Ic):
c
c     1 [  m  ] NABC
c     1 [  k  ] 3     {X,Y,Z}
c     0 [  Ip ] Ia+Ib+Ic
c     0 [  Ia ] La
c     0 [  Ib ] Lb
c     0 [  Ic ] Lc
c
c******************************************************************************

c Initialize the Hermite expansion coefficients.

      isz_e = NABC*3*(La+Lb+Lc+1)*(La+1)*(Lb+1)*(Lc+1)
      call dcopy(isz_e,0d0,0,E,1)
c Define E(Ip,Ia,Ib,Ic) for Ip=0, Ia=0, Ib=0, Ic=0.

      do 100 m = 1,NABC

       E(m,1,0,0,0,0) = alpha(4,m)
       E(m,2,0,0,0,0) = 1.D0
       E(m,3,0,0,0,0) = 1.D0

       ABC2I(m       ) = 0.5D0/(alpha(1,m) + alpha(2,m) + alpha(3,m))
       ABC2I(m+1*NABC) = ABC2I(m)
       ABC2I(m+2*NABC) = ABC2I(m)

       GT(m,1) = G(1,m) - Axyz(1)
       GT(m,2) = G(2,m) - Axyz(2)
       GT(m,3) = G(3,m) - Axyz(3)

  100 continue

c Define E(Ip,Ia,Ib,Ic) for Ip=0,Ia+Ib+Ic, Ia=1,La, Ib=0, Ic=0.

      do 250 Ia = 1,La

c                          ===>   Ip = 0   <===

       do 210 m = 1,3*NABC
        E(m,1,0,Ia,0,0) =   GT(m,1)*E(m,1,0,Ia-1,0,0)
     &                    +         E(m,1,1,Ia-1,0,0)
  210  continue

c                       ===>   Ip = 1,Ia-1   <===

        do 230 Ip = 1,Ia-1

         do 220 m = 1,3*NABC
          E(m,1,Ip,Ia,0,0) =   ABC2I(m)*E(m,1,Ip-1,Ia-1,0,0)
     &                       +  GT(m,1)*E(m,1,Ip  ,Ia-1,0,0)
     &                       +   (Ip+1)*E(m,1,Ip+1,Ia-1,0,0)
  220   continue

  230  continue

c                         ===>   Ip = Ia   <===

       Ip = Ia

       do 240 m = 1,3*NABC
        E(m,1,Ip,Ia,0,0) =   ABC2I(m)*E(m,1,Ip-1,Ia-1,0,0)
     &                     +  GT(m,1)*E(m,1,Ip  ,Ia-1,0,0)
  240  continue

  250 continue

c Define E(Ip,Ia,Ib,Ic) for Ip=0,Ia+Ib+Ic, Ia=0,La, Ib=1,Lb, Ic=0.

      do 300 m = 1,NABC
       GT(m,1) = G(1,m) - Bxyz(1)
       GT(m,2) = G(2,m) - Bxyz(2)
       GT(m,3) = G(3,m) - Bxyz(3)
  300 continue

      do 360 Ib = 1,Lb

       do 350 Ia = 0,La

c                          ===>   Ip = 0   <===

        do 310 m = 1,3*NABC
         E(m,1,0,Ia,Ib,0) =   GT(m,1)*E(m,1,0,Ia,Ib-1,0)
     &                      +         E(m,1,1,Ia,Ib-1,0)
  310   continue

c                       ===>   Ip = 1,Ia+Ib-1   <===

        do 330 Ip = 1,Ia+Ib-1

         do 320 m = 1,3*NABC
          E(m,1,Ip,Ia,Ib,0) =   ABC2I(m)*E(m,1,Ip-1,Ia,Ib-1,0)
     &                        +  GT(m,1)*E(m,1,Ip  ,Ia,Ib-1,0)
     &                        +   (Ip+1)*E(m,1,Ip+1,Ia,Ib-1,0)
  320    continue

  330   continue

c                        ===>   Ip = Ia+Ib   <===

        Ip = Ia + Ib

        do 340 m = 1,3*NABC
         E(m,1,Ip,Ia,Ib,0) =   ABC2I(m)*E(m,1,Ip-1,Ia,Ib-1,0)
     &                       +  GT(m,1)*E(m,1,Ip  ,Ia,Ib-1,0)
  340   continue

  350  continue

  360 continue

c Define E(Ip,Ia,Ib,Ic) for Ip=0,Ia+Ib+Ic, Ia=0,La, Ib=0,Lb, Ic=1,Lc.

      do 400 m = 1,NABC
       GT(m,1) = G(1,m) - Cxyz(1)
       GT(m,2) = G(2,m) - Cxyz(2)
       GT(m,3) = G(3,m) - Cxyz(3)
  400 continue

      do 470 Ic = 1,Lc

       do 460 Ib = 0,Lb

        do 450 Ia = 0,La

c                          ===>   Ip = 0   <===

         do 410 m = 1,3*NABC
          E(m,1,0,Ia,Ib,Ic) =   GT(m,1)*E(m,1,0,Ia,Ib,Ic-1)
     &                        +         E(m,1,1,Ia,Ib,Ic-1)
  410    continue

c                       ===>   Ip = 1,Ia+Ib+Ic-1   <===

         do 430 Ip = 1,Ia+Ib+Ic-1

          do 420 m = 1,3*NABC
           E(m,1,Ip,Ia,Ib,Ic) =   ABC2I(m)*E(m,1,Ip-1,Ia,Ib,Ic-1)
     &                          +  GT(m,1)*E(m,1,Ip  ,Ia,Ib,Ic-1)
     &                          +   (Ip+1)*E(m,1,Ip+1,Ia,Ib,Ic-1)
  420     continue

  430    continue

c                        ===>   Ip = Ia+Ib+Ic   <===

         Ip = Ia + Ib + Ic

         do 440 m = 1,3*NABC
          E(m,1,Ip,Ia,Ib,Ic) =   ABC2I(m)*E(m,1,Ip-1,Ia,Ib,Ic-1)
     &                         +  GT(m,1)*E(m,1,Ip  ,Ia,Ib,Ic-1)
  440    continue

  450   continue

  460  continue

  470 continue

      end
