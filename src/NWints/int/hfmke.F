      Subroutine hfmke(Axyz,Bxyz,alpha,ES,E,pf,nd,NPP,MXD,La,Lb)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

c--> Cartesian Coordinates

      Dimension Axyz(3),Bxyz(3)

c--> Exponents & Scaling Factors

      Dimension alpha(2,NPP),ES(3,NPP)

c--> Hermite Linear Expansion Coefficients

      Dimension E(3,NPP,0:MXD,0:(La+Lb),0:La,0:Lb)

c--> Scratch Space

      Dimension pf(2,NPP)
c
c Define the linear expansion coefficients for the product of two CGF in
c terms of HGF.
c
c   Recursion Formulas:
c
c      Ia+1,Ib;n            Ia,Ib;n            Ia,Ib;n           Ia,Ib;n
c     E           = (1/2p) E        - (b/p) R E        + (Ip+1) E
c      Ip                   Ip-1               Ip                Ip+1
c
c      Ia,Ib+1;n            Ia,Ib;n            Ia,Ib;n           Ia,Ib;n
c     E           = (1/2p) E        + (a/p) R E        + (Ip+1) E
c      Ip                   Ip-1               Ip                Ip+1
c
c   Initial Values:
c
c       0,0;0            0,0;0            0,0;0
c     Ex      = ESx,   Ey      = ESy,   Ez      = ESz,
c       0                0                0
c
c                               /  PI   \ 1/2      /   a b     2 \
c      where typically   ESx = | ------- |     EXP| - -----  Rx   |
c                               \ a + b /          \  a + b      /
c
c N.B. The prefactors for the overlap distribution of the two CGF 
c are typically used to scale the expansion coefficients. Products of 
c contraction coefficients may also be incorporated.
c
c   Indices for E(k,m,n,Ip,Ia,Ib):
c
c     1 [  k  ] 3            Cartesian Components {X,Y,Z}
c     1 [  m  ] NPP          CGF Product Pair Index
c     0 [  nd ] MXD          Derivative Index
c     0 [  Ip ] Ia+Ib        Angular Momentum Index for the HGF
c     0 [  Ia ] La           Angular Momentum Index for the CGF on Center "A"
c     0 [  Ib ] Lb           Angular Momentum Index for the CGF on Center "B"
c
c******************************************************************************

c Compute the vector between the centers.

      Rx = Axyz(1) - Bxyz(1)
      Ry = Axyz(2) - Bxyz(2)
      Rz = Axyz(3) - Bxyz(3)

      if( nd.eq.0 )then

c Define the expansion coefficients for the products of two CGF.

c Define E(Ip,Ia,Ib) for Ip=0, Ia=0, Ib=0.

       do 100 m = 1,NPP
        E(1,m,nd,0,0,0) = ES(1,m)
        E(2,m,nd,0,0,0) = ES(2,m)
        E(3,m,nd,0,0,0) = ES(3,m)
  100  continue

c Compute the prefactor for the 1st term of the recursion formulas.

       if( La.gt.0 .or. Lb.gt.0 )then
        do 110 m = 1,NPP
         pf(1,m) = 0.5D0/( alpha(1,m) + alpha(2,m) )
  110   continue
       end if

c Define E(Ip,Ia,Ib) for Ip = 0,Ia+Ib; Ia = 1,La; Ib = 0.

       if( La.gt.0 )then

c Compute the prefactor for the 2nd term of the recursion formula.

        do 200 m = 1,NPP
         pf(2,m) = -2.D0*alpha(2,m)*pf(1,m)
  200   continue

c            ===>   Ip = 0,Ia; Ia = 1; Ib = 0   <===

        do 210 m = 1,NPP

         Ex2 = Rx*E(1,m,nd,0,0,0)
         Ey2 = Ry*E(2,m,nd,0,0,0)
         Ez2 = Rz*E(3,m,nd,0,0,0)

         E(1,m,nd,0,1,0) = pf(2,m)*Ex2
         E(2,m,nd,0,1,0) = pf(2,m)*Ey2
         E(3,m,nd,0,1,0) = pf(2,m)*Ez2

         Ex1 = E(1,m,nd,0,0,0)
         Ey1 = E(2,m,nd,0,0,0)
         Ez1 = E(3,m,nd,0,0,0)

         E(1,m,nd,1,1,0) = pf(1,m)*Ex1
         E(2,m,nd,1,1,0) = pf(1,m)*Ey1
         E(3,m,nd,1,1,0) = pf(1,m)*Ez1

  210   continue

        do 260 Ia = 2,La

c            ===>   Ip = 0; Ia = 2,La; Ib = 0   <===

         do 220 m = 1,NPP

          Ex2 = Rx*E(1,m,nd,0,Ia-1,0)
          Ex3 =    E(1,m,nd,1,Ia-1,0)

          Ey2 = Ry*E(2,m,nd,0,Ia-1,0)
          Ey3 =    E(2,m,nd,1,Ia-1,0)

          Ez2 = Rz*E(3,m,nd,0,Ia-1,0)
          Ez3 =    E(3,m,nd,1,Ia-1,0)

          E(1,m,nd,0,Ia,0) = pf(2,m)*Ex2 + Ex3
          E(2,m,nd,0,Ia,0) = pf(2,m)*Ey2 + Ey3
          E(3,m,nd,0,Ia,0) = pf(2,m)*Ez2 + Ez3

  220    continue

c            ===>   Ip = 1,Ia-2; Ia = 2,La; Ib = 0   <===

          do 240 Ip = 1,Ia-2

           do m = 1,NPP

            Ex1 =    E(1,m,nd,Ip-1,Ia-1,0)
            Ex2 = Rx*E(1,m,nd,Ip  ,Ia-1,0)
            Ex3 =    E(1,m,nd,Ip+1,Ia-1,0)

            Ey1 =    E(2,m,nd,Ip-1,Ia-1,0)
            Ey2 = Ry*E(2,m,nd,Ip  ,Ia-1,0)
            Ey3 =    E(2,m,nd,Ip+1,Ia-1,0)

            Ez1 =    E(3,m,nd,Ip-1,Ia-1,0)
            Ez2 = Rz*E(3,m,nd,Ip  ,Ia-1,0)
            Ez3 =    E(3,m,nd,Ip+1,Ia-1,0)

            E(1,m,nd,Ip,Ia,0) = pf(1,m)*Ex1 + pf(2,m)*Ex2 + (Ip+1)*Ex3
            E(2,m,nd,Ip,Ia,0) = pf(1,m)*Ey1 + pf(2,m)*Ey2 + (Ip+1)*Ey3
            E(3,m,nd,Ip,Ia,0) = pf(1,m)*Ez1 + pf(2,m)*Ez2 + (Ip+1)*Ez3

         enddo

  240    continue

c            ===>   Ip = Ia-1,Ia; Ia = 2,La; Ib = 0   <===

         Ip = Ia-1
cbreakakakakaka
         do  m = 1,NPP
            Ex1=E(1,m,nd,Ip,Ia-1,0)
            Ey1=E(2,m,nd,Ip,Ia-1,0)
            Ez1=E(3,m,nd,Ip,Ia-1,0)
            E(1,m,nd,Ip+1,Ia,0) = pf(1,m)*Ex1
            E(2,m,nd,Ip+1,Ia,0) = pf(1,m)*Ey1
            E(3,m,nd,Ip+1,Ia,0) = pf(1,m)*Ez1
            E(1,m,nd,Ip,Ia,0) =   pf(2,m)*Rx*Ex1
            E(2,m,nd,Ip,Ia,0) =   pf(2,m)*Ry*Ey1
            E(3,m,nd,Ip,Ia,0) =   pf(2,m)*Rz*Ez1
            
         enddo

         do  m = 1,NPP
          E(1,m,nd,Ip,Ia,0) = E(1,m,nd,Ip,Ia,0) + 
     +           pf(1,m)*E(1,m,nd,Ip-1,Ia-1,0) 
          E(2,m,nd,Ip,Ia,0) = E(2,m,nd,Ip,Ia,0) +
     +         pf(1,m)*E(2,m,nd,Ip-1,Ia-1,0) 
          E(3,m,nd,Ip,Ia,0) = E(3,m,nd,Ip,Ia,0) +
     +         pf(1,m)*E(3,m,nd,Ip-1,Ia-1,0)  


        enddo

  260   continue

       end if

c Define E(Ip,Ia,Ib) for Ip=0,Ia+Ib, Ia=0,La, Ib=1,Lb.

       if( Lb.gt.0 )then

c Compute the prefactor for the 2nd term of the recursion formula.

        do 300 m = 1,NPP
         pf(2,m) = 2.D0*alpha(1,m)*pf(1,m)
  300   continue

c    ===>   Ip = 0,Ia+Ib; Ia = 0; Ib = 1   <===

        do 310 m = 1,NPP

         Ex2 = Rx*E(1,m,nd,0,0,0)
         Ey2 = Ry*E(2,m,nd,0,0,0)
         Ez2 = Rz*E(3,m,nd,0,0,0)

         E(1,m,nd,0,0,1) = pf(2,m)*Ex2
         E(2,m,nd,0,0,1) = pf(2,m)*Ey2
         E(3,m,nd,0,0,1) = pf(2,m)*Ez2

         Ex1 = E(1,m,nd,0,0,0)
         Ey1 = E(2,m,nd,0,0,0)
         Ez1 = E(3,m,nd,0,0,0)

         E(1,m,nd,1,0,1) = pf(1,m)*Ex1
         E(2,m,nd,1,0,1) = pf(1,m)*Ey1
         E(3,m,nd,1,0,1) = pf(1,m)*Ez1

  310   continue

        do 370 Ib = 1,Lb

         if( Ib.eq.1 )then
          Ia1 = 1
         else
          Ia1 = 0
         end if

         do 360 Ia = Ia1,La

c    ===>   Ip = 0; Ia = Ia1,La; Ib = 1,Lb   <===

          do 320 m = 1,NPP

           Ex2 = Rx*E(1,m,nd,0,Ia,Ib-1)
           Ex3 =    E(1,m,nd,1,Ia,Ib-1)

           Ey2 = Ry*E(2,m,nd,0,Ia,Ib-1)
           Ey3 =    E(2,m,nd,1,Ia,Ib-1)

           Ez2 = Rz*E(3,m,nd,0,Ia,Ib-1)
           Ez3 =    E(3,m,nd,1,Ia,Ib-1)

           E(1,m,nd,0,Ia,Ib) = pf(2,m)*Ex2 + Ex3
           E(2,m,nd,0,Ia,Ib) = pf(2,m)*Ey2 + Ey3
           E(3,m,nd,0,Ia,Ib) = pf(2,m)*Ez2 + Ez3

  320     continue

c    ===>   Ip = 1,Ia+Ib-2; Ia = Ia1,La; Ib = 1,Lb   <===

          do 340 Ip = 1,Ia+Ib-2

           do 330 m = 1,NPP

            Ex1 =    E(1,m,nd,Ip-1,Ia,Ib-1)
            Ex2 = Rx*E(1,m,nd,Ip  ,Ia,Ib-1)
            Ex3 =    E(1,m,nd,Ip+1,Ia,Ib-1)

            Ey1 =    E(2,m,nd,Ip-1,Ia,Ib-1)
            Ey2 = Ry*E(2,m,nd,Ip  ,Ia,Ib-1)
            Ey3 =    E(2,m,nd,Ip+1,Ia,Ib-1)

            Ez1 =    E(3,m,nd,Ip-1,Ia,Ib-1)
            Ez2 = Rz*E(3,m,nd,Ip  ,Ia,Ib-1)
            Ez3 =    E(3,m,nd,Ip+1,Ia,Ib-1)

            E(1,m,nd,Ip,Ia,Ib) = pf(1,m)*Ex1 + pf(2,m)*Ex2 + (Ip+1)*Ex3
            E(2,m,nd,Ip,Ia,Ib) = pf(1,m)*Ey1 + pf(2,m)*Ey2 + (Ip+1)*Ey3
            E(3,m,nd,Ip,Ia,Ib) = pf(1,m)*Ez1 + pf(2,m)*Ez2 + (Ip+1)*Ez3

  330      continue

  340     continue

c    ===>   Ip = Ia+Ib-1,Ia+Ib; Ia = Ia1,La; Ib = 1,Lb   <===

          Ip = Ia+Ib-1

          do 350 m = 1,NPP

           Ex1 =    E(1,m,nd,Ip-1,Ia,Ib-1)
           Ex2 = Rx*E(1,m,nd,Ip  ,Ia,Ib-1)

           Ey1 =    E(2,m,nd,Ip-1,Ia,Ib-1)
           Ey2 = Ry*E(2,m,nd,Ip  ,Ia,Ib-1)

           Ez1 =    E(3,m,nd,Ip-1,Ia,Ib-1)
           Ez2 = Rz*E(3,m,nd,Ip  ,Ia,Ib-1)

           E(1,m,nd,Ip,Ia,Ib) = pf(1,m)*Ex1 + pf(2,m)*Ex2
           E(2,m,nd,Ip,Ia,Ib) = pf(1,m)*Ey1 + pf(2,m)*Ey2
           E(3,m,nd,Ip,Ia,Ib) = pf(1,m)*Ez1 + pf(2,m)*Ez2

           Ex1 = E(1,m,nd,Ip,Ia,Ib-1)
           Ey1 = E(2,m,nd,Ip,Ia,Ib-1)
           Ez1 = E(3,m,nd,Ip,Ia,Ib-1)

           E(1,m,nd,Ip+1,Ia,Ib) = pf(1,m)*Ex1
           E(2,m,nd,Ip+1,Ia,Ib) = pf(1,m)*Ey1
           E(3,m,nd,Ip+1,Ia,Ib) = pf(1,m)*Ez1

  350     continue

  360    continue

  370   continue

       end if

      else

c Define the expansion coefficients for derivatives of the products of two CGF.

c Define E(Ip,Ia,Ib) for Ip=0, Ia=0, Ib=0.

       if( nd.eq.1 )then

        do 1100 m = 1,NPP

         c = -2.D0*( alpha(1,m)*alpha(2,m)/( alpha(1,m) + alpha(2,m) ) )

         E(1,m,nd,0,0,0) = c*Rx*E(1,m,nd-1,0,0,0)
         E(2,m,nd,0,0,0) = c*Ry*E(2,m,nd-1,0,0,0)
         E(3,m,nd,0,0,0) = c*Rz*E(3,m,nd-1,0,0,0)
 
 1100   continue

       else

        do 1105 m = 1,NPP

         c = -2.D0*( alpha(1,m)*alpha(2,m)/( alpha(1,m) + alpha(2,m) ) )

         E1 = c*(Rx*E(1,m,nd-1,0,0,0) + (nd-1)*E(1,m,nd-2,0,0,0))
         E2 = c*(Ry*E(2,m,nd-1,0,0,0) + (nd-1)*E(2,m,nd-2,0,0,0))
         E3 = c*(Rz*E(3,m,nd-1,0,0,0) + (nd-1)*E(3,m,nd-2,0,0,0))

         E(1,m,nd,0,0,0) = E1
         E(2,m,nd,0,0,0) = E2
         E(3,m,nd,0,0,0) = E3

 1105   continue

       end if

c Compute the prefactor for the 1st term of recursion formulas.

       if( La.gt.0 .or. Lb.gt.0 )then
        do 1110 m = 1,NPP
         pf(1,m) = 0.5D0/( alpha(1,m) + alpha(2,m) )
 1110   continue
       end if

c Define E(Ip,Ia,Ib) for Ip = 0,Ia+Ib; Ia = 1,La; Ib = 0.

       if( La.gt.0 )then

c Compute the prefactor for the 2nd term of the recursion formula.

        do 1200 m = 1,NPP
         pf(2,m) = -2.D0*alpha(2,m)*pf(1,m)
 1200   continue

c            ===>   Ip = 0,Ia; Ia = 1; Ib = 0   <===

        do 1210 m = 1,NPP

         Ex2 = Rx*E(1,m,nd,0,0,0) + nd*E(1,m,nd-1,0,0,0)
         Ey2 = Ry*E(2,m,nd,0,0,0) + nd*E(2,m,nd-1,0,0,0)
         Ez2 = Rz*E(3,m,nd,0,0,0) + nd*E(3,m,nd-1,0,0,0)

         E(1,m,nd,0,1,0) = pf(2,m)*Ex2
         E(2,m,nd,0,1,0) = pf(2,m)*Ey2
         E(3,m,nd,0,1,0) = pf(2,m)*Ez2

         Ex1 = E(1,m,nd,0,0,0)
         Ey1 = E(2,m,nd,0,0,0)
         Ez1 = E(3,m,nd,0,0,0)

         E(1,m,nd,1,1,0) = pf(1,m)*Ex1
         E(2,m,nd,1,1,0) = pf(1,m)*Ey1
         E(3,m,nd,1,1,0) = pf(1,m)*Ez1

 1210   continue

        do 1260 Ia = 2,La

c            ===>   Ip = 0; Ia = 2,La; Ib = 0   <===

         do 1220 m = 1,NPP

          Ex2 = Rx*E(1,m,nd,0,Ia-1,0) + nd*E(1,m,nd-1,0,Ia-1,0)
          Ex3 =    E(1,m,nd,1,Ia-1,0)

          Ey2 = Ry*E(2,m,nd,0,Ia-1,0) + nd*E(2,m,nd-1,0,Ia-1,0)
          Ey3 =    E(2,m,nd,1,Ia-1,0)

          Ez2 = Rz*E(3,m,nd,0,Ia-1,0) + nd*E(3,m,nd-1,0,Ia-1,0)
          Ez3 =    E(3,m,nd,1,Ia-1,0)

          E(1,m,nd,0,Ia,0) = pf(2,m)*Ex2 + Ex3
          E(2,m,nd,0,Ia,0) = pf(2,m)*Ey2 + Ey3
          E(3,m,nd,0,Ia,0) = pf(2,m)*Ez2 + Ez3

 1220    continue

c            ===>   Ip = 1,Ia-2; Ia = 2,La; Ib = 0   <===

          do 1240 Ip = 1,Ia-2

           do 1230 m = 1,NPP

            Ex1 =    E(1,m,nd,Ip-1,Ia-1,0)
            Ex2 = Rx*E(1,m,nd,Ip  ,Ia-1,0) + nd*E(1,m,nd-1,Ip,Ia-1,0)
            Ex3 =    E(1,m,nd,Ip+1,Ia-1,0)

            Ey1 =    E(2,m,nd,Ip-1,Ia-1,0)
            Ey2 = Ry*E(2,m,nd,Ip  ,Ia-1,0) + nd*E(2,m,nd-1,Ip,Ia-1,0)
            Ey3 =    E(2,m,nd,Ip+1,Ia-1,0)

            Ez1 =    E(3,m,nd,Ip-1,Ia-1,0)
            Ez2 = Rz*E(3,m,nd,Ip  ,Ia-1,0) + nd*E(3,m,nd-1,Ip,Ia-1,0)
            Ez3 =    E(3,m,nd,Ip+1,Ia-1,0)

            E(1,m,nd,Ip,Ia,0) = pf(1,m)*Ex1 + pf(2,m)*Ex2 + (Ip+1)*Ex3
            E(2,m,nd,Ip,Ia,0) = pf(1,m)*Ey1 + pf(2,m)*Ey2 + (Ip+1)*Ey3
            E(3,m,nd,Ip,Ia,0) = pf(1,m)*Ez1 + pf(2,m)*Ez2 + (Ip+1)*Ez3

 1230     continue

 1240    continue

c            ===>   Ip = Ia-1,Ia; Ia = 2,La; Ib = 0   <===

         Ip = Ia-1

cbreak2
         do m = 1,NPP

          Ex1 = E(1,m,nd,Ip,Ia-1,0)
          Ey1 = E(2,m,nd,Ip,Ia-1,0)
          Ez1 = E(3,m,nd,Ip,Ia-1,0)

          E(1,m,nd,Ip,Ia,0) =  pf(2,m) * (Rx*Ex1 + 
     +         nd*E(1,m,nd-1,Ip,Ia-1,0))
          E(2,m,nd,Ip,Ia,0) =  pf(2,m) * (Ry*Ey1 + 
     +         nd*E(2,m,nd-1,Ip,Ia-1,0))
          E(3,m,nd,Ip,Ia,0) =  pf(2,m) * (Rz*Ez1 + 
     +         nd*E(3,m,nd-1,Ip,Ia-1,0))


          E(1,m,nd,Ip+1,Ia,0) = pf(1,m)*Ex1
          E(2,m,nd,Ip+1,Ia,0) = pf(1,m)*Ey1
          E(3,m,nd,Ip+1,Ia,0) = pf(1,m)*Ez1

       enddo
         do  m = 1,NPP

            Ex1 =    E(1,m,nd,Ip-1,Ia-1,0) 
            Ey1 =    E(2,m,nd,Ip-1,Ia-1,0) 
            Ez1 =    E(3,m,nd,Ip-1,Ia-1,0) 
            

            E(1,m,nd,Ip,Ia,0) = E(1,m,nd,Ip,Ia,0)+pf(1,m)*Ex1 
            E(2,m,nd,Ip,Ia,0) = E(2,m,nd,Ip,Ia,0)+pf(1,m)*Ey1 
            E(3,m,nd,Ip,Ia,0) = E(3,m,nd,Ip,Ia,0)+pf(1,m)*Ez1 
            

         enddo

 1260   continue

       end if

c Define E(Ip,Ia,Ib) for Ip=0,Ia+Ib, Ia=0,La, Ib=1,Lb.

       if( Lb.gt.0 )then

c Compute the prefactor for the 2nd term of the recursion formula.

        do 1300 m = 1,NPP
         pf(2,m) = 2.D0*alpha(1,m)*pf(1,m)
 1300   continue

c    ===>   Ip = 0,Ia+Ib; Ia = 0; Ib = 1   <===

        do 1310 m = 1,NPP

         Ex2 = Rx*E(1,m,nd,0,0,0) + nd*E(1,m,nd-1,0,0,0)
         Ey2 = Ry*E(2,m,nd,0,0,0) + nd*E(2,m,nd-1,0,0,0)
         Ez2 = Rz*E(3,m,nd,0,0,0) + nd*E(3,m,nd-1,0,0,0)

         E(1,m,nd,0,0,1) = pf(2,m)*Ex2
         E(2,m,nd,0,0,1) = pf(2,m)*Ey2
         E(3,m,nd,0,0,1) = pf(2,m)*Ez2

         Ex1 = E(1,m,nd,0,0,0)
         Ey1 = E(2,m,nd,0,0,0)
         Ez1 = E(3,m,nd,0,0,0)

         E(1,m,nd,1,0,1) = pf(1,m)*Ex1
         E(2,m,nd,1,0,1) = pf(1,m)*Ey1
         E(3,m,nd,1,0,1) = pf(1,m)*Ez1

 1310   continue

        do 1370 Ib = 1,Lb

         if( Ib.eq.1 )then
          Ia1 = 1
         else
          Ia1 = 0
         end if

         do 1360 Ia = Ia1,La

c    ===>   Ip = 0; Ia = Ia1,La; Ib = 1,Lb   <===

          do 1320 m = 1,NPP

           Ex2 = Rx*E(1,m,nd,0,Ia,Ib-1) + nd*E(1,m,nd-1,0,Ia,Ib-1)
           Ex3 =    E(1,m,nd,1,Ia,Ib-1)

           Ey2 = Ry*E(2,m,nd,0,Ia,Ib-1) + nd*E(2,m,nd-1,0,Ia,Ib-1)
           Ey3 =    E(2,m,nd,1,Ia,Ib-1)

           Ez2 = Rz*E(3,m,nd,0,Ia,Ib-1) + nd*E(3,m,nd-1,0,Ia,Ib-1)
           Ez3 =    E(3,m,nd,1,Ia,Ib-1)

           E(1,m,nd,0,Ia,Ib) = pf(2,m)*Ex2 + Ex3
           E(2,m,nd,0,Ia,Ib) = pf(2,m)*Ey2 + Ey3
           E(3,m,nd,0,Ia,Ib) = pf(2,m)*Ez2 + Ez3

 1320     continue

c    ===>   Ip = 1,Ia+Ib-2; Ia = Ia1,La; Ib = 1,Lb   <===

          do 1340 Ip = 1,Ia+Ib-2

           do 1330 m = 1,NPP

            Ex1 =    E(1,m,nd,Ip-1,Ia,Ib-1)
            Ex2 = Rx*E(1,m,nd,Ip  ,Ia,Ib-1) + nd*E(1,m,nd-1,Ip,Ia,Ib-1)
            Ex3 =    E(1,m,nd,Ip+1,Ia,Ib-1)

            Ey1 =    E(2,m,nd,Ip-1,Ia,Ib-1)
            Ey2 = Ry*E(2,m,nd,Ip  ,Ia,Ib-1) + nd*E(2,m,nd-1,Ip,Ia,Ib-1)
            Ey3 =    E(2,m,nd,Ip+1,Ia,Ib-1)

            Ez1 =    E(3,m,nd,Ip-1,Ia,Ib-1)
            Ez2 = Rz*E(3,m,nd,Ip  ,Ia,Ib-1) + nd*E(3,m,nd-1,Ip,Ia,Ib-1)
            Ez3 =    E(3,m,nd,Ip+1,Ia,Ib-1)

            E(1,m,nd,Ip,Ia,Ib) = pf(1,m)*Ex1 + pf(2,m)*Ex2 + (Ip+1)*Ex3
            E(2,m,nd,Ip,Ia,Ib) = pf(1,m)*Ey1 + pf(2,m)*Ey2 + (Ip+1)*Ey3
            E(3,m,nd,Ip,Ia,Ib) = pf(1,m)*Ez1 + pf(2,m)*Ez2 + (Ip+1)*Ez3

 1330      continue

 1340     continue

c    ===>   Ip = Ia+Ib-1,Ia+Ib; Ia = Ia1,La; Ib = 1,Lb   <===

          Ip = Ia+Ib-1

          do 1350 m = 1,NPP

           Ex1 =    E(1,m,nd,Ip-1,Ia,Ib-1)
           Ex2 = Rx*E(1,m,nd,Ip  ,Ia,Ib-1) + nd*E(1,m,nd-1,Ip,Ia,Ib-1)

           Ey1 =    E(2,m,nd,Ip-1,Ia,Ib-1)
           Ey2 = Ry*E(2,m,nd,Ip  ,Ia,Ib-1) + nd*E(2,m,nd-1,Ip,Ia,Ib-1)

           Ez1 =    E(3,m,nd,Ip-1,Ia,Ib-1)
           Ez2 = Rz*E(3,m,nd,Ip  ,Ia,Ib-1) + nd*E(3,m,nd-1,Ip,Ia,Ib-1)

           E(1,m,nd,Ip,Ia,Ib) = pf(1,m)*Ex1 + pf(2,m)*Ex2
           E(2,m,nd,Ip,Ia,Ib) = pf(1,m)*Ey1 + pf(2,m)*Ey2
           E(3,m,nd,Ip,Ia,Ib) = pf(1,m)*Ez1 + pf(2,m)*Ez2

           Ex1 = E(1,m,nd,Ip,Ia,Ib-1)
           Ey1 = E(2,m,nd,Ip,Ia,Ib-1)
           Ez1 = E(3,m,nd,Ip,Ia,Ib-1)

           E(1,m,nd,Ip+1,Ia,Ib) = pf(1,m)*Ex1
           E(2,m,nd,Ip+1,Ia,Ib) = pf(1,m)*Ey1
           E(3,m,nd,Ip+1,Ia,Ib) = pf(1,m)*Ez1

 1350     continue

 1360    continue

 1370   continue

       end if

      end if

      end
