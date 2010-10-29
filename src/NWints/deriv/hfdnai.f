      Subroutine hfdnai(E,R0,IJK,Vab,NPP,Nint,La,Lb,Li,Lp,Lp3,ncenters,
     &                  MXD,inder,Nder,canAB)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB

c--> Hermite Linear Expansion Coefficients

      Dimension E(3,NPP,0:MXD,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> Auxiliary Function Integrals & Index

      Dimension R0(NPP,Lp3),IJK(0:Lp,0:Lp,0:Lp)

c--> Nuclear Attraction Integrals

      Dimension Vab(NPP,Nint,Nder)

c--> Derivative Index

      Dimension inder(6,Nder)

c--> Scratch Space

      Dimension Nxyz(3)
c
c Compute the nuclear attraction integrals.
c
c     formula:
c           __
c           \    Ia,Ib    Ja,Jb    Ka,Kb
c     Vab = /  Ex     * Ey     * Ez     * R
c           --   Ip       Jp       Kp      Ip,Jp,Kp
c        Ip,Jp,Kp
c
c******************************************************************************

c Initialize the block of derivative NAIs.

      do 10 id = 1,Nder
      do 10 nn = 1,Nint
      do 10 mp = 1,NPP
       Vab(mp,nn,id) = 0.D0
   10 continue

c Define the number of shell components on each center.

      La2 = ((La+1)*(La+2))/2
      Lb2 = ((Lb+1)*(Lb+2))/2

c Loop over shell components.

      nn = 0

      do 50 ma = 1,La2

c Define the angular momentum indices for shell "A".

       call getNxyz(La,ma,Nxyz)

       Ia = Nxyz(1)
       Ja = Nxyz(2)
       Ka = Nxyz(3)

       if( canAB )then
        mb_limit = ma
       else
        mb_limit = Lb2
       end if

       do 40 mb = 1,mb_limit

c Define the angular momentum indices for shell "B".

        call getNxyz(Lb,mb,Nxyz)

        Ib = Nxyz(1)
        Jb = Nxyz(2)
        Kb = Nxyz(3)

        nn = nn + 1

        do 30 Ip = 0,Ia+Ib
        do 30 Jp = 0,Ja+Jb
        do 30 Kp = 0,Ka+Kb

         do 25 id = 1,Nder

          n1 = inder(1,id)
          n2 = inder(2,id)
          n3 = inder(3,id)
          n4 = inder(4,id)
          n5 = inder(5,id)
          n6 = inder(6,id)

          np = IJK(Ip+n1,Jp+n2,Kp+n3)

          do 20 mp = 1,NPP
           Vab(mp,nn,id) = Vab(mp,nn,id) + E(1,mp,n4,Ip,Ia,Ib)*
     &                                     E(2,mp,n5,Jp,Ja,Jb)*
     &                                     E(3,mp,n6,Kp,Ka,Kb)*
     &                                     R0(mp,np)
   20    continue

   25    continue

   30   continue

   40  continue

   50 continue

      end
