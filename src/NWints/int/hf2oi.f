      Subroutine hf2oi(E,Sab,Nint,NPP,La,Lb,Li,canAB)
c $Id: hf2oi.f,v 1.2 1994-04-04 20:31:01 d3e129 Exp $

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB

c--> Hermite Linear Expansion Coefficients

      Dimension E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> 2-Center Overlap Integrals

      Dimension Sab(Nint)

      Dimension Nxyz(3)
c
c Compute the 2-ctr overlap integrals.
c
c     formula:
c
c             Ia,Ib   Ja,Jb   Ka,Kb
c     Sab = Ex      Ey      Ez
c             0       0       0
c
c******************************************************************************

c Initialize the block of NAIs.

      do 10 nn = 1,Nint
       Sab(nn) = 0.D0
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

        do 20 mp = 1,NPP
         Sab(nn) = Sab(nn) + E(1,mp,0,Ia,Ib)*
     &                       E(2,mp,0,Ja,Jb)*
     &                       E(3,mp,0,Ka,Kb)
   20   continue

   40  continue

   50 continue

      end
