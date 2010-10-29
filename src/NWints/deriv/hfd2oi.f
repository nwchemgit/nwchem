      Subroutine hfd2oi(E,Sab,NPP,Nint,La,Lb,Li,MXD,inder,Nder,canAB)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB

c--> Hermite Linear Expansion Coefficients

      Dimension E(3,NPP,0:MXD,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> 2-Center Overlap Integrals

      Dimension Sab(Nint,*)

c--> Derivative Indices

      Dimension inder(3,Nder)

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

c Initialize the block of O2Is.

      do 10 id = 1,Nder
      do 10 nn = 1,Nint
       Sab(nn,id) = 0.D0
   10 continue

c Define the number of shell components on each center.

      La2 = ((La+1)*(La+2))/2
      Lb2 = ((Lb+1)*(Lb+2))/2

c Loop over shell components.

      mm = 0

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

        mm = mm + 1

        do 25 id = 1,Nder

         n1 = inder(1,id)
         n2 = inder(2,id)
         n3 = inder(3,id)

         do 20 mp = 1,NPP
          Sab(mm,id) = Sab(mm,id) + E(1,mp,n1,0,Ia,Ib)*
     &                              E(2,mp,n2,0,Ja,Jb)*
     &                              E(3,mp,n3,0,Ka,Kb)
   20    continue

   25   continue

   40  continue

   50 continue

      end
