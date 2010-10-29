      Subroutine hfefi(E,R0C,IJK,bEFI,
     &       NPP,Nint,La,Lb,Li,Lp,Lp3,ncenters,
     &       MXD,canAB,ictrA,ictrB)
c $Id$
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB

c--> Hermite Linear Expansion Coefficients

      Dimension E(3,NPP,0:MXD,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> Auxiliary Function Integrals & Index

      Dimension R0C(ncenters,NPP,Lp3),IJK(0:Lp,0:Lp,0:Lp)

c--> Block of Electric Field Integrals

      Dimension bEFI(Nint,3,ncenters)

c--> Scratch Space

      Dimension Nxyz(3)
c
c Compute the electric field integrals.
c
c     formula:
c               __
c        Cx     \    Ia,Ib    Ja,Jb    Ka,Kb   C
c     EFI     = /  Ex     * Ey     * Ez     * R
c        ab     --   Ip       Jp       Kp      Ip+1,Jp,Kp
c            Ip,Jp,Kp
c
c******************************************************************************

c Initialize the block of EFIs.

      do 10 ic = 1,ncenters
        do 10 nn = 1,Nint
          bEFI(nn,1,ic) = 0.D0
          bEFI(nn,2,ic) = 0.D0
          bEFI(nn,3,ic) = 0.D0
10    continue

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

         npx = IJK(Ip+1,Jp  ,Kp  )
         npy = IJK(Ip  ,Jp+1,Kp  )
         npz = IJK(Ip  ,Jp  ,Kp+1)

         do 20 mp = 1,NPP

          E3 = E(1,mp,0,Ip,Ia,Ib)*
     &         E(2,mp,0,Jp,Ja,Jb)*
     &         E(3,mp,0,Kp,Ka,Kb)

          do 15 ic = 1,ncenters
              bEFI(nn,1,ic) = bEFI(nn,1,ic) - E3*R0C(ic,mp,npx)
              bEFI(nn,2,ic) = bEFI(nn,2,ic) - E3*R0C(ic,mp,npy)
              bEFI(nn,3,ic) = bEFI(nn,3,ic) - E3*R0C(ic,mp,npz)
15        continue
c
*          bEFI(nn,1,ictra) = bEFI(nn,1,ictra) - E3*R0C(ictra,mp,npx)
*          bEFI(nn,2,ictra) = bEFI(nn,2,ictra) - E3*R0C(ictra,mp,npy)
*          bEFI(nn,3,ictra) = bEFI(nn,3,ictra) - E3*R0C(ictra,mp,npz)
*          if (ictra.ne.ictrb) then
*            bEFI(nn,1,ictrb) =  bEFI(nn,1,ictrb) - E3*R0C(ictrb,mp,npx)
*            bEFI(nn,2,ictrb) =  bEFI(nn,2,ictrb) - E3*R0C(ictrb,mp,npy)
*            bEFI(nn,3,ictrb) =  bEFI(nn,3,ictrb) - E3*R0C(ictrb,mp,npz)
*          endif
20      continue

30    continue

40    continue

50    continue

      end
