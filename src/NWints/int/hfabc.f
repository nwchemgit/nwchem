      Subroutine hfabc(E,Sabc,NABC,La,Lb,Lc,TriDiag)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical TriDiag

      Dimension E(NABC,3,0:(La+Lb+Lc),0:La,0:Lb,0:Lc),Sabc(*)

      Dimension Nxyz(3)
c
c Compute a block of 3-ctr OIs.
c
c     formula:
c
c              Ia,Ib,Ic   Ja,Jb,Jc   Ka,Kb,Kc
c     Sabc = Ex         Ey         Ez
c              0          0          0
c
c******************************************************************************

c Define the number of shell components on each center.

      La2 = ((La+1)*(La+2))/2
      Lb2 = ((Lb+1)*(Lb+2))/2
      Lc2 = ((Lc+1)*(Lc+2))/2

c Loop over shell components.

      n = 0

      do 40 ica = 1,La2

       call getNxyz(La,ica,Nxyz)

       Ia = Nxyz(1)
       Ja = Nxyz(2)
       Ka = Nxyz(3)

       do 30 icb = 1,Lb2

        call getNxyz(Lb,icb,Nxyz)

        Ib = Nxyz(1)
        Jb = Nxyz(2)
        Kb = Nxyz(3)

        if( TriDiag )then
         icc_lim = icb
        else
         icc_lim = Lc2
        end if

        do 20 icc = 1,icc_lim

         call getNxyz(Lc,icc,Nxyz)

         Ic = Nxyz(1)
         Jc = Nxyz(2)
         Kc = Nxyz(3)

         n = n + 1

         Sabc(n) = 0.D0

         do 10 m = 1,NABC
          Sabc(n) = Sabc(n) + 
     &             E(m,1,0,Ia,Ib,Ic)*E(m,2,0,Ja,Jb,Jc)*E(m,3,0,Ka,Kb,Kc)

   10    continue

   20   continue

   30  continue

   40 continue

      end
