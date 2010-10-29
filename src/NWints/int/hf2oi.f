      Subroutine hf2oi(E,Sab,Nints,NPP,La,Lb,Li,canAB)
c $Id$

      Implicit none

      integer Nints,NPP,La,Lb,Li
      logical canAB

c--> Hermite Linear Expansion Coefficients

      double precision E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> 2-Center Overlap Integrals

      double precision Sab(Nints)

      integer Nxyz(3)

c--> Local variables

      integer nn,ma,mb,mb_limit,mp,La2,Lb2
      integer Ia,Ja,Ka, Ib,Jb,Kb
      double precision sabval
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

c Initialize the block of integrals


c Define the number of shell components on each center.

      La2 = ((La+1)*(La+2))/2
      Lb2 = ((Lb+1)*(Lb+2))/2

c Loop over shell components.

      nn = 0

      do ma = 1,La2

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

        do mb = 1,mb_limit

c Define the angular momentum indices for shell "B".

          call getNxyz(Lb,mb,Nxyz)

          Ib = Nxyz(1)
          Jb = Nxyz(2)
          Kb = Nxyz(3)


          sabval=0d0
          do mp = 1,NPP
            sabval = sabval + E(1,mp,0,Ia,Ib)*
     &                          E(2,mp,0,Ja,Jb)*
     &                          E(3,mp,0,Ka,Kb)
          end do
          nn = nn + 1
          Sab(nn) = sabval

        end do

      end do

      end
************************************************************************
      Subroutine hf2oi_gc(E,Sab,SabP,SabH,Acoefs,Bcoefs,ipairp,
     &    NPA,NPB,NCA,NCB,NPP,La,Lb,La2,Lb2,Li,canAB)
c
      implicit none

      integer NPA,NPB,NCA,NCB,NPP,La,Lb,La2,Lb2,Li
      logical canAB

c--> Hermite Linear Expansion Coefficients

      double precision E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> Index of primitives

      integer ipairp(2,NPP)

c--> Kinetic Energy Integrals

      double precision Sab(Lb2,ncb,La2,nca)
      double precision SabP(NPP)
      double precision SabH(NPA,NCB)

c--> general contraction matrices

      double precision Acoefs(NPA,NCA)
      double precision Bcoefs(NPB,NCB)

c--> Scratch Space

      integer Nxyz(3)

c--> Local variables

      integer ma,mb,mp, ica,icb,icb_limit,ipa,ipb
      integer Ia,Ja,Ka, Ib,Jb,Kb
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

c Initialize the block of integrals

      call dfill(La2*Lb2*nca*ncb,0.0d00,Sab,1)

c Loop over shell components.

      do ma = 1,La2

c Define the angular momentum indices for shell "A".

        call getNxyz(La,ma,Nxyz)

        Ia = Nxyz(1)
        Ja = Nxyz(2)
        Ka = Nxyz(3)

        do mb = 1,Lb2

c Define the angular momentum indices for shell "B".

          call getNxyz(Lb,mb,Nxyz)

          Ib = Nxyz(1)
          Jb = Nxyz(2)
          Kb = Nxyz(3)

          do mp = 1,NPP
            SabP(mp) = E(1,mp,0,Ia,Ib)*E(2,mp,0,Ja,Jb)*E(3,mp,0,Ka,Kb)
          end do

c Contract over B shell

          call dfill(NCB*NPA,0.0d00,SabH,1)
          do icb = 1,NCB
            do mp = 1,NPP
              ipa = ipairp(1,mp)
              ipb = ipairp(2,mp)
              SabH(ipa,icb) = SabH(ipa,icb)+SabP(mp)*Bcoefs(ipb,icb)
            end do
          end do

c Contract over A shell
        
          do ica = 1,NCA
            if( canAB )then
              icb_limit = ica
            else
              icb_limit = NCB
            end if
            do icb = 1,icb_limit
              do ipa = 1,NPA
                Sab(mb,icb,ma,ica) = Sab(mb,icb,ma,ica) +
     &              SabH(ipa,icb)*Acoefs(ipa,ica)
              end do
            end do
          end do

        end do

      end do

      if (canAB) call canon_ab(Sab,Sab,Lb2*NCB,La2*NCA)

      end
