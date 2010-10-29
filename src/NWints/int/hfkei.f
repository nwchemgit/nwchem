      Subroutine hfkei(alpha,E,Tab,Ti,Nints,NPP,La,Lb,Li,canAB)
c $Id$

      Implicit none

      integer La,Lb,Li,Nints,NPP
      logical canAB

c--> Hermite Linear Expansion Coefficients

      double precision E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> Exponents

      double precision alpha(2,NPP)

c--> Kinetic Energy Integrals

      double precision Tab(Nints)

c--> Scratch Space

      integer Nxyz(3)
      double precision Ti(NPP)

c--> Local variables

      integer nn,ma,mb,mb_limit,m,La2,Lb2
      integer Ia,Ja,Ka, Ib,Jb,Kb
      double precision dia,dja,dka,dib,djb,dkb
c
c Compute the kinetic energy integrals.
c
c     Formula:                                                       
c
c            1  /  Ia,Ib   Ja,Jb   Ka,Kb                       \
c     Tab =  - | T      Ey      Ez      + "Y-term" + "Z-term"   |
c            2  \  X       0       0                           /
c                                                                     
c      i,j         i-1,j-1      i-1,j+1       i+1,j-1       i+1,j+1   
c     T      = ijEx       - 2ibEx      - 2ajEx       + 4abEx          
c      X           0             0            0             0
c                                                                     
c******************************************************************************
  
c Initialize the block of KEIs.


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
        dia = dble(ia)
        dja = dble(ja)
        dka = dble(ka)

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
          dib = dble(ib)
          djb = dble(jb)
          dkb = dble(kb)

          nn = nn + 1
          tab(nn)=0d0
  
c Build Tx.

          if( Ia.gt.0 .and. Ib.gt.0 )then
            do m = 1,NPP
              Ti(m) =   0.5D0*(   dia*dib       )*E(1,m,0,Ia-1,Ib-1)
     &            -       (       dia*alpha(2,m))*E(1,m,0,Ia-1,Ib+1)
     &            -       (alpha(1,m)*dib       )*E(1,m,0,Ia+1,Ib-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
            end do
          else if( Ia.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (       dIa*alpha(2,m))*E(1,m,0,Ia-1,Ib+1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
            end do
          else if( Ib.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (alpha(1,m)*dIb       )*E(1,m,0,Ia+1,Ib-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
            end do
          else
            do m = 1,NPP
              Ti(m) = 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
            end do
          end if
  
c Add Tx*Ey*Ez to Tab
  
          do m = 1,NPP
            Tab(nn) = Tab(nn) + Ti(m)*E(2,m,0,Ja,Jb)*E(3,m,0,Ka,Kb)
          end do

c Build Ty.
  
          if( Ja.gt.0 .and. Jb.gt.0 )then
            do m = 1,NPP
              Ti(m) = 0.5D0*(     dJa*dJb       )*E(2,m,0,Ja-1,Jb-1)
     &            -       (       dJa*alpha(2,m))*E(2,m,0,Ja-1,Jb+1)
     &            -       (alpha(1,m)*dJb       )*E(2,m,0,Ja+1,Jb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
            end do
          else if( Ja.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (       dJa*alpha(2,m))*E(2,m,0,Ja-1,Jb+1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
            end do
          else if( Jb.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (alpha(1,m)*dJb       )*E(2,m,0,Ja+1,Jb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
            end do
          else
            do m = 1,NPP
              Ti(m) = 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
            end do
          end if
  
c Add Ex*Ty*Ez to Tab.
  
          do m = 1,NPP
            Tab(nn) = Tab(nn) + E(1,m,0,Ia,Ib)*Ti(m)*E(3,m,0,Ka,Kb)
          end do
  
c Build Tz.
  
          if( Ka.gt.0 .and. Kb.gt.0 )then
            do m = 1,NPP
              Ti(m) = 0.5D0*(     dKa*dKb       )*E(3,m,0,Ka-1,Kb-1)
     &            -       (       dKa*alpha(2,m))*E(3,m,0,Ka-1,Kb+1)
     &            -       (alpha(1,m)*dKb       )*E(3,m,0,Ka+1,Kb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
            end do
          else if( Ka.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (       dKa*alpha(2,m))*E(3,m,0,Ka-1,Kb+1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
            end do
          else if( Kb.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (alpha(1,m)*dKb       )*E(3,m,0,Ka+1,Kb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
            end do
          else
            do m = 1,NPP
              Ti(m) = 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
            end do
          end if
  
c Add Ex*Ey*Tz to Tab.
  
          do m = 1,NPP
            Tab(nn) = Tab(nn) + E(1,m,0,Ia,Ib)*E(2,m,0,Ja,Jb)*Ti(m)
          end do

        end do

      end do
  
      end
************************************************************************
      Subroutine hfkei_gc(alpha,E,Tab,TabP,TabH,Ti,Acoefs,Bcoefs,ipairp,
     &    NPA,NPB,NCA,NCB,NPP,La,Lb,La2,Lb2,Li,canAB)
c
      implicit none

      integer NPA,NPB,NCA,NCB,NPP,La,Lb,La2,Lb2,Li
      logical canAB

c--> Hermite Linear Expansion Coefficients

      double precision E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> Index of primitives

      integer ipairp(2,NPP)

c--> Exponents

      double precision alpha(2,NPP)

c--> Kinetic Energy Integrals

      double precision Tab(Lb2,ncb,La2,nca)
      double precision TabP(NPP)
      double precision TabH(NPA,NCB)

c--> general contraction matrices

      double precision Acoefs(NPA,NCA)
      double precision Bcoefs(NPB,NCB)

c--> Scratch Space

      double precision Ti(NPP)
      integer Nxyz(3)

c--> Local variables

      integer ma,mb,m, ica,icb,icb_limit,ipa,ipb
      integer Ia,Ja,Ka, Ib,Jb,Kb
      double precision dia,dja,dka,dib,djb,dkb
c
c Compute the kinetic energy integrals.
c
c     Formula:                                                       
c
c            1  /  Ia,Ib   Ja,Jb   Ka,Kb                       \
c     Tab =  - | T      Ey      Ez      + "Y-term" + "Z-term"   |
c            2  \  X       0       0                           /
c                                                                     
c      i,j         i-1,j-1      i-1,j+1       i+1,j-1       i+1,j+1   
c     T      = ijEx       - 2ibEx      - 2ajEx       + 4abEx          
c      X           0             0            0             0
c                                                                     
c******************************************************************************
  
c Initialize the block of KEIs.

      call dfill(La2*Lb2*nca*ncb,0.0d00,Tab,1)

c Loop over shell components.

      do ma = 1,La2

c Define the angular momentum indices for shell "A".

        call getNxyz(La,ma,Nxyz)

        Ia = Nxyz(1)
        Ja = Nxyz(2)
        Ka = Nxyz(3)
        dia = dble(ia)
        dja = dble(ja)
        dka = dble(ka)

        do mb = 1,Lb2

c Define the angular momentum indices for shell "B".

          call getNxyz(Lb,mb,Nxyz)

          Ib = Nxyz(1)
          Jb = Nxyz(2)
          Kb = Nxyz(3)
          dib = dble(ib)
          djb = dble(jb)
          dkb = dble(kb)

          call dfill(NPP,0.0d00,TabP,1)
  
c Build Tx.
  
          if( Ia.gt.0 .and. Ib.gt.0 )then
            do m = 1,NPP
              Ti(m) =   0.5D0*(   dia*dib       )*E(1,m,0,Ia-1,Ib-1)
     &            -       (       dia*alpha(2,m))*E(1,m,0,Ia-1,Ib+1)
     &            -       (alpha(1,m)*dib       )*E(1,m,0,Ia+1,Ib-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
            end do
          else if( Ia.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (       dIa*alpha(2,m))*E(1,m,0,Ia-1,Ib+1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
            end do
          else if( Ib.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (alpha(1,m)*dIb       )*E(1,m,0,Ia+1,Ib-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
            end do
          else
            do m = 1,NPP
              Ti(m) = 2.0D0*(alpha(1,m)*alpha(2,m))*E(1,m,0,Ia+1,Ib+1)
            end do
          end if
  
c Add Tx*Ey*Ez to Tab
  
          do m = 1,NPP
            TabP(m) = TabP(m) + Ti(m)*E(2,m,0,Ja,Jb)*E(3,m,0,Ka,Kb)
          end do
  
c Build Ty.
  
          if( Ja.gt.0 .and. Jb.gt.0 )then
            do m = 1,NPP
              Ti(m) = 0.5D0*(     dJa*dJb       )*E(2,m,0,Ja-1,Jb-1)
     &            -       (       dJa*alpha(2,m))*E(2,m,0,Ja-1,Jb+1)
     &            -       (alpha(1,m)*dJb       )*E(2,m,0,Ja+1,Jb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
            end do
          else if( Ja.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (       dJa*alpha(2,m))*E(2,m,0,Ja-1,Jb+1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
            end do
          else if( Jb.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (alpha(1,m)*dJb       )*E(2,m,0,Ja+1,Jb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
            end do
          else
            do m = 1,NPP
              Ti(m) = 2.0D0*(alpha(1,m)*alpha(2,m))*E(2,m,0,Ja+1,Jb+1)
            end do
          end if
  
c Add Ex*Ty*Ez to Tab.
  
          do m = 1,NPP
            TabP(m) = TabP(m) + E(1,m,0,Ia,Ib)*Ti(m)*E(3,m,0,Ka,Kb)
          end do
  
c Build Tz.
  
          if( Ka.gt.0 .and. Kb.gt.0 )then
            do m = 1,NPP
              Ti(m) = 0.5D0*(     dKa*dKb       )*E(3,m,0,Ka-1,Kb-1)
     &            -       (       dKa*alpha(2,m))*E(3,m,0,Ka-1,Kb+1)
     &            -       (alpha(1,m)*dKb       )*E(3,m,0,Ka+1,Kb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
            end do
          else if( Ka.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (       dKa*alpha(2,m))*E(3,m,0,Ka-1,Kb+1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
            end do
          else if( Kb.gt.0 )then
            do m = 1,NPP
              Ti(m) = -   (alpha(1,m)*dKb       )*E(3,m,0,Ka+1,Kb-1)
     &            + 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
            end do
          else
            do m = 1,NPP
              Ti(m) = 2.0D0*(alpha(1,m)*alpha(2,m))*E(3,m,0,Ka+1,Kb+1)
            end do
          end if
  
c Add Ex*Ey*Tz to Tab.
  
          do m = 1,NPP
            TabP(m) = TabP(m) + E(1,m,0,Ia,Ib)*E(2,m,0,Ja,Jb)*Ti(m)
          end do

c Contract over B shell

          call dfill(NCB*NPA,0.0d00,TabH,1)
          do icb = 1,NCB
            do m = 1,NPP
              ipa = ipairp(1,m)
              ipb = ipairp(2,m)
              TabH(ipa,icb) = TabH(ipa,icb)+TabP(m)*Bcoefs(ipb,icb)
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
                Tab(mb,icb,ma,ica) = Tab(mb,icb,ma,ica) +
     &              TabH(ipa,icb)*Acoefs(ipa,ica)
              end do
            end do
          end do

        end do

      end do

      if (canAB) call canon_ab(Tab,Tab,Lb2*NCB,La2*NCA)

      end
