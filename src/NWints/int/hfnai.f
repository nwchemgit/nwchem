      Subroutine hfnai(E,R0,IJK,Vab,Nints,NPP,La,Lb,Li,Lp,Lp3,canAB)
c $Id$
      
      Implicit none
      
      integer Nints,NPP,La,Lb,Li,Lp,Lp3
      Logical canAB
      
c--> Hermite Linear Expansion Coefficients
      
      double precision E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))
      
c--> Auxiliary Function Integrals & Index
      
      double precision R0(NPP,Lp3)
      integer IJK(0:Lp,0:Lp,0:Lp)
      
c--> Nuclear Attraction Integrals
      
      double precision Vab(Nints)
      
c--> Scratch Space
      
      integer Nxyz(3)

c--> Local variables

      integer nn,ma,mb,mb_limit,mp,np,La2,Lb2
      integer Ia,Ja,Ka, Ib,Jb,Kb, Ip,Jp,Kp
      double precision vab_int
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

c Initialize the block of NAIs.
      
!      call dcopy(Nints,0d0,0,Vab,1)
      
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
          
          nn = nn + 1
          
          vab_int = 0.0d00
          
          do Ip = 0,Ia+Ib
            do Jp = 0,Ja+Jb
              do Kp = 0,Ka+Kb
                
                np = IJK(Ip,Jp,Kp)
            
                do mp = 1,NPP
                  vab_int = vab_int +
     &                E(1,mp,Ip,Ia,Ib)*
     &                E(2,mp,Jp,Ja,Jb)*
     &                E(3,mp,Kp,Ka,Kb)*R0(mp,np)
              
                end do

              end do
            end do
          end do

          Vab(nn) = vab_int

        end do

      end do

      end
*----------------------------------------------------------------------
      Subroutine hfnai_gc(E,R0,IJK,Vab,VabP,VabH,Acoefs,Bcoefs,ipairp,
     &    NPA,NPB,NCA,NCB,NPP,La,Lb,La2,Lb2,Li,Lp,Lp3,canAB)

      implicit none

      logical canAB
      integer NPA,NPB,NCA,NCB,NPP,La,Lb,La2,Lb2,Li,Lp,Lp3

c--> Hermite Linear Expansion Coefficients

      double precision E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> Index of primitives

      integer ipairp(2,NPP)

c--> Auxiliary Function Integrals & Index

      double precision R0(NPP,Lp3)
      integer IJK(0:Lp,0:Lp,0:Lp)

c--> Nuclear Attraction Integrals

      double precision Vab(Lb2,NCB,La2,NCA)
      double precision VabP(NPP)
      double precision VabH(NPA,NCB)

c--> general contraction matrices

      double precision Acoefs(NPA,NCA)
      double precision Bcoefs(NPB,NCB)

c--> Scratch Space

      integer Nxyz(3)

c--> Local variables

      integer ma,mb,mp,np, ica,icb,icb_limit,ipa,ipb
      integer Ia,Ja,Ka, Ib,Jb,Kb, Ip,Jp,Kp
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

c Initialize the block of NAIs.

      call dcopy(La2*Lb2*nca*ncb,0d0,0,Vab,1)

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
          call dcopy(NPP,0d0,0,VabP,1)
          do Ip = 0,Ia+Ib
            do Jp = 0,Ja+Jb
              do Kp = 0,Ka+Kb

                np = IJK(Ip,Jp,Kp)

                do mp = 1,NPP
                  VabP(mp) = VabP(mp) + R0(mp,np)*
     &                E(1,mp,Ip,Ia,Ib)*
     &                E(2,mp,Jp,Ja,Jb)*
     &                E(3,mp,Kp,Ka,Kb)
                end do

              end do
            end do
          end do

c Contract over B shell

          call dcopy(NCB*NPA,0d0,0,VabH,1)
          do icb = 1,NCB
            do mp = 1,NPP
              ipa = ipairp(1,mp)
              ipb = ipairp(2,mp)
              VabH(ipa,icb) = VabH(ipa,icb)+VabP(mp)*Bcoefs(ipb,icb)
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
                Vab(mb,icb,ma,ica) = Vab(mb,icb,ma,ica) +
     &              VabH(ipa,icb)*Acoefs(ipa,ica)
              end do
            end do
          end do

        end do

      end do

      if (canAB) call canon_ab(Vab,Vab,Lb2*NCB,La2*NCA)
      
      end
************************************************************************
      subroutine canon_ab (Vcanon,Vsquare,NB,NA)
      implicit none
      integer NA,NB,i,j,k
      double precision Vcanon(*),Vsquare(NB,NA)

      k = 1
      do i = 1,NA
        do j = 1,i
          k = k+1
          Vcanon(k) = Vsquare(j,i)
        end do
      end do

      end
