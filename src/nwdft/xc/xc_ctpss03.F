
c    Tao,Perdew, Staroverov, Scuseria exchange functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density
c                              tauW - von Weiszacker kinetic energy density
c                              tauU - uniform-gas KE density


      Subroutine xc_ctpss03(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght,
     &                     tau, Amat, Cmat, Mmat)     
     

c
c$Id: xc_ctpss03.F,v 1.6 2006-01-25 01:03:33 edo Exp $
c

c     References:
c     [a] J. Tao, J.P. Perdew, V.N.Staroverov, G. Scuseria
c         PRL 91, 146401 (2003).

      Implicit none
c
c
c     Input and other parameters
c
      integer ipol, nq
      double precision dummy(1)

      double precision cfac
      logical lcfac, nlcfac
c      logical  ldew
c      double precision ffunc(*)

      logical lfac, nlfac
      double precision fac
      double precision tol_rho

c
c     Threshold parameters
c
      double precision TOLL, EXPTOL
      double precision EPS
      parameter (TOLL = 1.0D-40, EXPTOL = 40.0d0)
      parameter (EPS = 1.0e-8)
c
c     Correlation energy
c
      double precision Ec
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
 
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol), gammaval, gam12
      
c
c     Kinetic Energy Density
c
      double precision tau(nq,3,ipol), tauN
 
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision Mmat(nq,*)

      integer n
      double precision rhoval


c    call to the cPBE subroutine

      double precision  neGGA, dneGGAdn(2), dneGGAdg(3)
      double precision  facG, Gqwght, epGGA
      logical  lfacG, nlfacG
      integer Gnq
      double precision rho_t(3), delrho_t(3,2)
      double precision  gammavalA, gammavalB
      double precision  tauNA,tauNB

c     spectral representation terms
c      double precision aterma1, aterma2, aterma3
c      double precision atermb1, atermb2, atermb3
c      double precision aterm1, aterm2, aterm3
c      double precision cterm1, cterm2, cterm3
c      double precision mterm
c      double precision ctermb1, ctermb2, ctermb3

 

c  FULLY SPIN POLARIZED SECTION	
      double precision   rhovalF, epFSP
      double precision  neFSP, dneFSPdn(2), dneFSPdg(3)

      integer  ipolF, Fnq
      double precision  facF, Fqwght
      logical  lfacF, nlfacF
      double precision delrho_A(3,2), rho_A(3)



c     TPSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs

      double precision THRD, FACE, F73, F43
      double precision zeta, ccc, fer
      double precision tauwup, tauwdn, tauw
      double precision xx2,xx3, ter1,ter2,ter3
      double precision sumz, grz, szeta2, term1, term2
      double precision yyz, pbe, en
      double precision aa,bb,cc,dd,dnu
      double precision up,dn,grup, grdn, gr, grupdn
      double precision dens
      double precision pbeup, pbedn, delc, revsic1, revsic2

c     derivsssssssssssssssssssssssss

      double precision  decggadn, dtwt2dn, dtwt3dn
      double precision  dpbeupdn, drevdn
      double precision  atermn, btermn
      double precision  cpartn(3), apartn, bpartn, finaln
      double precision  drevdg, dpbeupdg
      double precision  dtwt2dg, dtwt3dg, decggadg
      double precision  atermg, btermg
      double precision  apartg, cpartg(2), finalg
      double precision  drevdgab, finalgab
      double precision drevdt, apartt, bpartt(2),finalt

      double precision dzetadna, dzetadnb
      double precision dc1dna, dc1dnb
      double precision dter1dna, dter2dna, dter3dna
      double precision dter1dnb, dter2dnb, dter3dnb
      double precision dlongterdna, dlongterdnb
      double precision dgrzdna, dgrzdnb
      double precision dxidna, dxidnb
      double precision dbaddna, dbaddnb
      double precision ddennydna, ddennydnb
      double precision dCCCdna, dCCCdnb, CCC1
      double precision dtwt2dna, dtwt2dnb
      double precision decggadna, decggadnb, delcnt
      double precision deFSPdna, dsmdn, dsmdna, dsmdnb
      double precision quima(4), quimb(4), drevdna, drevdnb
      double precision dtwt3dna, dtwt3dnb, psea(2), pseb(2)

      double precision dter1dgaa,dgrzdgaa, dxidgaa, dbaddgaa
      double precision ddennydgaa, dCCCdgaa, dtwt2dgaa
      double precision deggadgaa, dsmdgaa, quimga(4),drevdgaa
      double precision dtwt3dgaa, psega(2)
      double precision dter1dgbb,dgrzdgbb, dxidgbb
      double precision dbaddgbb,ddennydgbb,dCCCdgbb,dtwt2dgbb
      double precision deggadgbb,dsmdgbb,quimgb(4),drevdgbb
      double precision dtwt3dgbb,psegb(2)
      double precision dter3dgab,dgrzdgab,dxidgab,dbaddgab
      double precision ddennydgab, dCCCdgab,dtwt2dgab
      double precision deggadgab,dsmdgab,quimgab(4)
      double precision dtwt3dgab,psegab(2)
      double precision dtwt2dt, dtwt3dt, pset(2)

      parameter ( THRD = 0.33333333333333333333333333333d0 )
      parameter ( FACE = 3.09366772628013593097d0 )
c      FACE = (3*PI**2)**(1/3)
      parameter ( F73 = 2.3333333333333333333333d0)
      parameter ( F43 = 1.3333333333333333333333d0)
       

      aa = 0.87d0
      bb = 0.50d0
      cc = 2.26d0
c      dd = 0.d0   !set dd=0 to recover pkzb spin restricted cases   
      dd = 2.8d0

      dnu = 1.d0
c     TPSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs

      fac = cfac
      lfac = lcfac
      nlfac = nlcfac


      if (ipol.eq.1 )then
c        ======> SPIN-RESTRICTED <======

         do 22  n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 22

         rhoval = rho(n,1)

C   set up values to call PBE subroutine
         rho_t(1) = rho(n,1)
c do delrho
            delrho_t(1,1) = delrho(n,1,1)
            delrho_t(2,1) = delrho(n,2,1)
            delrho_t(3,1) = delrho(n,3,1)
     
         neGGA = 0.0d0  !Ec in PBE
         dneGGAdn(1) = 0.0d0   !Amat in PBE

         dneGGAdg(1) = 0.0d0  !Cmat in PBE
         dneGGAdg(2) = 0.0d0  !Cmat in PBE

         Gnq = 1
         Gqwght = 1.0d0

            facG = 1.0d0 !MUST BE ONE
            lfacG = .false.  !MUST BE F
            nlfacG = .true. !MUST BE a T/F 

            call xc_cMpbe96(tol_rho, facG, lfacG, nlfacG,
     &         rho_t, delrho_t,
     &         dneGGAdn,dneGGAdg, 
     &          Gnq, ipol, neGGA,Gqwght,
     &         .false.,dummy)


c
c        epGGA = n*(epsilon_c^GGA) / n =cor. energy per electron
c        epGGA= ec^LDA +H  = pbe


            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1)
            gam12=dsqrt(gammaval)

            tauN = tau(n,1,1) + 
     &             tau(n,2,1) + 
     &             tau(n,3,1)

        up=rhoval*0.5d0
        dn=rhoval*0.5d0
        grup=gam12*0.5d0
        grdn=gam12*0.5d0
        gr=gam12
        grupdn=gammaval*0.25d0
        pbe = neGGA/rhoval

        dens = rhoval

        zeta=0.d0

         CCC = 0.53d0 !since zeta=0
         fer = FACE*dens**THRD


         tauw = 0.125d0*gammaval/dens

c         tau = tauup + taudn = tauN

         xx2 = (tauw/tauN)**2.d0
         xx3 = (tauw/tauN)**3.d0

          en = pbe*(1.d0 + CCC*xx2)

c    set up values to call PBE subroutine as FSP system
c   n-up = n-dn = n/2   hence we can just call the subroutine once


         rhovalF  = (0.5d0)*rho(n,1)  ! set total density as half
         rho_A(1) = (0.5d0)*rho(n,1)   ! total   equals (1/2)n_tot
         rho_A(2) = (0.5d0)*rho(n,1)   ! alpha   equals (1/2)n_tot
         rho_A(3) = 0.d0               ! beta  equals zero
c do delrho
            delrho_A(1,1) = (0.5d0)*delrho(n,1,1)  ! nabla n_up x  
            delrho_A(2,1) = (0.5d0)*delrho(n,2,1)  ! nabla n_up y
            delrho_A(3,1) = (0.5d0)*delrho(n,3,1)  ! nabla n_up z
     
        
            delrho_A(1,2) = 0.d0    ! set beta gradient to zero
            delrho_A(2,2) = 0.d0    ! set beta gradient to zero
            delrho_A(3,2) = 0.d0    ! set beta gradient to zero


         neFSP = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE

         Fnq = 1
         Fqwght = 1.0d0
          ipolF=2

            facF = 1.0d0 
            lfacF = .false.  
            nlfacF = .true. 

            call xc_cMpbe96(tol_rho, facF, lfacF, nlfacF,
     &         rho_A, delrho_A,
     &         dneFSPdn,dneFSPdg, 
     &          Fnq, ipolF, neFSP, Fqwght,
     &         .false.,dummy)


c         if(abs(neFSP).lt.tol_rho) goto 22
c
c        epsilon_c = n*(epsilon_c^GGA) / n
c
         pbeup = neFSP/rhovalF

c        functional deriv info below fffffffffffff

          dtwt2dn = -2.d0*xx2/rhoval
          dtwt3dn = -3.d0*xx3/rhoval
          decggadn= dneGGAdn(1)/rhoval - neGGA/rhoval**2 


          dtwt2dg = 2.d0*0.125d0*tauw/(rhoval*tauN**2)


          dtwt3dg = 3.d0*xx2*0.125d0/(rhoval*tauN)
          decggadg= dneGGAdg(1)/dens

       if (pbeup.lt.pbe) then
              delc= xx2*pbe
C  eps-tilda is eps_c
C  functional deriv info below fffffffffffffffff

          drevdn= -pbe*dtwt2dn+decggadn*(1.d0 - xx2)
          drevdg= -pbe*dtwt2dg+decggadg*(1.d0 -xx2)
          drevdgab=(1.d0-xx2)*dneGGAdg(2)/dens

          drevdt= 2.d0*pbe*xx2/tauN
     
       else
              delc= xx2*pbeup

C  eps-tilda is eps^FSP
C  functional deriv info below fffffffffffffffff

corig        dpbeupdn = 0.5d0*dneFSPdn(1)/rhovalF- 2.d0*neFSP/rhovalF**2
        dpbeupdn = 0.5d0*dneFSPdn(1)/rhovalF- 2.d0*neFSP/rhoval**2
c  above note the .5's.  you are taking deriv wrt total density n
c                        not deriv wrt n_up
        dpbeupdg = 0.25d0*dneFSPdg(1)/rhovalF
c  note .25 above is because you want gamma=deln_tot*deln_tot


           atermn=pbe*CCC*dtwt2dn+(1.d0+CCC*xx2)*decggadn
           btermn=(1.d0+CCC)*xx2*dpbeupdn +
     &             pbeup*(1.d0+CCC)*dtwt2dn
          drevdn=atermn - btermn

          atermg=pbe*ccc*dtwt2dg+(1.d0+CCC*xx2)*decggadg
          btermg=(1.d0+CCC)*(xx2*dpbeupdg+pbeup*dtwt2dg)
          drevdg=atermg-btermg
          drevdgab=(1.d0+CCC*xx2)*dneGGAdg(2)/dens

         drevdt=(CCC*pbe-(1.d0+CCC)*pbeup)*xx2*(-2.d0/tauN)


       endif

         delc = -(1.d0 + CCC)*delc
         revsic1 = en + delc
         revsic2 = revsic1*(1.d0 + dd*revsic1*xx3)

             Ec = Ec + rhoval*revsic2*qwght(n)*fac

c  functional derivative eqns nnnnnnnnnnnnnnnn
c          derivs wrt n
         cpartn(1)=dens*revsic1*revsic1*dtwt3dn
         cpartn(2)=dens*xx3*2.d0*revsic1*drevdn
         cpartn(3)=revsic1*revsic1*xx3
         apartn=dens*drevdn
         bpartn=revsic1
         finaln=apartn+bpartn+
     &    dd*(cpartn(1)+cpartn(2)+cpartn(3) )

          Amat(n,1)=Amat(n,1)+(finaln)*fac


c  functional derivative eqns gggggggggggggg
c          derivs wrt g
         apartg=dens*drevdg
         cpartg(1)=revsic1*revsic1*dtwt3dg
         cpartg(2)=xx3*2.d0*revsic1*drevdg
         finalg=apartg+dens*dd*( cpartg(1)+cpartg(2) )

         Cmat(n,1)=Cmat(n,1)+ finalg*fac

         finalgab=dens*drevdgab+dd*dens*xx3*2.d0*revsic1*drevdgab
         Cmat(n,2)=Cmat(n,2)+finalgab*fac

c  functional derivative eqns tttttttttttttttt
c          derivs wrt t
         apartt=dens*drevdt
         bpartt(1)=-revsic1*revsic1*xx3*3.d0/tauN
         bpartt(2)=xx3*2.d0*revsic1*drevdt


         finalt=apartt+dd*dens*(bpartt(1)+bpartt(2))



        Mmat(n,1)=Mmat(n,1)+0.5d0*finalt*fac


c   analytically, from the equation (17) in the tpss paper
c    sum(n-sigma/n-tot)* eps^tilda
c    = n-up*eps^tilda-up/n  +  n-dn*eps^tilda-dn/n
c    but eps^tilda-up = eps^tilda-dn and n-up=n-dn
c    so the sum part of (17) reduces to eps^tilda



22     continue




      else !ipol=2



c        ======> SPIN-UNRESTRICTED <======
c                     or
c                SPIN-POLARIZED
c
c
c         ***** SPIN POLARIZATION  *****

         do 20 n = 1, nq
c            

             if (rho(n,1).lt.tol_rho) goto 20

         rhoval = rho(n,1)
         rho_t(1) = rho(n,1)
         rho_t(2) = rho(n,2)
         rho_t(3) = rho(n,3)
         
c do delrho
            delrho_t(1,1) = delrho(n,1,1)
            delrho_t(2,1) = delrho(n,2,1)
            delrho_t(3,1) = delrho(n,3,1)
     
        
            delrho_t(1,2) = delrho(n,1,2)
            delrho_t(2,2) = delrho(n,2,2)
            delrho_t(3,2) = delrho(n,3,2)


         neGGA = 0.0d0  !Ec in PBE
         dneGGAdn(1) = 0.0d0   !Amat in PBE (n,1)
         dneGGAdn(2) = 0.0d0   !Amat in PBE (n,2)


         dneGGAdg(1) = 0.0d0  !Cmat in PBE--aa
         dneGGAdg(2) = 0.0d0  !Cmat in PBE--ab
         dneGGAdg(3) = 0.0d0  !Cmat in PBE--bb


         Gnq = 1
         Gqwght = 1.0d0

            facG = 1.0d0 
            lfacG = .false.  
            nlfacG = .true. 

            call xc_cMpbe96(tol_rho, facG, lfacG, nlfacG,
     &         rho_t, delrho_t,
     &         dneGGAdn,dneGGAdg, 
     &          Gnq, ipol, neGGA,Gqwght,
     &         .false.,dummy)



            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,1,2)*delrho(n,1,2) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,2,2)*delrho(n,2,2) +
     &                 delrho(n,3,1)*delrho(n,3,1) +
     &                 delrho(n,3,2)*delrho(n,3,2) +
     &           2.d0*(delrho(n,1,1)*delrho(n,1,2) +
     &                 delrho(n,2,1)*delrho(n,2,2) +
     &                 delrho(n,3,1)*delrho(n,3,2))

              gam12=dsqrt(gammaval)

            tauNA = tau(n,1,1) + 
     &              tau(n,2,1) + 
     &              tau(n,3,1)


            tauNB = tau(n,1,2) + 
     &              tau(n,2,2) + 
     &              tau(n,3,2)   

        tauN=tauNA + tauNB


        up=rho(n,2)
        dn=rho(n,3)

        grup= dsqrt(delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1) )

        grdn= dsqrt(  delrho(n,1,2)*delrho(n,1,2) +
     &                 delrho(n,2,2)*delrho(n,2,2) +
     &                 delrho(n,3,2)*delrho(n,3,2) )

        gr=gam12
        grupdn=    delrho(n,1,1)*delrho(n,1,2) +
     &             delrho(n,2,1)*delrho(n,2,2) +
     &             delrho(n,3,1)*delrho(n,3,2)
        pbe = neGGA/rhoval
        dens = rhoval
C jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj

         zeta = (up - dn)/dens

         fer = FACE*dens**THRD
         
c         if(dens.lt.tol_rho) then
c            tauw=0.d0
c         else
c            tauw = 0.125d0*gammaval/dens
c         endif

            tauw = 0.125d0*gammaval/dens
         

         xx2 = (tauw/tauN)**2.d0
         xx3 = (tauw/tauN)**3.d0

c----------------------------
         ter1 = (1.d0 - zeta)**2.d0*grup**2.d0
         ter2 = (1.d0 + zeta)**2.d0*grdn**2.d0
         ter3 = 2.d0*(1.d0 - zeta**2.d0)*grupdn
         sumz = ter1 + ter2 - ter3
         sumz = dabs(sumz)
         grz = dsqrt(sumz)/dens

         szeta2 = (grz/(2.d0*fer))**2.d0

         if (dabs(zeta).lt.0.9999d0) then
            term1 = 1.d0/(1.d0 + zeta)**(4.d0/3.d0)
            term2 = 1.d0/(1.d0 - zeta)**(4.d0/3.d0) 
            yyz = szeta2*(term1 + term2)/2.d0
         else
            yyz = 0.d0
         end if

c           do we need below, or is it taken care of above
          if(zeta.lt.-1.0d0) zeta=-1.0d0
          if(zeta.gt. 1.0d0) zeta= 1.0d0


         CCC = 0.53d0 + aa*zeta**2.d0 + bb*zeta**4.d0 + cc*zeta**6.d0
         CCC1 =0.53d0 + aa*zeta**2.d0 + bb*zeta**4.d0 + cc*zeta**6.d0

         CCC = CCC/(1.d0 + dnu*yyz)**4.d0

          en = pbe*(1.d0 + CCC*xx2)
         
         dzetadna= 2.d0*rho(n,3)/(rhoval**2)
         dzetadnb=-2.d0*rho(n,2)/(rhoval**2)

         dc1dna=(1.74d0*zeta+2.d0*zeta**3 
     &         +13.56d0*zeta**5)*dzetadna
         dc1dnb=(1.74d0*zeta+2.d0*zeta**3 
     &         +13.56d0*zeta**5)*dzetadnb


         dter1dna=-2.d0*(1.d0-zeta)*grup*grup*dzetadna
         dter2dna= 2.d0*(1.d0+zeta)*grdn*grdn*dzetadna
         dter3dna=-4.d0*zeta*grupdn*dzetadna

         dter1dnb=-2.d0*(1.d0-zeta)*grup*grup*dzetadnb
         dter2dnb= 2.d0*(1.d0+zeta)*grdn*grdn*dzetadnb
         dter3dnb=-4.d0*zeta*grupdn*dzetadnb
c  sumz = ter1 + ter2 - ter3
c  sumz = dabs(sumz)
c  grz = dsqrt(sumz)/dens

 

         if(sumz.lt.tol_rho) then
             dgrzdna=0.d0
         else
         dlongterdna=0.5d0*(dter1dna+dter2dna-dter3dna)/
     &               sumz**(.5d0)
         dgrzdna=(dlongterdna/dens)-
     &       (sumz**(.5d0))/dens**2
         end if
    

          if(sumz.lt.tol_rho) then
             dgrzdnb=0.d0
          else
         dlongterdnb=0.5d0*(dter1dnb+dter2dnb-dter3dnb)/
     &                (sumz**(.5d0))
         dgrzdnb=(dlongterdnb/dens)-
     &       (sumz**(.5d0))/dens**2
          end if
 

C         fer = FACE*dens**THRD
         dxidna=(-THRD*grz/( 2.d0*FACE*dens**(F43) ) ) +
     &        dgrzdna/(2.d0*fer)
         dxidnb=(-THRD*grz/( 2.d0*FACE*dens**(F43) ) ) +
     &        dgrzdnb/(2.d0*fer)
c let zterm=(1+z)^-4/3 + (1-z)^-4/3
c let  szeta2 = (  DELz/(2.d0*(3Pi*Pi dens)^1/3)  )**2.d0



c         szeta2 = (grz/(2.d0*fer))**2.d0

       if(dabs(1.d0-zeta).lt.tol_rho) then
       dbaddna= -F43*szeta2*dzetadna/(1.d0+zeta)**(F73) +
     &         2.d0*(grz/(2.d0*fer))*dxidna/(1.d0+zeta)**(F43)
  
       dbaddnb= -F43*szeta2*dzetadnb/(1.d0+zeta)**(F73) +
     &         2.d0*(grz/(2.d0*fer))*dxidnb/(1.d0+zeta)**(F43)


       else if (dabs(1.d0+zeta).lt.tol_rho) then
       dbaddna=  F43*szeta2*dzetadna/(1.d0-zeta)**(F73)    
     &  +      2.d0*(grz/(2.d0*fer))*dxidna/(1.d0-zeta)**(F43)

       dbaddnb=  F43*szeta2*dzetadnb/(1.d0-zeta)**(F73)    
     &  +      2.d0*(grz/(2.d0*fer))*dxidnb/(1.d0-zeta)**(F43)

        else
       dbaddna= -F43*szeta2*dzetadna/(1.d0+zeta)**(F73) +
     &         2.d0*(grz/(2.d0*fer))*dxidna/(1.d0+zeta)**(F43)
     &  +        F43*szeta2*dzetadna/(1.d0-zeta)**(F73)    
     &  +      2.d0*(grz/(2.d0*fer))*dxidna/(1.d0-zeta)**(F43)

       dbaddnb= -F43*szeta2*dzetadnb/(1.d0+zeta)**(F73) +
     &         2.d0*(grz/(2.d0*fer))*dxidnb/(1.d0+zeta)**(F43)
     &  +        F43*szeta2*dzetadnb/(1.d0-zeta)**(F73)    
     &  +      2.d0*(grz/(2.d0*fer))*dxidnb/(1.d0-zeta)**(F43)

        endif


c yyz = szeta2*(term1 + term2)/2.d0

       ddennydna=-2.d0*dbaddna/(1.d0+yyz)**(3)
       ddennydnb=-2.d0*dbaddnb/(1.d0+yyz)**(3)


       dCCCdna=CCC1*ddennydna+
     &         dc1dna*(1.d0+yyz)**(-4)
       dCCCdnb=CCC1*ddennydnb+
     &         dc1dnb*(1.d0+yyz)**(-4)



          dtwt2dna = -2.d0*xx2/rhoval
c          dtwt2dnb = -2.d0*xx2/rhoval
          dtwt2dnb = dtwt2dna

          dtwt3dna = -3.d0*xx3/rhoval
c          dtwt3dnb = -3.d0*xx3/rhoval
          dtwt3dnb = dtwt3dna

          decggadna= dneGGAdn(1)/rhoval - neGGA/rhoval**2 
          decggadnb= dneGGAdn(2)/rhoval - neGGA/rhoval**2 


c below is derivs  wrt gaa and gbb gab below ggggggggggggggggggggg
c         fer = FACE*dens**THRD

          dter1dgaa=(1.d0-zeta)**2

          if(sumz.lt.tol_rho) then
              dgrzdgaa=0.d0
          else
             dgrzdgaa =0.5d0*dter1dgaa/(dens*sumz**(0.5d0))
          end if
          dxidgaa  =dgrzdgaa/(2.d0*fer)

           
          dter1dgbb=(1.d0+zeta)**2
          if(sumz.lt.tol_rho) then
             dgrzdgbb=0.d0
          else
          dgrzdgbb =0.5d0*dter1dgbb/(dens*sumz**(0.5d0))
          end if
          dxidgbb  =dgrzdgbb/(2.d0*fer)
  

          dter3dgab=-2.d0*(1.d0-zeta**2)
          if(sumz.lt.tol_rho) then
             dgrzdgab=0.d0
          else
          dgrzdgab=0.5d0*dter3dgab/(dens*sumz**(0.5d0))
          end if
          dxidgab = dgrzdgab/(2.d0*fer)

        
       if(dabs(1.d0-zeta).lt.tol_rho) then
          dbaddgaa=(grz/(2.d0*fer))*dxidgaa/(1.d0+zeta)**(F43)
          dbaddgbb=(grz/(2.d0*fer))*dxidgbb/(1.d0+zeta)**(F43)
          dbaddgab=(grz/(2.d0*fer))*dxidgab/(1.d0+zeta)**(F43)

       else if(dabs(1.d0+zeta).lt.tol_rho) then
          dbaddgaa=(grz/(2.d0*fer))*dxidgaa/(1.d0-zeta)**(F43)
          dbaddgbb=(grz/(2.d0*fer))*dxidgbb/(1.d0-zeta)**(F43)
          dbaddgab=(grz/(2.d0*fer))*dxidgab/(1.d0-zeta)**(F43)


       else
          dbaddgaa=(grz/(2.d0*fer))*dxidgaa/(1.d0+zeta)**(F43)+
     &             (grz/(2.d0*fer))*dxidgaa/(1.d0-zeta)**(F43)
          dbaddgbb=(grz/(2.d0*fer))*dxidgbb/(1.d0+zeta)**(F43)+
     &             (grz/(2.d0*fer))*dxidgbb/(1.d0-zeta)**(F43)
          dbaddgab=(grz/(2.d0*fer))*dxidgab/(1.d0+zeta)**(F43)+
     &             (grz/(2.d0*fer))*dxidgab/(1.d0-zeta)**(F43)



       endif

c yyz = szeta2*(term1 + term2)/2.d0

          ddennydgaa=-4.d0*dbaddgaa/(1.d0+yyz)**(3)
          dCCCdgaa = CCC1*ddennydgaa
          dtwt2dgaa = 2.d0*(tauw/tauN)/(8.d0*dens*tauN)
          deggadgaa = dneGGAdg(1)/dens

          ddennydgbb=-4.d0*dbaddgbb/(1.d0+yyz)**(3)
          dCCCdgbb = CCC1*ddennydgbb
          dtwt2dgbb = dtwt2dgaa
          deggadgbb = dneGGAdg(3)/dens

          ddennydgab=-4.d0*dbaddgab/(1.d0+yyz)**(3)
          dCCCdgab = CCC1*ddennydgab
          dtwt2dgab = 4.d0*(tauw/tauN)/(8.d0*dens*tauN)
          deggadgab = dneGGAdg(2)/dens



c UP UP UP UP UP UP UP UP UP UP UP UP UP UP UP \/

c    set up values to call PBE subroutine as FSP system
c      n_up
 

c do we need to initialize delc? yes
c do we need to initialize dsmdna/b? yes

         delc=0.d0
         dsmdna=0.d0
         dsmdnb=0.d0



      if (rho(n,2).gt.tol_rho) then  !check n_up density

         rhovalF  = rho(n,2)  ! set total density to n_up
         rho_A(1) = rho(n,2)   ! total   equals alpha
         rho_A(2) = rho(n,2)   ! alpha   equals n_up
         rho_A(3) = 0.d0               ! beta  equals zero
c do delrho
            delrho_A(1,1) = delrho(n,1,1)  ! nabla n_up x  
            delrho_A(2,1) = delrho(n,2,1)  ! nabla n_up y
            delrho_A(3,1) = delrho(n,3,1)  ! nabla n_up z
     
        
            delrho_A(1,2) = 0.d0    ! set beta gradient to zero
            delrho_A(2,2) = 0.d0    ! set beta gradient to zero
            delrho_A(3,2) = 0.d0    ! set beta gradient to zero


         neFSP = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE n_a*eps^FSP
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE

         Fnq = 1
         Fqwght = 1.0d0
         ipolF=2

            facF = 1.0d0 
            lfacF = .false.  
            nlfacF = .true. 

            call xc_cMpbe96(tol_rho, facF, lfacF, nlfacF,
     &         rho_A, delrho_A,
     &         dneFSPdn,dneFSPdg, 
     &          Fnq, ipolF, neFSP, Fqwght,
     &         .false.,dummy)


c         if(abs(neFSP).lt.tol_rho) goto 22
c
c        epsilon_c = n*(epsilon_c^GGA) / n
c
         pbeup = neFSP/rhovalF

C  functional deriv information for na, nb



c iffffffffffffffffffffffffffffffffffffffffffffffff
         if (pbeup.lt.pbe) then

            delc=  xx2*up*pbe !pbe
            delcnt= up*pbe/dens
             dsmdna=(up*decggadna+pbe-up*pbe/dens)/dens
             dsmdnb=(up*decggadnb-up*pbe/dens)/dens
             dsmdgaa=up*deggadgaa/dens
             dsmdgbb=up*deggadgbb/dens
             dsmdgab=up*deggadgab/dens
         else

            delc=  xx2*up*pbeup
            delcnt= up*pbeup/dens
              dsmdna=(dneFSPdn(1) - neFSP/dens)/dens
              dsmdnb= -neFSP/dens**2
              dsmdgaa=dneFSPdg(1)/dens
              dsmdgbb=dneFSPdg(3)/dens
              dsmdgab=dneFSPdg(2)/dens
         endif

c iffffffffffffffffffffffffffffffffffffffffffffffff

      endif  !end here if n_up is too small


c  DOWN DOWN  DOWN  DOWN  DOWN  DOWN  DOWN

c    set up values to call PBE subroutine as FSP system
c      n_down

       if (rho(n,3).gt.tol_rho) then  !check n_dn density


         rhovalF = rho(n,3)  ! set n_beta to total density 
         rho_A(1) = rho(n,3)    ! total   equals n_up
         rho_A(2) = rho(n,3)    ! alpha   equals n_up
         rho_A(3) = 0.d0        ! beta  equals zero

         
c do delrho
            delrho_A(1,1) = delrho(n,1,2)  ! nabla n_down x  
            delrho_A(2,1) = delrho(n,2,2)  ! nabla n_down y
            delrho_A(3,1) = delrho(n,3,2)  ! nabla n_down z
     
        
            delrho_A(1,2) = 0.d0    ! set beta gradient to zero
            delrho_A(2,2) = 0.d0    ! set beta gradient to zero
            delrho_A(3,2) = 0.d0    ! set beta gradient to zero


         neFSP = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0

         Fnq = 1
         Fqwght = 1.0d0
         ipolF=2

            facF = 1.0d0 
            lfacF = .false.  
            nlfacF = .true. 

            call xc_cMpbe96(tol_rho, facF, lfacF, nlfacF,
     &         rho_A, delrho_A,
     &         dneFSPdn,dneFSPdg, 
     &          Fnq, ipolF, neFSP, Fqwght,
     &         .false.,dummy)

c
c        epsilon_c = n*(epsilon_c^GGA) / n
c
         pbedn = neFSP/rhovalF

c iffffffffffffffffffffffffffffffffffffffffffffffff

         if (pbedn.lt.pbe) then
              delc= delc + xx2*dn*pbe 
              delcnt=delcnt + dn*pbe/dens
             dsmdna=dsmdna+(decggadna-pbe/dens)*up/dens  
             dsmdnb=dsmdnb+(dn*decggadnb+pbe-dn*pbe/dens)/dens
             dsmdgaa=dsmdgaa + dn*deggadgaa/dens
             dsmdgbb=dsmdgbb + dn*deggadgbb/dens
             dsmdgab=dsmdgab + dn*deggadgab/dens
         else
              delc= delc + xx2*dn*pbedn
              delcnt=delcnt+ dn*pbedn/dens
              dsmdna=dsmdna-dn*pbedn/dens**2
              dsmdnb=dsmdnb+(dneFSPdn(1) - dn*pbedn/dens)/dens
              dsmdgaa= dsmdgaa + dneFSPdg(1)/dens
              dsmdgbb= dsmdgbb + dneFSPdg(3)/dens
              dsmdgab= dsmdgab + dneFSPdg(2)/dens
         endif
c iffffffffffffffffffffffffffffffffffffffffffffffff



       endif ! end here if n_dn is small

         delc = -(1.d0 + CCC)*delc/dens
         revsic1 = en + delc
         revsic2 = revsic1*(1.d0 + dd*revsic1*xx3)

           Ec = Ec + revsic2*rhoval*qwght(n)*fac

c below n_a and n_b nnnnnnnnnnnnnnnnnnnnnnnnnnnn
         quima(1)=decggadna
         quima(2)=pbe*CCC*dtwt2dna + CCC*xx2*decggadna +
     &            xx2*pbe*dCCCdna
         quima(3)=xx2*dsmdna+delcnt*dtwt2dna
         quima(4)=CCC*xx2*dsmdna+xx2*delcnt*dCCCdna +
     &            delcnt*CCC*dtwt2dna

         quimb(1)=decggadnb
         quimb(2)=pbe*CCC*dtwt2dnb + CCC*xx2*decggadnb +
     &            xx2*pbe*dCCCdnb
         quimb(3)=xx2*dsmdnb+delcnt*dtwt2dnb
         quimb(4)=CCC*xx2*dsmdnb+xx2*delcnt*dCCCdnb +
     &            delcnt*CCC*dtwt2dnb


         drevdna=quima(1)+quima(2)-quima(3)-quima(4)
         drevdnb=quimb(1)+quimb(2)-quimb(3)-quimb(4)

         dtwt3dna= -3.d0*xx3/dens
         dtwt3dnb=dtwt3dna

         psea(1)=dens*drevdna+revsic1
         psea(2)=(dens*dtwt3dna*revsic1**2 +
     &    xx3*dens*2.d0*revsic1*drevdna + xx3*revsic1**2)*dd

         pseb(1)=dens*drevdnb+revsic1
         pseb(2)=(dens*dtwt3dnb*revsic1**2 +
     &    xx3*dens*2.d0*revsic1*drevdnb + xx3*revsic1**2)*dd

         Amat(n,1)=Amat(n,1) +(psea(1)+psea(2))*fac



         Amat(n,2)=Amat(n,2) +(pseb(1)+pseb(2))*fac


c below gamma_AA  ggggggggggggggggggggggggggg

         quimga(1)=deggadgaa
         quimga(2)=pbe*CCC*dtwt2dgaa + CCC*xx2*deggadgaa +
     &            xx2*pbe*dCCCdgaa
         quimga(3)=xx2*dsmdgaa + delcnt*dtwt2dgaa
         quimga(4)=CCC*xx2*dsmdgaa+xx2*delcnt*dCCCdgaa +
     &            delcnt*CCC*dtwt2dgaa


         drevdgaa=quimga(1)+quimga(2)-quimga(3)-quimga(4)
          dtwt3dgaa = (3.d0*(tauw/tauN)**2)/(8.d0*dens*tauN)

         psega(1)=dens*drevdgaa
         psega(2)=dd*dens*(revsic1*revsic1*dtwt3dgaa +
     &            xx3*2.d0*revsic1*drevdgaa)

         
         Cmat(n,1)=Cmat(n,1) +(psega(1)+psega(2))*fac


c below gamma_BB  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

         quimgb(1)=deggadgbb
         quimgb(2)=pbe*CCC*dtwt2dgbb + CCC*xx2*deggadgbb +
     &            xx2*pbe*dCCCdgbb
         quimgb(3)=xx2*dsmdgbb + delcnt*dtwt2dgbb
         quimgb(4)=CCC*xx2*dsmdgbb+xx2*delcnt*dCCCdgbb +
     &            delcnt*CCC*dtwt2dgbb

         drevdgbb=quimgb(1)+quimgb(2)-quimgb(3)-quimgb(4)
          dtwt3dgbb = dtwt3dgaa
         psegb(1)=dens*drevdgbb
         psegb(2)=dd*dens*(revsic1*revsic1*dtwt3dgbb +
     &            xx3*2.d0*revsic1*drevdgbb)
         
         Cmat(n,3)=Cmat(n,3) +(psegb(1)+psegb(2))*fac



c below gamma_AB  gggggggggggggggggggggggggggAAAAAAAAABBBBBBBBB
         quimgab(1)=deggadgab
         quimgab(2)=pbe*CCC*dtwt2dgab + CCC*xx2*deggadgab +
     &            xx2*pbe*dCCCdgab
         quimgab(3)=xx2*dsmdgab + delcnt*dtwt2dgab
         quimgab(4)=CCC*xx2*dsmdgab+xx2*delcnt*dCCCdgab +
     &            delcnt*CCC*dtwt2dgab

         drevdgab=quimgab(1)+quimgab(2)-quimgab(3)-quimgab(4)
          dtwt3dgab = 6.d0*xx2/(8.d0*dens*tauN)
         psegab(1)=dens*drevdgab
         psegab(2)=dd*dens*(revsic1*revsic1*dtwt3dgab +
     &            xx3*2.d0*revsic1*drevdgab)


        Cmat(n,2)=Cmat(n,2) +(psegab(1)+psegab(2))*fac


c below is mmat d/dtau tttttttttttttttttttttttttttttttttttttttt
 
        dtwt2dt=-2.d0*(tauw**2)/tauN**3
        dtwt3dt=-3.d0*(tauw**3)/tauN**4

        drevdt=pbe*CCC*dtwt2dt-(1.d0+CCC)*dtwt2dt*delcnt

        pset(1)=dens*drevdt
        pset(2)=dd*dens*(xx3*2.d0*revsic1*drevdt +
     &                   revsic1*revsic1*dtwt3dt)

        Mmat(n,1)=Mmat(n,1)+0.5d0*(pset(1)+pset(2))*fac

        Mmat(n,2)=Mmat(n,2)+0.5d0*(pset(1)+pset(2))*fac


     


20     continue


      endif  !end ipol=2 case

      return
      end

c

      Subroutine xc_ctpss03_d2()
      call errquit(' not coded ',0,0)
      return
      end



