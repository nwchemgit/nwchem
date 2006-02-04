#include "dft2drv.fh"
      Subroutine xc_cpkzb99(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat)     
     

c
c$Id: xc_cpkzb99.F,v 1.6 2006-02-04 04:56:08 edo Exp $
c

c     References:
c     [a] J.P. Perdew, S. Kurth, A. Zupan and P. Blaha,
c         PRL 82, 2544 (1999).

      Implicit none
c
c
c     Input and other parameters
c
      
      integer ipol, nq
      double precision dummy(1)

      double precision cfac
      logical lcfac, nlcfac
      logical  ldew
      double precision func(*)

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
      double precision delrho(nq,3,ipol), gammaval
      
c
c     Kinetic Energy Density
c
      double precision tau(nq,3,ipol)
 
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

c    first sigma term
      double precision  tauN,  Weiz
      double precision OE, CCC, mt
      parameter (OE = 0.125d0) !one eighth
      parameter (CCC = 0.53d0) !cpkzb empirical parameter

c   Second call to the cPBE subroutine

      double precision  neGGA, dneGGAdn(2), dneGGAdg(3)
      double precision  facG, Gqwght, epGGA
      logical  lfacG, nlfacG
      integer Gnq
      double precision rho_t(3), delrho_t(3,2)
      double precision  gammavalA, gammavalB
      double precision  tauNA,tauNB,  WeizA, WeizB
      double precision  sumTau, sumWeiz

c     spectral representation terms
      double precision aterma1, aterma2, aterma3
      double precision atermb1, atermb2, atermb3
      double precision aterm1, aterm2, aterm3
      double precision cterm1, cterm2, cterm3
      double precision mterm
      double precision ctermb1, ctermb2, ctermb3
      double precision gam12,pbe,tauw,xx2,xx3,en
      double precision pbeup,dtwt2dn,dtwt3dn,decggadn,dtwt2dg,
     ,     delc,dtwt3dg,decggadg,drevdn,drevdg,drevdt,
     ,     dpbeupdn,dpbeupdg,atermn,btermn,atermg,btermg,
     ,     erevc,erevsic,finaln,apartg,finalg,apartt,finalt
 

c  FULLY SPIN POLARIZED SECTION	
      double precision   rhovalF, epFSP
      double precision  neFSP, dneFSPdn(2), dneFSPdg(3)

      double precision  facF, Fqwght
      logical  lfacF, nlfacF
      integer Fnq
      double precision delrho_A(3,2), rho_A(3)


      fac = cfac
      lfac = lcfac
      nlfac = nlcfac


      if (ipol.eq.1 )then
c        ======> SPIN-RESTRICTED <======

         do 12  n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 12

         rhoval = rho(n,1)

C   set up values to call PBE subroutine
         rho_t(1) = rho(n,1)
c do delrho
         delrho_t(1,1) = delrho(n,1,1)
         delrho_t(2,1) = delrho(n,2,1)
         delrho_t(3,1) = delrho(n,3,1)
         gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &        delrho(n,2,1)*delrho(n,2,1) +
     &        delrho(n,3,1)*delrho(n,3,1)
         gam12=dsqrt(gammaval)
c
c     get E_GGA[rho,gamma]
c
         neGGA = 0.0d0  !Ec in PBE
         dneGGAdn(1) = 0.0d0   !Amat in PBE
         dneGGAdg(1) = 0.0d0  !Cmat in PBE
         dneGGAdg(2) = 0.0d0  !Cmat in PBE

         call xc_cMpbe96(tol_rho,
     &        rho_t, delrho_t,
     &        dneGGAdn,dneGGAdg, 
     &        1, ipol, neGGA)
         pbe = neGGA

         tauN = tau(n,1,1) + tau(n,2,1) + tau(n,3,1)

         tauw = 0.125d0*gammaval/rhoval

         xx2 = (tauw/tauN)**2.d0
         xx3 = (tauw/tauN)**3.d0

         en = pbe*(1.d0 + CCC*xx2)
c
c    set up values to call PBE subroutine as 
c     Fully SpinPolarized system
c

         rho_A(1) = (0.5d0)*rho(n,1)   ! total   equals (1/2)n_tot
         rho_A(2) = (0.5d0)*rho(n,1)   ! alpha   equals (1/2)n_tot
         rho_A(3) = 0.d0               ! beta  equals zero
         delrho_A(1,1) = (0.5d0)*delrho_t(1,1) ! nabla n_up x  
         delrho_A(2,1) = (0.5d0)*delrho_t(2,1) ! nabla n_up y
         delrho_A(3,1) = (0.5d0)*delrho_t(3,1) ! nabla n_up z
     
         delrho_A(1,2) = 0.d0   ! set beta gradient to zero
         delrho_A(2,2) = 0.d0   ! set beta gradient to zero
         delrho_A(3,2) = 0.d0   ! set beta gradient to zero

         neFSP = 0.0d0  !Ec in PBE
         dneFSPdn(1) = 0.0d0   !Amat in PBE
         dneFSPdn(2) = 0.0d0   !Amat in PBE

         dneFSPdg(1) = 0.0d0  !Cmat in PBE
         dneFSPdg(2) = 0.0d0  !Cmat in PBE
         dneFSPdg(3) = 0.0d0  !Cmat in PBE

c
c     get E_GGA[rho_alpha,0,gamma_alpha,0]
c
         call xc_cMpbe96(tol_rho, rho_A, delrho_A,
     &        dneFSPdn,dneFSPdg, 
     &        1, 2, neFSP)

         pbeup = neFSP

c        functional deriv info below fffffffffffff

         dtwt2dn = -2.d0*xx2/rhoval
         dtwt3dn = -3.d0*xx3/rhoval
         decggadn= dneGGAdn(1)

         dtwt2dg = 2.d0*0.125d0*tauw/(rhoval*tauN**2)

         dtwt3dg = 3.d0*xx2*0.125d0/(rhoval*tauN)
         decggadg= dneGGAdg(1)

         delc= xx2*pbeup

C  eps-tilda is eps^FSP
C  functional deriv info below fffffffffffffffff

        dpbeupdn = 0.5d0*dneFSPdn(1)
c  above note the .5's.  you are taking deriv wrt total density n
c                        not deriv wrt n_up
        dpbeupdg = 0.25d0*dneFSPdg(1)
c  note .25 above is because you want gamma=deln_tot*deln_tot


        atermn=pbe*CCC*dtwt2dn+(1.d0+CCC*xx2)*decggadn
        btermn=(1.d0+CCC)*(xx2*dpbeupdn +
     &             pbeup*dtwt2dn)
        drevdn=atermn - btermn

        atermg=pbe*ccc*dtwt2dg+(1.d0+CCC*xx2)*decggadg
        btermg=(1.d0+CCC)*(xx2*dpbeupdg+pbeup*dtwt2dg)
        drevdg=atermg-btermg
        
        drevdt=(CCC*pbe-(1.d0+CCC)*pbeup)*xx2*(-2.d0/tauN)

        
       delc = -(1.d0 + CCC)*delc
       erevc = en + delc
       erevsic = erevc
         
       if(ldew) func(n) = func(n) + rhoval*erevsic*fac
       Ec = Ec + rhoval*erevsic*qwght(n)*fac

c     derivs wrt n
       finaln= rhoval*drevdn + erevc 

       Amat(n,1)=Amat(n,1)+(finaln)*fac

c     derivs wrt g
         apartg=rhoval*drevdg
         finalg=apartg

         Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+ 2d0*finalg*fac

c     derivs wrt t
         apartt=rhoval*drevdt

         finalt=apartt

         Mmat(n,1)=Mmat(n,1)+0.5d0*finalt*fac

12     continue

c
c                                        HERE
c         ***** SPIN POLARIZATION  *****
       else   !ipol=2 and do alpha beta cases

         do 20 n = 1, nq
c             tolA
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
         dneGGAdn(1) = 0.0d0   !Amat in PBE
         dneGGAdn(2) = 0.0d0   !Amat in PBE


         dneGGAdg(1) = 0.0d0  !Cmat in PBE
         dneGGAdg(2) = 0.0d0  !Cmat in PBE
         dneGGAdg(3) = 0.0d0  !Cmat in PBE


         Gnq = 1
         Gqwght = 1.0d0

            facG = 1.0d0 
            lfacG = .false.  
            nlfacG = .true. 

            call xc_cMpbe96(tol_rho,
     &         rho_t, delrho_t,
     &         dneGGAdn,dneGGAdg, 
     &          Gnq, ipol, neGGA)

          if(abs(neGGA).lt.tol_rho) goto 20
c
c        epsilon_c = n*(epsilon_c^GGA) / n
c
         epGGA = neGGA

c         build GGA part
            Ec = Ec + neGGA*rhoval*qwght(n)*fac
            if(ldew) func(n)=func(n)+ negga*rhoval*fac

c
c     Alpha            ALPHA               ALPHA

            sumTau=0.d0
            sumWeiz=0.d0
            if (rho(n,2).lt.tol_rho) goto 25

c   tolB

            gammavalA = delrho(n,1,1)*delrho(n,1,1) +
     &                  delrho(n,2,1)*delrho(n,2,1) +
     &                  delrho(n,3,1)*delrho(n,3,1)

           WeizA = OE*gammavalA/rho(n,2)
 

            tauNA =  tau(n,1,1) + 
     &               tau(n,2,1) + 
     &               tau(n,3,1)   


            sumTau=tauNA
            sumWeiz=WeizA
c
c     Beta               BETA           BETA
c
25         continue

c take care of case with no BETA electron, ie hydrogen

c        tolC

        if (rho(n,3).lt.tol_rho) goto 30  !use this for things like H

            gammavalB = delrho(n,1,2)*delrho(n,1,2) +
     &                  delrho(n,2,2)*delrho(n,2,2) +
     &                  delrho(n,3,2)*delrho(n,3,2)
 
          WeizB = OE*gammavalB/rho(n,3)



            tauNB = tau(n,1,2) + 
     &              tau(n,2,2) + 
     &              tau(n,3,2)   

            
           sumTau = sumTau + tauNB
           sumWeiz = sumWeiz + WeizB


30      continue

             mt=rhoval*epGGA*CCC*( sumWeiz/sumTau )**2

            Ec = Ec + mt*qwght(n)*fac
            if(ldew) func(n)=func(n)+ mt*fac

c           if(rho(n,2).lt.tol_rho) goto199

        aterma1=dneGGAdn(1)
        aterma2=dneGGAdn(1)*CCC*( sumWeiz/sumTau )**2 
      aterma3=-neGGA*2.d0*CCC*(sumWeiz/sumTau)*(WeizA/sumTau)/rho(n,2)  

        cterm1=dneGGAdg(1)
        cterm2=dneGGAdg(1)*CCC*( sumWeiz/sumTau )**2
      cterm3=neGGA*2.d0*CCC*(sumWeiz/sumTau)/(sumTau*8.d0*rho(n,2))
         mterm=-2.d0*CCC*neGGA*(sumWeiz**2)/sumTau**3
      
           Amat(n,1)= Amat(n,1)+(aterma1+aterma2+aterma3)*fac  
           Cmat(n,1)= Cmat(n,1)+(cterm1+cterm2+cterm3)*fac
           Mmat(n,1)=Mmat(n,1)+mterm*fac

          Cmat(n,2)=Cmat(n,2) + (dneGGAdg(2) +
     &       CCC*dneGGAdg(2)*( sumWeiz/sumTau )**2)*fac



          if(rho(n,3).gt.tol_rho) then

        atermb1=dneGGAdn(2)
        atermb2=dneGGAdn(2)*CCC*( sumWeiz/sumTau )**2 
      atermb3=-neGGA*2.d0*CCC*(sumWeiz/sumTau)*(WeizB/sumTau)/rho(n,3)    
      
           Amat(n,2)= Amat(n,2)+(atermb1+atermb2+atermb3)*fac  

         ctermb1=dneGGAdg(3)
        ctermb2=dneGGAdg(3)*CCC*( sumWeiz/sumTau )**2
      ctermb3=neGGA*2.d0*CCC*(sumWeiz/sumTau)/(sumTau*8.d0*rho(n,3))
           Cmat(n,3)= Cmat(n,3)+(ctermb1+ctermb2+ctermb3)*fac

        mterm=-2.d0*CCC*neGGA*(sumWeiz**2)/sumTau**3
           Mmat(n,2)=Mmat(n,2)+mterm*fac


         else
           Amat(n,2)=0.d0
           Cmat(n,3)=0.d0
           Mmat(n,2)=0.d0
         endif

20      continue

       endif

C  EEEEEEEEEEEEEEEEEE  End of first half of equation (15)
C                      no matter the system, still need to run a
C                      Fully Spin Polarized (FSP) calc below
C                      on both unpolarized and polarized systems


c PPPPPPPPPPPPPPPPPPPP POLARIZED PPPPPPPPPPPPPPPPPPPPPP
         if (ipol.eq.2) then

C  FFFFFFFFFFFFFFF  FULLY SPIN POLARIZED PART FFFFFFFFFFFFFFFFFFFFFF

c  ipol still equals 2 for these cases, 
c  but explicitly  set things here

        do 50 n = 1, nq
c           tolD

             if (rho(n,2).lt.tol_rho) goto 50  


c     Alpha            ALPHA               ALPHA
         rhovalF = rho(n,2)  ! set n_alpha to total density 
         rho_A(1) = rho(n,2)    ! total   equals n_up
         rho_A(2) = rho(n,2)    ! alpha   equals n_up
         rho_A(3) = 0.d0        ! beta  equals zero

         
c do delrho
            delrho_A(1,1) = delrho(n,1,1)  ! nabla n_up x  
            delrho_A(2,1) = delrho(n,2,1)  ! nabla n_up y
            delrho_A(3,1) = delrho(n,3,1)  ! nabla n_up z
     
        
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

            facF = 1.0d0 !MUST BE a number now
            lfacF = .false.  !MUST BE a T/F now  always false
            nlfacF = .true. !MUST BE a T/F now

            call xc_cMpbe96(tol_rho, 
     &         rho_A, delrho_A,
     &         dneFSPdn,dneFSPdg, 
     &          1, 2, neFSP)


          if(abs(neFSP).lt.tol_rho) goto 50
c
c        epsilon_c = n*(epsilon_c^GGA) / n
c
         epFSP = neFSP
c
            gammavalA = delrho(n,1,1)*delrho(n,1,1) +
     &                  delrho(n,2,1)*delrho(n,2,1) +
     &                  delrho(n,3,1)*delrho(n,3,1)

           WeizA = OE*gammavalA/rho(n,2)

c   
           tauNA =  tau(n,1,1) + tau(n,2,1) + tau(n,3,1)   
            
            mt=(CCC + 1.d0)*rhovalF*epFSP*( WeizA/TauNA )**2

            Ec = Ec - mt*qwght(n)*fac
            if(ldew) func(n)=func(n) - mt*fac

            aterma1=dneFSPdn(1)*( WeizA/TauNA )**2
            aterma2= -2.d0*neFSP*((WeizA/TauNA )**2)/rhovalF

            cterm1=dneFSPdg(1)*( WeizA/TauNA )**2
            cterm2=neFSP*2.d0*((WeizA/TauNA )**2)/gammavalA

            mterm= -2.d0*neFSP*((WeizA/TauNA )**2)/TauNA

           Amat(n,1)=Amat(n,1)-(1.d0+CCC)*(aterma1+aterma2)*fac
          Cmat(n,1)=Cmat(n,1)-(1.d0+CCC)*(cterm1+cterm2)*fac
          Mmat(n,1)=Mmat(n,1)-(1.d0+CCC)*mterm*fac

50      continue

c
c     Beta               BETA           BETA
c

         do 75 n = 1, nq

        if (rho(n,3).lt.tol_rho) goto 75  !use this for things like H


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

         call xc_cMpbe96(tol_rho, 
     &        rho_A, delrho_A,
     &        dneFSPdn,dneFSPdg, 
     &        1, 2, neFSP)


         if(abs(neFSP).lt.tol_rho) goto 75  ! check me
c
c        epFSP = n*(epsilon_c^GGA) / n
c
         epFSP = neFSP


            gammavalB = delrho(n,1,2)*delrho(n,1,2) +
     &                  delrho(n,2,2)*delrho(n,2,2) +
     &                  delrho(n,3,2)*delrho(n,3,2)

 
          WeizB = OE*gammavalB/rho(n,3)

c   do we have to take care of cases where tauN blows up in denom??
c    must this be protected against

            tauNB = tau(n,1,2) + tau(n,2,2) + tau(n,3,2)   


             mt=(CCC + 1.d0)*rhovalF*epFSP*( WeizB/TauNB )**2

            Ec = Ec - mt*qwght(n)*fac
            if(ldew) func(n)=func(n) - mt*fac

            atermb1=dneFSPdn(1)*( WeizB/TauNB )**2
            atermb2= -2.d0*neFSP*((WeizB/TauNB )**2)/rhovalF

           Amat(n,2)=Amat(n,2)-(1.d0+CCC)*(atermb1+atermb2)*fac

            cterm1=dneFSPdg(1)*( WeizB/TauNB )**2
            cterm2=neFSP*2.d0*((WeizB/TauNB )**2)/gammavalB

          Cmat(n,3)=Cmat(n,3)-(1.d0+CCC)*(cterm1+cterm2)*fac

            mterm= -2.d0*neFSP*((WeizB/TauNB )**2)/TauNB
          Mmat(n,2)=Mmat(n,2)-(1.d0+CCC)*mterm*fac

75      continue

          endif



      return
      end

c

      Subroutine xc_cpkzb99_d2()
      call errquit(' not coded ',0,0)
      return
      end




