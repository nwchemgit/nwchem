

      Subroutine xc_cpkzb99(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght,
     &                     tau, Amat, Cmat, Mmat)     
     

c
c$Id: xc_cpkzb99.F,v 1.3 2003-10-15 19:20:04 pollack Exp $
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
      double precision OE, bigC, mt
      parameter (OE = 0.125d0) !one eighth
      parameter (bigC = 0.53d0) !cpkzb empirical parameter

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
 

c  FULLY SPIN POLARIZED SECTION	
      double precision   rhovalF, epFSP
      double precision  neFSP, dneFSPdn(2), dneFSPdg(3)

      integer  ipolF
      double precision  facF, Fqwght
      logical  lfacF, nlfacF
      integer Fnq
      double precision delrho_A(3,2), rho_A(3)


      fac = cfac
      lfac = lcfac
      nlfac = nlcfac


      if (ipol.eq.1 )then
c              NO SPIN POLARIZATION	
c        ======> SPIN-RESTRICTED <======
c                     or
c                SPIN-UNPOLARIZED
c
c

         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10

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
c
            facG = 1.0d0 !MUST BE ONE
            lfacG = .false.  !MUST BE F
            nlfacG = .true. !MUST BE a T/F 

            call xc_cMpbe96(tol_rho, facG, lfacG, nlfacG,
     &         rho_t, delrho_t,
     &         dneGGAdn,dneGGAdg, 
     &          Gnq, ipol, neGGA,Gqwght,
     &         .false.,dummy)


         if(abs(neGGA).lt.tol_rho) goto 10
c
c        epGGA = n*(epsilon_c^GGA) / n =cor. energy per electron
c
         epGGA = neGGA/rhoval

c         build GGA part
             Ec = Ec + neGGA*qwght(n)*fac


            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1)

           Weiz = OE*gammaval/rho(n,1)



            tauN = tau(n,1,1) + 
     &             tau(n,2,1) + 
     &             tau(n,3,1)


             aterm1=dneGGAdn(1)
             aterm2=bigC*dneGGAdn(1)*( Weiz/tauN )**2
          aterm3=-neGGA*2.d0*bigC*(( Weiz/tauN )**2)/rhoval

             cterm1=dneGGAdg(1)
             cterm2=dneGGAdg(1)*bigC*( Weiz/tauN )**2
          cterm3=2.d0*bigC*neGGA*(Weiz/tauN)/(8.d0*rhoval*tauN)



             mt=rhoval*epGGA*bigC*( Weiz/tauN )**2

             Ec = Ec + mt*qwght(n)*fac
          Amat(n,1)=Amat(n,1)+(aterm1+aterm2+aterm3)*fac

          Cmat(n,1)=Cmat(n,1)+(cterm1+cterm2+cterm3)*fac
          Cmat(n,2)=Cmat(n,2) + (dneGGAdg(2) +
     &        (( Weiz/tauN )**2)*bigC*dneGGAdg(2))*fac


         Mmat(n,1)=Mmat(n,1)-(2.d0*bigC*neGGA*(Weiz**2)/
     &             (tauN**3))*fac



10      continue
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

            call xc_cMpbe96(tol_rho, facG, lfacG, nlfacG,
     &         rho_t, delrho_t,
     &         dneGGAdn,dneGGAdg, 
     &          Gnq, ipol, neGGA,Gqwght,
     &         .false.,dummy)

          if(abs(neGGA).lt.tol_rho) goto 20
c
c        epsilon_c = n*(epsilon_c^GGA) / n
c
         epGGA = neGGA/rhoval

c         build GGA part
            Ec = Ec + neGGA*qwght(n)*fac


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

             mt=rhoval*epGGA*bigC*( sumWeiz/sumTau )**2

            Ec = Ec + mt*qwght(n)*fac

c           if(rho(n,2).lt.tol_rho) goto199

        aterma1=dneGGAdn(1)
        aterma2=dneGGAdn(1)*bigC*( sumWeiz/sumTau )**2 
      aterma3=-neGGA*2.d0*bigC*(sumWeiz/sumTau)*(WeizA/sumTau)/rho(n,2)  

        cterm1=dneGGAdg(1)
        cterm2=dneGGAdg(1)*bigC*( sumWeiz/sumTau )**2
      cterm3=neGGA*2.d0*bigC*(sumWeiz/sumTau)/(sumTau*8.d0*rho(n,2))
         mterm=-2.d0*bigC*neGGA*(sumWeiz**2)/sumTau**3
      
           Amat(n,1)= Amat(n,1)+(aterma1+aterma2+aterma3)*fac  
           Cmat(n,1)= Cmat(n,1)+(cterm1+cterm2+cterm3)*fac
           Mmat(n,1)=Mmat(n,1)+mterm*fac

          Cmat(n,2)=Cmat(n,2) + (dneGGAdg(2) +
     &       bigC*dneGGAdg(2)*( sumWeiz/sumTau )**2)*fac



c199           if(rho(n,3).lt.tol_rho) goto20

        atermb1=dneGGAdn(2)
        atermb2=dneGGAdn(2)*bigC*( sumWeiz/sumTau )**2 
      atermb3=-neGGA*2.d0*bigC*(sumWeiz/sumTau)*(WeizB/sumTau)/rho(n,3)    
      
           Amat(n,2)= Amat(n,2)+(atermb1+atermb2+atermb3)*fac  

         ctermb1=dneGGAdg(3)
        ctermb2=dneGGAdg(3)*bigC*( sumWeiz/sumTau )**2
      ctermb3=neGGA*2.d0*bigC*(sumWeiz/sumTau)/(sumTau*8.d0*rho(n,3))
           Cmat(n,3)= Cmat(n,3)+(ctermb1+ctermb2+ctermb3)*fac

        mterm=-2.d0*bigC*neGGA*(sumWeiz**2)/sumTau**3
           Mmat(n,2)=Mmat(n,2)+mterm*fac


20      continue

       endif

C  EEEEEEEEEEEEEEEEEE  End of first half of equation (15)
C                      no matter the system, still need to run a
C                      Fully Spin Polarized (FSP) calc below
C                      on both unpolarized and polarized systems


C  FFFFFFFFFFFFFFF  FULLY SPIN POLARIZED PART BELOW FFFFFFFFFFFFFFFFFFFFFF



         if (ipol.eq.1) then

c UNNNNNNNNNNNNNNNNNN UNPOLARIZED UUUNNNNNNNNNNNNNNNNN
c  just run for alpha case
c  since beta case is exactly equal,  double final answer

C   ipol=1  ex Helium etc...


c  tho ipol=1, ipolF=2 for these cases, so 
c   we want to explicitly  set things here

            ipolF = 2


         do 22 n = 1, nq
             if (rho(n,1).lt.tol_rho) goto 22  

c     Alpha            ALPHA               ALPHA
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

            facF = 1.0d0 
            lfacF = .false.  
            nlfacF = .true. 

            call xc_cMpbe96(tol_rho, facF, lfacF, nlfacF,
     &         rho_A, delrho_A,
     &         dneFSPdn,dneFSPdg, 
     &          Fnq, ipolF, neFSP, Fqwght,
     &         .false.,dummy)


         if(abs(neFSP).lt.tol_rho) goto 22
c
c        epsilon_c = n*(epsilon_c^GGA) / n
c
         epFSP = neFSP/rhovalF


            gammavalA = 0.25d0*(delrho(n,1,1)*delrho(n,1,1) +
     &                         delrho(n,2,1)*delrho(n,2,1) +
     &                         delrho(n,3,1)*delrho(n,3,1)   )


           WeizA = OE*gammavalA/rhovalF


c   
            tauNA =0.5d0*(tau(n,1,1) + 
     &                   tau(n,2,1) + 
     &                   tau(n,3,1)   )

C FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF sprep

             mt=(bigC + 1.d0)*rhovalF*epFSP*( WeizA/TauNA )**2

            aterm1=dneFSPdn(1)*( WeizA/TauNA )**2
            aterm2=-neFSP*2.d0*((WeizA/TauNA )**2)/rhovalF

            cterm1= dneFSPdg(1)*( WeizA/TauNA )**2
            cterm2=2.d0*neFSP*( WeizA/TauNA )/(8.d0*rhovalF*TauNA)
            mterm= -2.d0*neFSP*(WeizA **2)/TauNA**3

            Ec = Ec - (2.d0)*mt*qwght(n)*fac
         Amat(n,1)=Amat(n,1)-2.d0*(1.d0+bigC)*(aterm1+aterm2)*fac
        Cmat(n,1)=Cmat(n,1)-2.d0*(1.d0+bigC)*(cterm1+cterm2)*fac

        Mmat(n,1)=Mmat(n,1)-2.d0*(1.d0+bigC)*mterm*fac

c      the 2 is because we know n_up = n_dn = n/2
c      so we can just double the results     


c  because these systems are spin unpolarized, n_up = n_down
c  and we only have to calculate e_cGGA once but multiply by 2


22      continue

c                                           /\   
c  fully spin data for unpolarized systems  ||
c                                           ||

          endif

c PPPPPPPPPPPPPPPPPPPP POLARIZED PPPPPPPPPPPPPPPPPPPPPP
         if (ipol.eq.2) then


C   ipol=2  ex hydrogen, li, nitrogen etc...
C  FFFFFFFFFFFFFFF  FULLY SPIN POLARIZED PART FFFFFFFFFFFFFFFFFFFFFF


c  ipol still equals 2 for these cases, 
c  but explicitly  set things here

            ipolF = 2

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

            call xc_cMpbe96(tol_rho, facF, lfacF, nlfacF,
     &         rho_A, delrho_A,
     &         dneFSPdn,dneFSPdg, 
     &          Fnq, ipolF, neFSP, Fqwght,
     &         .false.,dummy)


          if(abs(neFSP).lt.tol_rho) goto 50
c
c        epsilon_c = n*(epsilon_c^GGA) / n
c
         epFSP = neFSP/rhovalF



c

            gammavalA = delrho(n,1,1)*delrho(n,1,1) +
     &                  delrho(n,2,1)*delrho(n,2,1) +
     &                  delrho(n,3,1)*delrho(n,3,1)


           WeizA = OE*gammavalA/rho(n,2)


c   
            tauNA =  tau(n,1,1) + 
     &               tau(n,2,1) + 
     &               tau(n,3,1)   




             mt=(bigC + 1.d0)*rhovalF*epFSP*( WeizA/TauNA )**2

            Ec = Ec - mt*qwght(n)*fac

            aterma1=dneFSPdn(1)*( WeizA/TauNA )**2
            aterma2= -2.d0*neFSP*((WeizA/TauNA )**2)/rhovalF

            cterm1=dneFSPdg(1)*( WeizA/TauNA )**2
            cterm2=neFSP*2.d0*((WeizA/TauNA )**2)/gammavalA

            mterm= -2.d0*neFSP*((WeizA/TauNA )**2)/TauNA

           Amat(n,1)=Amat(n,1)-(1.d0+bigC)*(aterma1+aterma2)*fac
          Cmat(n,1)=Cmat(n,1)-(1.d0+bigC)*(cterm1+cterm2)*fac
          Mmat(n,1)=Mmat(n,1)-(1.d0+bigC)*mterm*fac

50      continue

C                                                /\
c                                                ||
c     Alpha            ALPHA               ALPHA ||




c
c     Beta               BETA           BETA
c


c take care of case with no BETA electron, ie hydrogen


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

         Fnq = 1
         Fqwght = 1.0d0

            facF = 1.0d0 
            lfacF = .false.  
            nlfacF = .true. 

            call xc_cMpbe96(tol_rho, facF, lfacF, nlfacF,
     &         rho_A, delrho_A,
     &         dneFSPdn,dneFSPdg, 
     &          Fnq, ipolF, neFSP, Fqwght,
     &         .false.,dummy)


         if(abs(neFSP).lt.tol_rho) goto 75  ! check me
c
c        epFSP = n*(epsilon_c^GGA) / n
c
         epFSP = neFSP/rhovalF


            gammavalB = delrho(n,1,2)*delrho(n,1,2) +
     &                  delrho(n,2,2)*delrho(n,2,2) +
     &                  delrho(n,3,2)*delrho(n,3,2)

 
          WeizB = OE*gammavalB/rho(n,3)

c   do we have to take care of cases where tauN blows up in denom??
c    must this be protected against

            tauNB = tau(n,1,2) + 
     &              tau(n,2,2) + 
     &              tau(n,3,2)   


             mt=(bigC + 1.d0)*rhovalF*epFSP*( WeizB/TauNB )**2

            Ec = Ec - mt*qwght(n)*fac

            atermb1=dneFSPdn(1)*( WeizB/TauNB )**2
            atermb2= -2.d0*neFSP*((WeizB/TauNB )**2)/rhovalF

           Amat(n,2)=Amat(n,2)-(1.d0+bigC)*(atermb1+atermb2)*fac

            cterm1=dneFSPdg(1)*( WeizB/TauNB )**2
            cterm2=neFSP*2.d0*((WeizB/TauNB )**2)/gammavalB

          Cmat(n,3)=Cmat(n,3)-(1.d0+bigC)*(cterm1+cterm2)*fac

            mterm= -2.d0*neFSP*((WeizB/TauNB )**2)/TauNB
          Mmat(n,2)=Mmat(n,2)-(1.d0+bigC)*mterm*fac

75      continue

C                                                /\
c                                                ||
c     BETA            BETA               BETA    ||
c          fully spin polarized above


          endif



      return
      end

c

      Subroutine xc_cpkzb99_d2()
      call errquit(' not coded ',0)
      return
      end




