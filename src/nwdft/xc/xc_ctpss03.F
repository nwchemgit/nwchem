
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
c$Id: xc_ctpss03.F,v 1.1 2003-12-18 21:53:35 pollack Exp $
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

      double precision THRD, FACE
      double precision zeta, ccc, fer
      double precision tauwup, tauwdn, tauw
      double precision xx2,xx3, ter1,ter2,ter3
      double precision sumz, grz, szeta2, term1, term2
      double precision yyz, pbe, en
      double precision aa,bb,cc,dd,dnu
      double precision up,dn,grup, grdn, gr, grupdn
      double precision dens
      double precision pbeup, pbedn, delc, revsic

      parameter ( THRD = 0.33333333333333333333333333333d0 )
      parameter ( FACE = 3.09366772628013593097d0 )
c      FACE = (3*PI**2)**(1/3)
 
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
c              NO SPIN POLARIZATION	
c        ======> SPIN-RESTRICTED <======
c                     or
c                SPIN-UNPOLARIZED
c
c

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
c
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


c         zeta = (up - dn)/dens =0 for ipol=1
        zeta=0.d0

c         CCC = 0.53d0 + aa*zeta**2.d0 + bb*zeta**4.d0 + cc*zeta**6.d0
         CCC = 0.53d0 !since zeta=0
         fer = FACE*dens**THRD


         tauw = 0.125d0*gammaval/dens

c         tau = tauup + taudn = tauN

         xx2 = (tauw/tauN)**2.d0
         xx3 = (tauw/tauN)**3.d0
c----------------------------
c         ter1 = (1.d0 - zeta)**2.d0*grup**2.d0
c         ter2 = (1.d0 + zeta)**2.d0*grdn**2.d0
c         ter3 = 2.d0*(1.d0 - zeta**2.d0)*grupdn
c         sumz = ter1 + ter2 - ter3
c         sumz = dabs(sumz)
c         grz = dsqrt(sumz)/dens
          grz=0.d0 !zeta =0 so necessarily follows
         szeta2 = (grz/(2.d0*fer))**2.d0

         if (dabs(zeta).lt.0.9999d0) then
            term1 = 1.d0/(1.d0 + zeta)**(4.d0/3.d0)
            term2 = 1.d0/(1.d0 - zeta)**(4.d0/3.d0) 
            yyz = szeta2*(term1 + term2)/2.d0
         else
            yyz = 0.d0
         end if
         CCC = CCC/(1.d0 + dnu*yyz)**4.d0

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



c do we need to initialize delc? NO
c         delc=0.d0

         if (pbeup.lt.pbe) then
              delc= xx2*pbe
         else
              delc= xx2*pbeup
         endif

         delc = -(1.d0 + CCC)*delc
         revsic = en + delc
         revsic = revsic*(1.d0 + dd*revsic*xx3)

             Ec = Ec + rhoval*revsic*qwght(n)*fac

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
c         if(zeta.lt.-1.0d0) zeta=-1.0d0
c         if(zeta.gt. 1.0d0) zeta= 1.0d0


         CCC = 0.53d0 + aa*zeta**2.d0 + bb*zeta**4.d0 + cc*zeta**6.d0

         CCC = CCC/(1.d0 + dnu*yyz)**4.d0

          en = pbe*(1.d0 + CCC*xx2)

c UP UP UP UP UP UP UP UP UP UP UP UP UP UP UP \/

c    set up values to call PBE subroutine as FSP system
c      n_up
 

c do we need to initialize delc? yes
         delc=0.d0



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


         if (pbeup.lt.pbe) then
              delc=  xx2*up*pbe !pbe
c              delc= delc + xx2*up*pbe !pbe
         else
              delc=  xx2*up*pbeup
c              delc= delc + xx2*up*pbeup
         endif

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


         if (pbedn.lt.pbe) then
              delc= delc + xx2*dn*pbe   
         else
              delc= delc + xx2*dn*pbedn
         endif


       endif ! end here if n_dn is small

         delc = -(1.d0 + CCC)*delc/dens
         revsic = en + delc
         revsic = revsic*(1.d0 + dd*revsic*xx3)

           Ec = Ec + revsic*rhoval*qwght(n)*fac

20     continue


      endif  !end ipol=2 case

      return
      end

c

      Subroutine xc_ctpss03_d2()
      call errquit(' not coded ',0,0)
      return
      end



