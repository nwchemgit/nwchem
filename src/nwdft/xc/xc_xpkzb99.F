c     Perdew-Kurth-Zupan-Blaha '99 exchange functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                 cor only:    tauW - von Weiszacker kinetic energy density
c     References:
c     [a] J.P. Perdew, S. Kurth, A. Zupan and P. Blaha,
c         PRL 82, 2544 (1999).


      Subroutine xc_xpkzb99(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, tau,Mmat)


c   
c$Id$
c
      implicit none
c
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision tol_rho, pi
c
      integer n
      double precision rrho, rho43, rho13, gamma
c
c     kinetic energy density   or  tau
c
      double precision tau(nq,ipol), Mmat(nq,*)
      double precision tauN

      double precision  p, qtil, x
      double precision   rho53, rho83,  mt
      double precision   F83, F23, F53, F43, F13
      double precision   G920 
      double precision   bigD, uk, ruk
      double precision    CX1, CX2, CX3, CX4
      double precision   P32, Ax
c     functional derivatives below
      double precision   d1p(3), d1x(3)
      double precision   rg2, d1mt(3)
      double precision   d1qtil(3), T1
c     functional derivatives above

      parameter(uk=0.8040d0, bigD=0.113d0, ruk=1.d0/uk)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
      parameter (F83=8.d0/3.d0, F23=2.d0/3.d0, F53=5.d0/3.d0)
      parameter (G920  =9.d0/20.d0 )
      parameter (CX1  =  10.d0/81.d0, 
     &     CX2  = 146.d0/2025.d0,
     &     CX3  = -73.d0/405.d0 )
      parameter (CX4= (bigD + ruk*CX1**2) )
c     
      pi=acos(-1d0)
      Ax = (-0.75d0)*(3d0/pi)**F13
      P32 = 1.d0/( 2.d0*(3.d0*pi**2)**(F23) ) 
      
c
      if (ipol.eq.1 )then
c
c     ======> SPIN-RESTRICTED <======
c     

         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10

c  rho43= n*e_x^unif=exchange energy per electron for uniform electron gas
c       = n* Ax*n^(1/3)   or n*C*n^(1/3) 

            rho43 = Ax*rho(n,1)**F43  ! Ax*n^4/3
            rrho = 1d0/rho(n,1)       ! reciprocal of rho
            rho13 = F43*rho43*rrho   !functional deriv of rho43

            rho53 = rho(n,1)**F53
            rho83 = rho(n,1)**F83


C  Below we just sum up the LDA contribution to the functional
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac
                Amat(n,1) = Amat(n,1) + rho13*fac
               if(ldew)func(n) = func(n) + rho43*fac              
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
c            gam12 = dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 10
            tauN = tau(n,1)

            p=P32*gamma/(rho83*2.d0)
            qtil=(3.d0*tauN*P32/rho53)-G920-(p/12.d0)
c
c     Evaluate the GC part of Fx, i.e. mt = Fx(p,qtil) - 1
c

            x= CX1*p + CX2*qtil*qtil + CX3*qtil*p+ CX4*p*p

             if (.not.(nlfac.and.x.gt.tol_rho)) goto 10

            
            mt = uk - uk/(1.d0 + x*ruk)


C      functional derivatives

             rg2=1.d0/( (1.d0 + x*ruk)*(1.d0 + x*ruk) )
             
c    deriv wrt n, the density (for Amat)
             d1p(1) = -F83*rrho*p

             T1=3.d0*P32*tauN*(1.d0/rho53)
             d1qtil(1) = -F53*T1*rrho - d1p(1)/12.d0

             d1x(1) = CX1*d1p(1) + CX2*2.d0*qtil*d1qtil(1) +
     &       CX3*(qtil*d1p(1) + p*d1qtil(1)) + CX4*2d0*p*d1p(1)

             d1mt(1) = rg2*d1x(1)


c     deriv wrt gamma, the gradient (for Cmat)
             d1p(2) = 0.5d0*P32/rho83
             d1qtil(2) = -d1p(2)/12.d0

             d1x(2) = CX1*d1p(2) + CX2*2d0*qtil*d1qtil(2) +
     &       CX3*(qtil*d1p(2) + p*d1qtil(2)) + CX4*2d0*p*d1p(2)

             d1mt(2) = rg2*d1x(2)


c     deriv wrt tau, the Kinetic Energy Density (for Mmat)

             d1p(3) = 0.d0
c             d1qtil(3) = -d1p(2)/12.d0  am sure this is wrong

             d1qtil(3) = 3.d0*P32/rho53

             d1x(3) = CX2*2.d0*qtil*d1qtil(3) +
     &       CX3*p*d1qtil(3)

             d1mt(3) = rg2*d1x(3)



C    Below we add the MetaGGA correction to the LDA part from above

             if(ldew)func(n) = func(n) + rho43*mt*fac              
             Ex = Ex + rho43*mt*qwght(n)*fac

             Amat(n,1) =Amat(n,1) + (mt*rho13 + rho43*d1mt(1))*fac

             Cmat(n,1) = Cmat(n,1) + 2.d0*(rho43*d1mt(2)*fac)
c                  check on this two or one

             Mmat(n,1) = Mmat(n,1) +0.5d0*(rho43*d1mt(3)*fac)


 10      continue

      else
c
c        ======> SPIN-UNRESTRICTED <======

c
c  use spin density functional theory ie n-->2n
c                               Ex=(1/2)Ex[2*alpha] + (1/2)Ex[2*beta]

         do 20 n = 1, nq
             if (rho(n,1).lt.tol_rho) goto 20
c
c     Alpha            ALPHA               ALPHA
c
            if (rho(n,2).lt.tol_rho) goto 25
            rho43 = Ax*(2d0*rho(n,2))**F43 ! spin scaled
            rrho = 0.5d0/rho(n,2)          ! spin scaled
            rho13 = F43*rho43*rrho  !spin scaled & (1/2)factor

            rho53 = (2.d0*rho(n,2))**F53   ! spin scaled
            rho83 = (2.d0*rho(n,2))**F83   ! spin scaled
c  note all the "rho"  quantities ARE spin scaled, for later use

            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
                Amat(n,1) = Amat(n,1) + rho13*fac
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
c NOTE  gamma above  is not spin scaled.  that is why
c                           -there is 4.d0*gamma in p
c                           -there is 2.0 in the gam12 term



c            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 25

c   below note factor of two for spin scaling
            tauN = 2.d0* tau(n,1)
            
c
c     Evaluate the GC part of Fx, i.e. mt(x) = Fx - 1
c           

            p=0.5d0*P32*(4.d0*gamma)/rho83
            qtil=(3.d0*tauN*P32/rho53) - G920 - (p/12.d0)
 

           x= CX1*p + CX2*qtil*qtil + CX3*qtil*p+ CX4*p*p
            if (.not.(nlfac.and.x.gt.tol_rho)) goto 25

             rg2=1.d0/( (1.d0 + x*ruk)*(1.d0 + x*ruk) )

c  ccccccc   deriv wrt n, the density

            d1p(1) = p*(-F83)*(2.d0*rrho)  ! spin scaled

            T1=3.d0*P32*tauN/rho53        ! spin scaled
            d1qtil(1) = -F53*T1*2.d0*rrho - d1p(1)/12.d0

             d1x(1) = CX1*d1p(1) + CX2*2.d0*qtil*d1qtil(1) +
     &       CX3*(qtil*d1p(1) + p*d1qtil(1)) + CX4*2d0*p*d1p(1)

             d1mt(1) = rg2*d1x(1)


c     deriv wrt gamma, the gradient
             d1p(2) = 0.5d0*P32*4.d0/rho83     ! spin scaled
             d1qtil(2) = -d1p(2)/12.d0        !spin scaled

             d1x(2) = CX1*d1p(2) + CX2*2d0*qtil*d1qtil(2) +
     &       CX3*(qtil*d1p(2) + p*d1qtil(2)) + CX4*2d0*p*d1p(2)

             d1mt(2) = rg2*d1x(2)


c     deriv wrt tau, the Kinetic Energy Density

c             d1p(3) = 0.d0  term shown for completeness
             d1qtil(3) = 3.d0*P32/rho53

             d1x(3) = CX2*2.d0*qtil*d1qtil(3) +
     &       CX3*p*d1qtil(3)

             d1mt(3) = rg2*d1x(3)

            mt = uk - uk/(1.d0 + x*ruk)

            Ex = Ex + rho43*mt*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*mt*fac*0.5d0

            Amat(n,1)=Amat(n,1)+(mt*rho13 + 0.5d0*rho43*d1mt(1))*fac
c       note that the (.5) is built into the rho13 term already
c       hence we only need to put it onto the second term in Amat


            Cmat(n,1) = Cmat(n,1) + (0.5d0*rho43*d1mt(2)*fac)
            Mmat(n,1) = Mmat(n,1) +1.0d0*(0.5d0*rho43*d1mt(3)*fac)

          
c
c     Beta               BETA           BETA
c
 25         continue

            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = Ax*(2d0*rho(n,3))**F43  ! Ax (2 nBeta)^4/3
            rrho = 0.5d0/rho(n,3)           ! 1/(2 nBeta)
            rho13 = F43*rho43*rrho   !spin scaled func deriv of rho43  

            rho53 = (2.d0*rho(n,3))**F53
            rho83 = (2.d0*rho(n,3))**F83

C  note all "rho" quantities above are spin scaled for later use

            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
                Amat(n,2) = Amat(n,2) + rho13*fac
              if(ldew)func(n) = func(n) + rho43*fac*0.5d0
            endif
c

            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
c NOTE  gamma above  is not spin scaled.  that is why
c                           -there is 4.d0*gamma in p term
c                           -there is 2.0 in the gam12 term


c            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 20

c   below note factor of two for spin scaling
            tauN = 2.d0* tau(n,2)

            
c
c     Evaluate the GC part of F(x), i.e. mt(x) = Fx - 1
c
         

            p=0.5d0*P32*(4.d0*gamma)/rho83
            qtil=(3.d0*tauN*P32/rho53) - G920 - (p/12.d0)


           x= CX1*p + CX2*qtil*qtil + CX3*qtil*p+ CX4*p*p


            if (.not.(nlfac.and.x.gt.tol_rho)) goto 20


             rg2=1.d0/( (1.d0 + x*ruk)*(1.d0 + x*ruk) )

c  ccccccc   deriv wrt n, the density

            d1p(1) = p*(-F83)*(2.d0*rrho)  ! spin scaled

            T1=3.d0*P32*tauN/rho53        ! spin scaled
            d1qtil(1) = -F53*T1*2.d0*rrho - d1p(1)/12.d0

             d1x(1) = CX1*d1p(1) + CX2*2.d0*qtil*d1qtil(1) +
     &       CX3*(qtil*d1p(1) + p*d1qtil(1)) + CX4*2d0*p*d1p(1)

             d1mt(1) = rg2*d1x(1)


c     deriv wrt gamma, the gradient
             d1p(2) = 0.5d0*P32*4.d0/rho83     ! spin scaled
             d1qtil(2) = -d1p(2)/12.d0        !spin scaled

             d1x(2) = CX1*d1p(2) + CX2*2d0*qtil*d1qtil(2) +
     &       CX3*(qtil*d1p(2) + p*d1qtil(2)) + CX4*2d0*p*d1p(2)

             d1mt(2) = rg2*d1x(2)



c     deriv wrt tau, the Kinetic Energy Density

C             d1p(3) = 0.d0  included for completeness
             d1qtil(3) = 3.d0*P32/rho53


             d1x(3) = CX2*2.d0*qtil*d1qtil(3) +
     &       CX3*p*d1qtil(3)

             d1mt(3) = rg2*d1x(3)
             
             mt = uk - uk/(1.d0 + x*ruk)
            
             Ex = Ex + rho43*mt*qwght(n)*fac*0.5d0
             if(ldew)func(n) = func(n) + rho43*mt*fac*0.5d0

             Amat(n,2)=Amat(n,2)+(mt*rho13 + 0.5d0*rho43*d1mt(1))*fac
c       note that the (.5) is built into the rho13 term already
c       hence we only need to put it onto the second term in Amat

            Cmat(n,3) = Cmat(n,3) + 0.5d0*rho43*d1mt(2)*fac
            Mmat(n,2) = Mmat(n,2) + 0.5d0*rho43*d1mt(3)*fac

20      continue
      endif
c
      return
      end




      Subroutine xc_xpkzb99_d2()
      call errquit(' xpkzb99: d2 not coded ',0,0)
      return
      end


