c     Perdew-Kurth-Zupan-Blaha '99 exchange functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                      cor:    tauW - von Weiszacker kinetic energy density
c     References:
c     [a] J.P. Perdew, S. Kurth, A. Zupan and P. Blaha,
c         PRL 82, 2544 (1999).


      Subroutine xc_xpkzb99(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, tau)


c   
c$Id: xc_xpkzb99.F,v 1.1 2003-02-18 18:40:57 pollack Exp $
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
      double precision amat(nq,ipol), cmat(nq,*)
 

c
      double precision tol_rho, pi, um, uk, umk
      double precision C, Cs
      double precision F43, F13

      parameter(um=0.2195149727645171d0, uk=0.8040d0, umk=um/uk)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)

c
      integer n
      double precision rrho, rho43, rho13, gamma
         

C below are additional declarations needed for this Meta GGA
C  note this  meta-GGA is similar in form to PBE, hence the
C  PBE X backbone is used for this code

c
c     kinetic energy density   or  tau
c
      double precision tau(nq,3,ipol)

      double precision tauN

       double precision  p, qtil, x
       double precision   Cp, rho53, rho83,  mt
       double precision   F83, F23, F53
       double precision   G920, G1081, G146, G73
       double precision   bigD

        parameter (F83=8.d0/3.d0, F23=2.d0/3.d0, F53=5.d0/3.d0)
        parameter (G920=9.d0/20.d0, G1081=10.d0/81.d0,
     &             G146=146.d0/2025.d0,G73=73.d0/405.d0 )
        parameter (bigD=0.113d0)


         
c
      pi = acos(-1.d0)
      C = -3d0/(4d0*pi)*(3d0*pi*pi)**F13
      Cs = 0.5d0/(3d0*pi*pi)**F13
      Cp=Cs*Cs
      Cs = Cs * C               ! account for including C in rho43




c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c                     or
c                SPIN-UNPOLARIZED
c
c

         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10

c  rho43= e_x^unif=exchange energy per electron for uniform electron gas
c       = Ax*n^(1/3)   or C*n^(1/3)

            rho43 = C*rho(n,1)**F43
            rrho = 1d0/rho(n,1)
            rho13 = F43*rho43*rrho

            rho53 = rho(n,1)**F53
            rho83 = rho(n,1)**F83


C  Below we just sum up the LDA contribution to the functional
            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*fac              
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
c            gam12 = dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 10


            tauN = tau(n,1,1) + 
     &             tau(n,2,1) + 
     &             tau(n,3,1)

            p=Cp*gamma/rho83
            qtil=(3.d0*tauN*Cp*2.0/rho53)-G920-(p/12.d0)


c
c     Evaluate the GC part of Fx, i.e. mt = Fx(p,qtil) - 1
c


             x= G1081*p + G146*qtil*qtil
            x= x - G73*qtil*p
            x= x + (bigD + G1081*G1081/uk)*p*p

             if (.not.(nlfac.and.x.gt.tol_rho)) goto 10

            
            mt = uk - uk/(1.d0 + x/uk)

  
C    Below we add the MetaGGA correction to the LDA part from above



             Ex = Ex + rho43*mt*qwght(n)*fac
            if(ldew)func(n) = func(n) + rho43*mt*fac


 10      continue
c
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
            rho43 = C*(2d0*rho(n,2))**F43
            rrho = 0.5d0/rho(n,2)
            rho13 = F43*rho43*rrho


            rho53 = (2.d0*rho(n,2))**F53
            rho83 = (2.d0*rho(n,2))**F83

            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
               if(ldew)func(n) = func(n) + rho43*fac*0.5d0
            endif
c



            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)

c            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 25

c   below note factor of two for spin scaling
            tauN = 2.d0*( tau(n,1,1) + 
     &                    tau(n,2,1) + 
     &                    tau(n,3,1)   )

            
c
c     Evaluate the GC part of Fx, i.e. mt(x) = Fx - 1
c           

            p=Cp*(4.d0*gamma)/rho83
            qtil=(3.d0*tauN*Cp*2.0/rho53) - G920 - (p/12.d0)
 
            x= G1081*p + G146*qtil*qtil
            x= x - G73*qtil*p
            x= x + (bigD + G1081*G1081/uk)*p*p


            if (.not.(nlfac.and.x.gt.tol_rho)) goto 25

            mt = uk - uk/(1.d0 + x/uk)

            Ex = Ex + rho43*mt*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*mt*fac*0.5d0
           
c
c     Beta               BETA           BETA
c
 25         continue

            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = C*(2d0*rho(n,3))**F43
            rrho = 0.5d0/rho(n,3)
            rho13 = F43*rho43*rrho

            rho53 = (2.d0*rho(n,3))**F53
            rho83 = (2.d0*rho(n,3))**F83


            if (lfac) then
               Ex = Ex + rho43*qwght(n)*fac*0.5d0
              if(ldew)func(n) = func(n) + rho43*fac*0.5d0
            endif
c

            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
c            gam12 = 2d0*dsqrt(gamma)
c            if (.not.(nlfac.and.gam12.gt.tol_rho)) goto 20

c   below note factor of two for spin scaling
            tauN = 2.d0*( tau(n,1,2) + 
     &                    tau(n,2,2) + 
     &                    tau(n,3,2)   )

            
c
c     Evaluate the GC part of F(x), i.e. mt(x) = Fx - 1
c
         

            p=Cp*(4.d0*gamma)/rho83
            qtil=(3.d0*tauN*Cp*2.0/rho53)-G920-(p/12.d0)


            x= G1081*p + G146*qtil*qtil
            x= x - G73*qtil*p
            x= x + (bigD + G1081*G1081/uk)*p*p

            if (.not.(nlfac.and.x.gt.tol_rho)) goto 20
   
            mt = uk - uk/(1.d0 + x/uk)
            
            Ex = Ex + rho43*mt*qwght(n)*fac*0.5d0
            if(ldew)func(n) = func(n) + rho43*mt*fac*0.5d0

c
20      continue
      endif
c
      return
      end




      Subroutine xc_xpkzb99_d2()
      call errquit(' not coded ',0)
      return
      end


