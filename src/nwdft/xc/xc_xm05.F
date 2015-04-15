c   M05 or M05-2X  exchange functional  
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                              tauU - uniform-gas KE density
c                              ijzy - 1  M05
c                              ijzy - 2  M05-2X  
c     References: 
c     [a]	Zhao, Y.; Schultz, N. E.; Truhlar, D. G. J. Chem. Phys. 2005, 123, 161103;
c      Note that in this communication we interchanged cCab,i and cCss,i in Table 1. 
c     [b]       Zhao, Y.; Schultz, N. E.; Truhlar, D. G. J. Chem. Theory Comput. 2006, in press.    


      Subroutine xc_xm05(tol_rho, fac,lfac,nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, tau, Mmat, ijzy)

c   
c$Id$
c
      implicit none
c
#include "errquit.fh"
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
      integer n, ijzy
      double precision at1, at2, at3, at4, at5, at6, at7, at8, at9
      double precision at, at10, at11, C1, C2, fL, fNL
      double precision rrho, rho43, rho13, rhoo, rho53
      double precision Gamma2, Gamma
      double precision TauUEG, Tsig, Wsig, W1, W2, W3, W4, W5, W6
      double precision W7, W8, W9, W10, W11, Fsig
c
c     kinetic energy density   or  tau
c
      double precision tau(nq,ipol), Mmat(nq,*)

      double precision tauN,tauu,DTol
      double precision F83, F23, F53, F43, F13, F1o3
      double precision F1o4, F2o3, F3o2, F4o3, F4o9, F3o5
      double precision One, Two, Three, Four, Five, Six, Seven, Eight
      double precision Nine, F10, F11
      double precision Cs, Ax, P32, s, s2, En, Ed, E, dE, dEn, dEd

c      functional derivatives below FFFFFFFFFFFF

       double precision dFdW, dWdT, dTdR, dTdTau, dGGAdR, dFdR
       double precision dFdTau, dGGAdG,tauW

c      functional derivatives above FFFFFFFFFFFF


       parameter( pi = 3.1415926535897932384626433832795d0 )
         
        
        parameter (F1o3=1.d0/3.d0, F1o4=1.d0/4.d0, F2o3=2.d0/3.d0, 
     &             F3o2=3.d0/2.d0,F13=1.d0/3.d0)
        parameter (F4o3=4.d0/3.d0, F4o9=4.d0/9.d0, F3o5=3.d0/5.d0)
        parameter (F83=8.d0/3.0d0, F23=2.0d0/3.d0, F53=5.d0/3.d0)
        parameter (One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0, 
     &             Five=5.0d0,Six=6.0d0, Seven=7.0d0,
     &             Eight=8.0d0, Nine=9.0d0,F10=10.d0, F11=11.d0)


      if (ijzy.eq.1) then
c     Parameters for M05 
        at1=    0.08151d0
        at2=    -0.43956d0
        at3=    -3.22422d0
        at4=    2.01819d0
        at5=    8.79431d0
        at6=    -0.00295d0
        at7=    9.82029d0
        at8=    -4.82351d0
        at9=    -48.17574d0
        at10=   3.64802d0
        at11=   34.02248d0
      elseif (ijzy.eq.2) then
c     Parameters for M05-2X
        at1=    -0.56833d0
        at2=    -1.30057d0
        at3=    5.50070d0
        at4=    9.06402d0
        at5=    -32.21075d0
        at6=    -23.73298d0
        at7=    70.22996d0
        at8=    29.88614d0
        at9=    -60.25778d0
        at10=   -13.22205d0
        at11=   15.23694d0
      else
        call errquit("xc_xm05: illegal ijzy",ijzy,UERR)
      endif
      
      at=1.0d0
      C1 = 0.2195149727645171d0
      C2 = C1/0.804d0 
      DTol=tol_rho
C
C     Scale factors for local and non-local contributions.
C
      fL  =  fac
      fNL =  fac
      Cs = 0.5d0/(3.0d0*pi*pi)**F13
      P32 = (3.d0*pi**2)**F23
         
c     
       Ax = (-0.75d0)*(3.0d0/pi)**F13


c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c                     or
c                SPIN-UNPOLARIZED
c
c
         do 10 n = 1, nq
           
            if (rho(n,1).lt.DTol) goto 10
            rhoo = rho(n,1)
            rho43 = rhoo**F4o3  
            rrho = 1d0/rhoo       ! reciprocal of rho
            rho13 = rho43*rrho
            rho53 = rhoo**F53

c
             
            tauN = tau(n,1)
            tauu=tauN 
            Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 10
            TauUEG=0.3d0*P32*rho53
            Tsig =TauUEG/tauu
            Wsig =(Tsig-One)/(Tsig+One)
            W1=Wsig 
            W2=Wsig*W1
            W3=Wsig*W2
            W4=Wsig*W3
            W5=Wsig*W4
            W6=Wsig*W5
            W7=Wsig*W6
            W8=Wsig*W7
            W9=Wsig*W8
            W10=Wsig*W9
            W11=Wsig*W10
            Fsig =at*(at1*W1+ at2*W2 + at3*W3
     &          + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &          + at8*W8 + at9*W9 + at10*W10 + at11*W11)

            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
            Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*qwght(n)
            if(ldew) func(n)=func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)
c
c     functional derivatives 
c
            dEn   = Two*C1*s
            dEd   = Two*C2*s
            dE    = (dEn*Ed-En*dEd)/(Ed*Ed)
            dFdW = at*(at1 + Two*at2*W1 + Three*at3*W2
     &             + Four*at4*W3 + Five*at5*W4
     &             + Six*at6*W5 + Seven*at7*W6
     &             + Eight*at8*W7 + Nine*at9*W8
     &             + F10*at10*W9 + F11*at11*W10)
            dWdT = Two/((One + Tsig)**2)
            dTdR = (0.5d0*P32*rho13*rho13)/tauu
            dTdTau = -TauUEG/tauu**2
            dGGAdR = F4o3*rho13*Ax*(fL+fNL*(E-s*dE))
            dFdR = dFdW*dWdT*dTdR
            dFdTau=dFdW*dWdT*dTdTau
            dGGAdG =(fNL*dE*s/(Two*Gamma2))
            Amat(n,1) = Amat(n,1) + dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR
            Cmat(n,1)=  Cmat(n,1) + 
     &                    Two*dGGAdG*Ax*rho43*(One+Fsig) 
            Mmat(n,1)=  Mmat(n,1) + 0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau



10      continue


c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======

c
c  use spin density functional theory ie n-->2n
c                               Ex=(1/2)Ex[2*alpha] + (1/2)Ex[2*beta]

         do 20 n = 1, nq
           if (rho(n,1).lt.DTol) goto 20
c
c     Alpha            ALPHA               ALPHA
c
            if (rho(n,2).lt.DTol) goto 25
             rhoo = Two*rho(n,2)
             rho43 = rhoo**F4o3
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**F53

c

             tauN = tau(n,1)
             tauu = Two*tauN
             TauUEG=0.3d0*P32*rho53
             Tsig =TauUEG/tauu
             Wsig =(Tsig-One)/(Tsig+One)
             W1=Wsig
             W2=Wsig*W1
             W3=Wsig*W2
             W4=Wsig*W3
             W5=Wsig*W4
             W6=Wsig*W5
             W7=Wsig*W6
             W8=Wsig*W7
             W9=Wsig*W8
             W10=Wsig*W9
             W11=Wsig*W10
             Fsig =at*(at1*W1+ at2*W2 + at3*W3
     &           + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &           + at8*W8 + at9*W9 + at10*W10 + at11*W11)


            Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 25

            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
            Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*0.5d0*qwght(n)
            if(ldew) func(n)=
     =           func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)*.5d0
c
c     functional derivatives
c
            dEn   = Two*C1*s
            dEd   = Two*C2*s
            dE    = (dEn*Ed-En*dEd)/(Ed*Ed)
            dFdW = at*(at1 + Two*at2*W1 + Three*at3*W2
     &             + Four*at4*W3 + Five*at5*W4
     &             + Six*at6*W5 + Seven*at7*W6
     &             + Eight*at8*W7 + Nine*at9*W8
     &             + F10*at10*W9 + F11*at11*W10)
            dWdT = Two/((One + Tsig)**2)
            dTdR = Two*(0.5d0*P32*rho13*rho13)/tauu
            dTdTau = -Two*TauUEG/tauu**2
            dGGAdR = Two*F4o3*rho13*Ax*(fL+fNL*(E-s*dE))
            dFdR = dFdW*dWdT*dTdR
            dFdTau=dFdW*dWdT*dTdTau
            dGGAdG =Four*(fNL*dE*s/(Two*Gamma2))

            Amat(n,1) = Amat(n,1) + (dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,1)=  Cmat(n,1) + 
     &                      dGGAdG*Ax*rho43*(One+Fsig)*0.5d0
            Mmat(n,1)=  Mmat(n,1) + 
     &                  0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau

c             write (*,*) "Ex,Amat(n,1),Cmat(n,1)",
c     &        Ex,Amat(n,1),Cmat(n,1)

c
c     Beta               BETA           BETA
c

25         continue

c
c     Beta
c
            if (rho(n,3).lt.DTol) goto 20
             rhoo = Two*rho(n,3)
             rho43 = rhoo**F4o3
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**F53

c

             tauN = tau(n,2)
             tauu = Two*tauN
             TauUEG=0.3d0*P32*rho53
             Tsig =TauUEG/tauu
             Wsig =(Tsig-One)/(Tsig+One)
             W1=Wsig
             W2=Wsig*W1
             W3=Wsig*W2
             W4=Wsig*W3
             W5=Wsig*W4
             W6=Wsig*W5
             W7=Wsig*W6
             W8=Wsig*W7
             W9=Wsig*W8
             W10=Wsig*W9
             W11=Wsig*W10
             Fsig =at*(at1*W1+ at2*W2 + at3*W3
     &           + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &           + at8*W8 + at9*W9 + at10*W10 + at11*W11)


            Gamma2 = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 20
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
            Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*0.5d0*qwght(n)
            if(ldew) func(n)=
     =           func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)*.5d0
c
c     functional derivatives
c
            dEn   = Two*C1*s
            dEd   = Two*C2*s
            dE    = (dEn*Ed-En*dEd)/(Ed*Ed)
            dFdW = at*(at1 + Two*at2*W1 + Three*at3*W2
     &             + Four*at4*W3 + Five*at5*W4
     &             + Six*at6*W5 + Seven*at7*W6
     &             + Eight*at8*W7 + Nine*at9*W8
     &             + F10*at10*W9 + F11*at11*W10)
            dWdT = Two/((One + Tsig)**2)
            dTdR = Two*(0.5d0*P32*rho13*rho13)/tauu
            dTdTau = -Two*TauUEG/tauu**2
            dGGAdR = Two*F4o3*rho13*Ax*(fL+fNL*(E-s*dE))
            dFdR = dFdW*dWdT*dTdR
            dFdTau=dFdW*dWdT*dTdTau
            dGGAdG =Four*(fNL*dE*s/(Two*Gamma2))

            Amat(n,2) = Amat(n,2) + (dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,3)=  Cmat(n,3) + 
     &                   dGGAdG*Ax*rho43*(One+Fsig)*0.5d0
            Mmat(n,2)=  Mmat(n,2) +
     &                  0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau
     

c
20      continue
      endif
c
      return
      end




      Subroutine xc_xm05_d2()
      call errquit(' xm05: d2 not coded ',0,0)
      return
      end


c----------------------------------------------------------------------
c   dlDF  exchange functional  
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                              tauU - uniform-gas KE density
c                              
c     References: 
c     [a]	Pernal,Podeszwa,Patkowski,Szalewicz, PRL 103 263201 (2009)

      Subroutine xc_xdldf(tol_rho, fac,lfac,nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, tau, Mmat)
c   
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
      double precision at1, at2, at3, at4, at5, at6, at7, at8, at9
      double precision at, at10, at11, C1, C2, fL, fNL
      double precision rrho, rho43, rho13, rhoo, rho53
      double precision Gamma2, Gamma
      double precision TauUEG, Tsig, Wsig, W1, W2, W3, W4, W5, W6
      double precision W7, W8, W9, W10, W11, Fsig
c
c     kinetic energy density   or  tau
c
      double precision tau(nq,ipol), Mmat(nq,*)

      double precision tauN,tauu,DTol
      double precision F83, F23, F53, F43, F13, F1o3
      double precision F1o4, F2o3, F3o2, F4o3, F4o9, F3o5
      double precision One, Two, Three, Four, Five, Six, Seven, Eight
      double precision Nine, F10, F11
      double precision Cs, Ax, P32, s, s2, En, Ed, E, dE, dEn, dEd

c      functional derivatives below FFFFFFFFFFFF

       double precision dFdW, dWdT, dTdR, dTdTau, dGGAdR, dFdR
       double precision dFdTau, dGGAdG,tauW

c      functional derivatives above FFFFFFFFFFFF


       parameter( pi = 3.1415926535897932384626433832795d0 )
         
        
        parameter (F1o3=1.d0/3.d0, F1o4=1.d0/4.d0, F2o3=2.d0/3.d0, 
     &             F3o2=3.d0/2.d0,F13=1.d0/3.d0)
        parameter (F4o3=4.d0/3.d0, F4o9=4.d0/9.d0, F3o5=3.d0/5.d0)
        parameter (F83=8.d0/3.0d0, F23=2.0d0/3.d0, F53=5.d0/3.d0)
        parameter (One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0, 
     &             Five=5.0d0,Six=6.0d0, Seven=7.0d0,
     &             Eight=8.0d0, Nine=9.0d0,F10=10.d0, F11=11.d0)
   
c     Parameters for dlDF
        at1=    -1.637571d-01
        at2=    -1.880028d-01
        at3=    -4.490609d-01
        at4=    -8.2359d-03
        at5=     0.0d0
        at6=     0.0d0
        at7=     0.0d0
        at8=     0.0d0
        at9=     0.0d0
        at10=    0.0d0
        at11=    0.0d0
      
      
      at=1.0d0
      C1 = 0.3511128d0
      C2 = C1/4.8827323d0 
      DTol=tol_rho
C
C     Scale factors for local and non-local contributions.
C
      fL  =  fac
      fNL =  fac
      Cs = 0.5d0/(3.0d0*pi*pi)**F13
      P32 = (3.d0*pi**2)**F23
         
c     
       Ax = (-0.75d0)*(3.0d0/pi)**F13


c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c                     or
c                SPIN-UNPOLARIZED
c
c
         do 10 n = 1, nq
           
            if (rho(n,1).lt.DTol) goto 10
            rhoo = rho(n,1)
            rho43 = rhoo**F4o3  
            rrho = 1d0/rhoo       ! reciprocal of rho
            rho13 = rho43*rrho
            rho53 = rhoo**F53

c
             
            tauN = tau(n,1)
            tauu=tauN 
cedo            if (taun.lt.sqrt(DTol)) goto 10
            Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 10
            TauUEG=0.3d0*P32*rho53
            Tsig =TauUEG/tauu
            Wsig =(Tsig-One)/(Tsig+One)
            W1=Wsig 
            W2=Wsig*W1
            W3=Wsig*W2
            W4=Wsig*W3
            W5=Wsig*W4
            W6=Wsig*W5
            W7=Wsig*W6
            W8=Wsig*W7
            W9=Wsig*W8
            W10=Wsig*W9
            W11=Wsig*W10
            Fsig =at*(at1*W1+ at2*W2 + at3*W3
     &          + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &          + at8*W8 + at9*W9 + at10*W10 + at11*W11)

            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
            Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*qwght(n)
            if(ldew) func(n)=func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)
c
c     functional derivatives 
c
            dEn   = Two*C1*s
            dEd   = Two*C2*s
            dE    = (dEn*Ed-En*dEd)/(Ed*Ed)
            dFdW = at*(at1 + Two*at2*W1 + Three*at3*W2
     &             + Four*at4*W3 + Five*at5*W4
     &             + Six*at6*W5 + Seven*at7*W6
     &             + Eight*at8*W7 + Nine*at9*W8
     &             + F10*at10*W9 + F11*at11*W10)
            dWdT = Two/((One + Tsig)**2)
            dTdR = (0.5d0*P32*rho13*rho13)/tauu
            dTdTau = -TauUEG/tauu**2
            dGGAdR = F4o3*rho13*Ax*(fL+fNL*(E-s*dE))
            dFdR = dFdW*dWdT*dTdR
            dFdTau=dFdW*dWdT*dTdTau
            dGGAdG =(fNL*dE*s/(Two*Gamma2))
            Amat(n,1) = Amat(n,1) + dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR
            Cmat(n,1)=  Cmat(n,1) + 
     &                    Two*dGGAdG*Ax*rho43*(One+Fsig) 
            Mmat(n,1)=  Mmat(n,1) + 0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau



10      continue


c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======

c
c  use spin density functional theory ie n-->2n
c                               Ex=(1/2)Ex[2*alpha] + (1/2)Ex[2*beta]

         do 20 n = 1, nq
           if (rho(n,1).lt.DTol) goto 20
c
c     Alpha            ALPHA               ALPHA
c
            if (rho(n,2).lt.DTol) goto 25
             rhoo = Two*rho(n,2)
             rho43 = rhoo**F4o3
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**F53

c

             tauN = tau(n,1)
             tauu = Two*tauN
             TauUEG=0.3d0*P32*rho53
             Tsig =TauUEG/tauu
             Wsig =(Tsig-One)/(Tsig+One)
             W1=Wsig
             W2=Wsig*W1
             W3=Wsig*W2
             W4=Wsig*W3
             W5=Wsig*W4
             W6=Wsig*W5
             W7=Wsig*W6
             W8=Wsig*W7
             W9=Wsig*W8
             W10=Wsig*W9
             W11=Wsig*W10
             Fsig =at*(at1*W1+ at2*W2 + at3*W3
     &           + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &           + at8*W8 + at9*W9 + at10*W10 + at11*W11)


            Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 25

            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
            Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*0.5d0*qwght(n)
            if(ldew) func(n)=
     =           func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)*.5d0
c
c     functional derivatives
c
            dEn   = Two*C1*s
            dEd   = Two*C2*s
            dE    = (dEn*Ed-En*dEd)/(Ed*Ed)
            dFdW = at*(at1 + Two*at2*W1 + Three*at3*W2
     &             + Four*at4*W3 + Five*at5*W4
     &             + Six*at6*W5 + Seven*at7*W6
     &             + Eight*at8*W7 + Nine*at9*W8
     &             + F10*at10*W9 + F11*at11*W10)
            dWdT = Two/((One + Tsig)**2)
            dTdR = Two*(0.5d0*P32*rho13*rho13)/tauu
            dTdTau = -Two*TauUEG/tauu**2
            dGGAdR = Two*F4o3*rho13*Ax*(fL+fNL*(E-s*dE))
            dFdR = dFdW*dWdT*dTdR
            dFdTau=dFdW*dWdT*dTdTau
            dGGAdG =Four*(fNL*dE*s/(Two*Gamma2))

            Amat(n,1) = Amat(n,1) + (dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,1)=  Cmat(n,1) + 
     &                      dGGAdG*Ax*rho43*(One+Fsig)*0.5d0
            Mmat(n,1)=  Mmat(n,1) + 
     &                  0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau

c             write (*,*) "Ex,Amat(n,1),Cmat(n,1)",
c     &        Ex,Amat(n,1),Cmat(n,1)

c
c     Beta               BETA           BETA
c

25         continue

c
c     Beta
c
            if (rho(n,3).lt.DTol) goto 20
             rhoo = Two*rho(n,3)
             rho43 = rhoo**F4o3
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**F53

c

             tauN = tau(n,2)
             tauu = Two*tauN
             TauUEG=0.3d0*P32*rho53
             Tsig =TauUEG/tauu
             Wsig =(Tsig-One)/(Tsig+One)
             W1=Wsig
             W2=Wsig*W1
             W3=Wsig*W2
             W4=Wsig*W3
             W5=Wsig*W4
             W6=Wsig*W5
             W7=Wsig*W6
             W8=Wsig*W7
             W9=Wsig*W8
             W10=Wsig*W9
             W11=Wsig*W10
             Fsig =at*(at1*W1+ at2*W2 + at3*W3
     &           + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &           + at8*W8 + at9*W9 + at10*W10 + at11*W11)


            Gamma2 = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 20
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
            Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*0.5d0*qwght(n)
            if(ldew) func(n)=
     =           func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)*.5d0
c
c     functional derivatives
c
            dEn   = Two*C1*s
            dEd   = Two*C2*s
            dE    = (dEn*Ed-En*dEd)/(Ed*Ed)
            dFdW = at*(at1 + Two*at2*W1 + Three*at3*W2
     &             + Four*at4*W3 + Five*at5*W4
     &             + Six*at6*W5 + Seven*at7*W6
     &             + Eight*at8*W7 + Nine*at9*W8
     &             + F10*at10*W9 + F11*at11*W10)
            dWdT = Two/((One + Tsig)**2)
            dTdR = Two*(0.5d0*P32*rho13*rho13)/tauu
            dTdTau = -Two*TauUEG/tauu**2
            dGGAdR = Two*F4o3*rho13*Ax*(fL+fNL*(E-s*dE))
            dFdR = dFdW*dWdT*dTdR
            dFdTau=dFdW*dWdT*dTdTau
            dGGAdG =Four*(fNL*dE*s/(Two*Gamma2))

            Amat(n,2) = Amat(n,2) + (dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,3)=  Cmat(n,3) + 
     &                   dGGAdG*Ax*rho43*(One+Fsig)*0.5d0
            Mmat(n,2)=  Mmat(n,2) +
     &                  0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau
     

c
20      continue
      endif
c
      return
      end



      Subroutine xc_xdldf_d2()
      call errquit('xdldf: d2 not coded',0,0)
      return
      end


