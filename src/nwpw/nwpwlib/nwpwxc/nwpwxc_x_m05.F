C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_m05.F
C> Implementation of the M05 correlation functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The M05 correlation functional
C>
C> The M05 functional [1,2] is a meta-GGA of which this evaluates
C> the correlation component.
C>
C> ### References ###
C>
C> [1] Y Zhao, NE Schultz, DG Truhlar,
C>     "Exchange-correlation functional with broad accuracy for 
C>     metallic and nonmetallic compounds, kinetics, and 
C>     noncovalent interactions",
C>     J.Chem.Phys. <b>123</b>, 161103-161106 (2005), DOI:
C>     <a href="https://doi.org/10.1063/1.2126975">
C>     10.1063/1.2126975</a>.
C>
C> [2] Y Zhao, NE Schultz, DG Truhlar,
C>     "Design of density functionals by combining the method of
C>     constraint satisfaction parametrization for thermochemistry,
C>     thermochemical kinetics, and noncovalent interactions",
C>     J.Chem.Theory.Comput. <b>2</b>, 364-382 (2006), DOI:
C>     <a href="https://doi.org/10.1021/ct0502763">
C>     10.1021/ct0502763</a>.
C>

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


      Subroutine nwpwxc_x_m05(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      tau, func, Amat, Cmat, Mmat)
c   
c$Id$
c
      implicit none
c
#include "nwpwxc_param.fh"
c
c
c     Input and other parameters
c
      double precision param(*) !< [Input] Parameters of functional
                                !< - param(1): \f$ a_{1} \f$
                                !< - param(2): \f$ a_{2} \f$
                                !< - param(3): \f$ a_{3} \f$
                                !< - param(4): \f$ a_{4} \f$
                                !< - param(5): \f$ a_{5} \f$
                                !< - param(6): \f$ a_{6} \f$
                                !< - param(7): \f$ a_{7} \f$
                                !< - param(8): \f$ a_{8} \f$
                                !< - param(9): \f$ a_{9} \f$
                                !< - param(10): \f$ a_{10} \f$
                                !< - param(11): \f$ a_{11} \f$
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient Norm
c
      double precision rgamma(nq,*) !< [Input] The density gradient norm
c
c     Kinetic Energy Density
c
      double precision tau(nq,*) !< [Input] The kinetic energy density
c
c     Functional values
c
      double precision func(*) !< [Output] The functional value
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,*) !< [Output] Derivative wrt density
      double precision Cmat(nq,*) !< [Output] Derivative wrt rgamma
      double precision Mmat(nq,*) !< [Output] Derivative wrt tau
c      
      double precision pi
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


      at1=  param(1)
      at2=  param(2)
      at3=  param(3)
      at4=  param(4)
      at5=  param(5)
      at6=  param(6)
      at7=  param(7)
      at8=  param(8)
      at9=  param(9)
      at10= param(10)
      at11= param(11)
c     if (ijzy.eq.1) then
c     Parameters for M05 
c       at1=    0.08151d0
c       at2=    -0.43956d0
c       at3=    -3.22422d0
c       at4=    2.01819d0
c       at5=    8.79431d0
c       at6=    -0.00295d0
c       at7=    9.82029d0
c       at8=    -4.82351d0
c       at9=    -48.17574d0
c       at10=   3.64802d0
c       at11=   34.02248d0
c     elseif (ijzy.eq.2) then
c     Parameters for M05-2X
c       at1=    -0.56833d0
c       at2=    -1.30057d0
c       at3=    5.50070d0
c       at4=    9.06402d0
c       at5=    -32.21075d0
c       at6=    -23.73298d0
c       at7=    70.22996d0
c       at8=    29.88614d0
c       at9=    -60.25778d0
c       at10=   -13.22205d0
c       at11=   15.23694d0
c     endif
      
      at=1.0d0
      C1 = 0.2195149727645171d0
      C2 = C1/0.804d0 
      DTol=tol_rho
C
C     Scale factors for local and non-local contributions.
C
      fL  =  wght
      fNL =  wght
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
           
            if (rho(n,R_T).lt.DTol) goto 10
            rhoo = rho(n,R_T)
            rho43 = rhoo**F4o3  
            rrho = 1d0/rhoo       ! reciprocal of rho
            rho13 = rho43*rrho
            rho53 = rhoo**F53

c
             
            tauN = tau(n,T_T)
            tauu=tauN 
c           Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = rgamma(n,G_TT)
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
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*qwght(n)
            func(n)=func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)
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
            Amat(n,D1_RA) = Amat(n,D1_RA) + dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR
            Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + 
     &                    Two*dGGAdG*Ax*rho43*(One+Fsig) 
            Mmat(n,D1_TA)=  Mmat(n,D1_TA)
     &                   + 0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau



10      continue


c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======

c
c  use spin density functional theory ie n-->2n
c                               Ex=(1/2)Ex[2*alpha] + (1/2)Ex[2*beta]

         do 20 n = 1, nq
           if (rho(n,R_A)+rho(n,R_B).lt.DTol) goto 20
c
c     Alpha            ALPHA               ALPHA
c
            if (rho(n,R_A).lt.0.5d0*DTol) goto 25
             rhoo = Two*rho(n,R_A)
             rho43 = rhoo**F4o3
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**F53

c

             tauN = tau(n,T_A)
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


c           Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = rgamma(n,G_AA)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 25

            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*0.5d0*qwght(n)
            func(n)= func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)*.5d0
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

            Amat(n,D1_RA) = Amat(n,D1_RA) + (dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + 
     &                      dGGAdG*Ax*rho43*(One+Fsig)*0.5d0
            Mmat(n,D1_TA)=  Mmat(n,D1_TA) + 
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
            if (rho(n,R_B).lt.0.5d0*DTol) goto 20
             rhoo = Two*rho(n,R_B)
             rho43 = rhoo**F4o3
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**F53

c

             tauN = tau(n,T_B)
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


c           Gamma2 = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            Gamma2 = rgamma(n,G_BB)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 20
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*0.5d0*qwght(n)
            func(n)= func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)*.5d0
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

            Amat(n,D1_RB) = Amat(n,D1_RB) + (dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,D1_GBB)=  Cmat(n,D1_GBB) + 
     &                   dGGAdG*Ax*rho43*(One+Fsig)*0.5d0
            Mmat(n,D1_TB)=  Mmat(n,D1_TB) +
     &                  0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau
     

c
20      continue
      endif
c
      return
      end




      Subroutine nwpwxc_x_m05_d2()
      call errquit(' xm05: d2 not coded ',0,0)
      return
      end


c----------------------------------------------------------------------
C> \brief Calculate the dlDF correlation functional
C>
C> Calculate the dlDF correlation functional [1].
C>
C> ### References ###
C>
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz,
C> "Dispersionless density functional theory",
C> Phys.Rev.Lett. <b>103</b>, 263201-263204 (2009), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.103.263201">
C> 10.1103/PhysRevLett.103.263201</a>.
C>
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

      Subroutine nwpwxc_x_dldf(tol_rho,ipol,nq,wght,rho,rgamma,tau,
     &                       func, Amat, Cmat, Mmat)
c
      implicit none
c
#include "nwpwxc_param.fh"
c      
c
c     Input and other parameters
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient Norm
c
      double precision rgamma(nq,*) !< [Input] The density gradient norm
c
c     Kinetic Energy Density
c
      double precision tau(nq,*) !< [Input] The kinetic energy density
c
c     Functional values
c
      double precision func(*) !< [Output] The functional value
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,*) !< [Output] Derivative wrt density
      double precision Cmat(nq,*) !< [Output] Derivative wrt rgamma
      double precision Mmat(nq,*) !< [Output] Derivative wrt tau
c
      double precision pi
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
      fL  =  wght
      fNL =  wght
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
           
            if (rho(n,R_T).lt.DTol) goto 10
            rhoo = rho(n,R_T)
            rho43 = rhoo**F4o3  
            rrho = 1d0/rhoo       ! reciprocal of rho
            rho13 = rho43*rrho
            rho53 = rhoo**F53

c
             
            tauN = tau(n,T_T)
            tauu=tauN 
cedo            if (taun.lt.sqrt(DTol)) goto 10
c           Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = rgamma(n,G_TT)
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
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*qwght(n)
            func(n)=func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)
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
            Amat(n,D1_RA) = Amat(n,D1_RA) + dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR
            Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + 
     &                    Two*dGGAdG*Ax*rho43*(One+Fsig) 
            Mmat(n,D1_TA)=  Mmat(n,D1_TA)
     &                   + 0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau



10      continue


c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======

c
c  use spin density functional theory ie n-->2n
c                               Ex=(1/2)Ex[2*alpha] + (1/2)Ex[2*beta]

         do 20 n = 1, nq
           if (rho(n,R_A)+rho(n,R_B).lt.DTol) goto 20
c
c     Alpha            ALPHA               ALPHA
c
            if (rho(n,R_A).lt.0.5d0*DTol) goto 25
             rhoo = Two*rho(n,R_A)
             rho43 = rhoo**F4o3
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**F53

c

             tauN = tau(n,T_A)
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


c           Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = rgamma(n,G_AA)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 25

            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*0.5d0*qwght(n)
            func(n)= func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)*.5d0
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

            Amat(n,D1_RA) = Amat(n,D1_RA) + (dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + 
     &                      dGGAdG*Ax*rho43*(One+Fsig)*0.5d0
            Mmat(n,D1_TA)=  Mmat(n,D1_TA) + 
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
            if (rho(n,R_B).lt.0.5d0*DTol) goto 20
             rhoo = Two*rho(n,R_B)
             rho43 = rhoo**F4o3
             rrho = 1.0d0/rhoo       ! reciprocal of rho
             rho13 = rho43*rrho
             rho53 = rhoo**F53

c

             tauN = tau(n,T_B)
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


c           Gamma2 = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            Gamma2 = rgamma(n,G_BB)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if (gamma.lt.DTol) goto 20
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(One+Fsig)*0.5d0*qwght(n)
            func(n)=func(n)+rho43*Ax*(fL+fNL*E)*(One+Fsig)*.5d0
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

            Amat(n,D1_RB) = Amat(n,D1_RB) + (dGGAdR*(One+Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,D1_GBB)=  Cmat(n,D1_GBB) + 
     &                   dGGAdG*Ax*rho43*(One+Fsig)*0.5d0
            Mmat(n,D1_TB)=  Mmat(n,D1_TB) +
     &                  0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau
     

c
20      continue
      endif
c
      return
      end



      Subroutine nwpwxc_x_dldf_d2()
      call errquit('xdldf: d2 not coded',0,0)
      return
      end
C> @}


