C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_m06.F
C> Implementation of the M06 exchange functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The M06 exchange functional
C>
C> The M06 functional [1,2] is a meta-GGA of which this evaluates
C> the exchange component.
C>
C> ### References ###
C>
C> [1] Y Zhao, DG Truhlar,
C> "A new local density functional for main-group thermochemistry,
C> transition metal bonding, thermochemical kinetics, and noncovalent
C> interactions",
C> J. Chem. Phys. <b>125</b>, 194101 (2006), DOI:
C> <a href="https://doi.org/10.1063/1.2370993">
C> 10.1063/1.2370993</a>.
C>
C> [2] Y Zhao, DG Truhlar,
C> "Density functional for spectroscopy: No long-range self-interaction
C> error, good performance for Rydberg and charge-transfer states,
C> and better performance on average than B3LYP for ground states",
C> J. Phys. Chem. A <b>110</b>, 13126-13130 (2006), DOI:
C> <a href="https://doi.org/10.1021/jp066479k">
C> 10.1021/jp066479k</a>.
C>
c   M06 suite  exchange functional  
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                              tauU - uniform-gas KE density
c                              ijzy - 1  M06-L  
c                              ijzy - 2  M06-HF
c                              ijzy - 3  M06
c                              ijzy - 4  M06-2X
c     References: 
c     [a]	Zhao, Y. and  Truhlar, D. G. J. Chem. Phys. 2006, 125, 194101;
c     [b]       Zhao, Y. and  Truhlar, D. G. J. Phys. Chem. A (2006),110(49),13126-13130.    


      Subroutine nwpwxc_x_m06(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      tau, func, Amat, Cmat, Mmat)
c   
c$Id$
c
      implicit none
c
#include "nwpwxc_param.fh"
c
      double precision param(*) !< [Input] Parameters of the functional
                                !< - param( 1): \f$ d_0 \f$
                                !< - param( 2): \f$ d_1 \f$
                                !< - param( 3): \f$ d_2 \f$
                                !< - param( 4): \f$ d_3 \f$
                                !< - param( 5): \f$ d_4 \f$
                                !< - param( 6): \f$ d_5 \f$
                                !< - param( 7): \f$ a_0 \f$
                                !< - param( 8): \f$ a_1 \f$
                                !< - param( 9): \f$ a_2 \f$
                                !< - param(10): \f$ a_3 \f$
                                !< - param(11): \f$ a_4 \f$
                                !< - param(12): \f$ a_5 \f$
                                !< - param(13): \f$ a_6 \f$
                                !< - param(14): \f$ a_7 \f$
                                !< - param(15): \f$ a_8 \f$
                                !< - param(16): \f$ a_9 \f$
                                !< - param(17): \f$ a_10 \f$
                                !< - param(18): \f$ a_11 \f$
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
      double precision at, at10, at11, at0, C1, C2, fL, fNL
      double precision rrho, rho43, rho13, rhoo, rho53
      double precision Gamma2, Gamma
      double precision TauUEG, Tsig, Wsig, W1, W2, W3, W4, W5, W6
      double precision W7, W8, W9, W10, W11, Fsig
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

c     functional derivatives above FFFFFFFFFFFF


cedo       parameter( pi = 3.1415926535897932384626433832795d0 )
         
      
      parameter (F1o3=1.d0/3.d0, F1o4=1.d0/4.d0, F2o3=2.d0/3.d0, 
     &     F3o2=3.d0/2.d0,F13=1.d0/3.d0)
      parameter (F4o3=4.d0/3.d0, F4o9=4.d0/9.d0, F3o5=3.d0/5.d0)
      parameter (F83=8.d0/3.0d0, F23=2.0d0/3.d0, F53=5.d0/3.d0)
      parameter (One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0, 
     &     Five=5.0d0,Six=6.0d0, Seven=7.0d0,
     &     Eight=8.0d0, Nine=9.0d0,F10=10.d0, F11=11.d0)
      pi=acos(-1d0)      

      at0  = param(7)
      at1  = param(8)
      at2  = param(9)
      at3  = param(10)
      at4  = param(11)
      at5  = param(12)
      at6  = param(13)
      at7  = param(14)
      at8  = param(15)
      at9  = param(16)
      at10 = param(17)
      at11 = param(18)
c     if (ijzy.eq.1) then
c       at0=    3.987756D-01
c       at1=    2.548219D-01
c       at2=    3.923994D-01
c       at3=    -2.103655D+00
c       at4=    -6.302147D+00
c       at5=    1.097615D+01
c       at6=    3.097273D+01
c       at7=    -2.318489D+01
c       at8=    -5.673480D+01
c       at9=    2.160364D+01
c       at10=   3.421814D+01
c       at11=   -9.049762D+00
c      elseif (ijzy.eq.2) then
c     Parameters for M06-HF
c       at0=    1.179732D-01
c       at1=    -1.066708D+00
c       at2=    -1.462405D-01
c       at3=    7.481848D+00
c       at4=    3.776679D+00
c       at5=    -4.436118D+01
c       at6=    -1.830962D+01
c       at7=    1.003903D+02
c       at8=    3.864360D+01
c       at9=    -9.806018D+01
c       at10=   -2.557716D+01
c       at11=   3.590404D+01
c      elseif (ijzy.eq.3) then
c     Parameters for M06
c       at0=    5.877943D-01
c       at1=    -1.371776D-01
c       at2=    2.682367D-01
c       at3=    -2.515898D+00
c       at4=    -2.978892D+00
c       at5=    8.710679D+00
c       at6=    1.688195D+01
c       at7=    -4.489724D+00
c       at8=    -3.299983D+01
c       at9=    -1.449050D+01
c       at10=   2.043747D+01
c       at11=   1.256504D+01
c      elseif (ijzy.eq.4) then
c     Parameters for M06-2X
c       at0=    4.600000D-01
c       at1=    -2.206052D-01
c       at2=    -9.431788D-02
c       at3=    2.164494D+00
c       at4=    -2.556466D+00
c       at5=    -1.422133D+01
c       at6=    1.555044D+01
c       at7=    3.598078D+01
c       at8=    -2.722754D+01
c       at9=    -3.924093D+01
c       at10=   1.522808D+01
c       at11=   1.522227D+01
c     endif

      at=1.0d0
      call nwpwxc_x_vs98(param,tol_rho,ipol,nq,wght,rho,rgamma,tau,
     &                 func, Amat, Cmat, Mmat)


      C1 = 0.2195149727645171d0
      C2 = C1/0.804d0 
cedo      DTol=1.0D-8
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
         if(taun.lt.dtol) goto 10
            tauu=tauN 
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
            Fsig =at*(at0 + at1*W1+ at2*W2 + at3*W3
     &          + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &          + at8*W8 + at9*W9 + at10*W10 + at11*W11)

c           Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = rgamma(n,G_TT)
            Gamma = dsqrt(Gamma2)
         if(gamma.lt.dtol) goto 10
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*Fsig*qwght(n)
            func(n)=func(n)+rho43*Ax*(fL+fNL*E)*Fsig
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
            Amat(n,D1_RA)  = Amat(n,D1_RA) + dGGAdR*Fsig
     &                     + (fL+fNL*E)*Ax*rho43*dFdR
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA)
     &                     + Two*dGGAdG*Ax*rho43*Fsig 
            Mmat(n,D1_TA)  = Mmat(n,D1_TA)
     &                     + 0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau

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
             Fsig =at*(at0 + at1*W1+ at2*W2 + at3*W3
     &           + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &           + at8*W8 + at9*W9 + at10*W10 + at11*W11)


c           Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = rgamma(n,G_AA)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)

            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(Fsig)*0.5d0*qwght(n)
            func(n)=func(n)+rho43*Ax*(fL+fNL*E)*(Fsig)*.5d0
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

            Amat(n,D1_RA)  = Amat(n,D1_RA) + (dGGAdR*(Fsig)
     &                     + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA)
     &                     + dGGAdG*Ax*rho43*(Fsig)*0.5d0
            Mmat(n,D1_TA)  = Mmat(n,D1_TA)
     &                     + 0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau

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
             Fsig =at*(at0+at1*W1+ at2*W2 + at3*W3
     &           + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &           + at8*W8 + at9*W9 + at10*W10 + at11*W11)


c           Gamma2 = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            Gamma2 = rgamma(n,G_BB)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(Fsig)*0.5d0*qwght(n)
            func(n)= func(n)+rho43*Ax*(fL+fNL*E)*(Fsig)*.5d0
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

            Amat(n,D1_RB)  = Amat(n,D1_RB) + (dGGAdR*(Fsig)
     &                     + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB)
     &                     + dGGAdG*Ax*rho43*(Fsig)*0.5d0
            Mmat(n,D1_TB)  = Mmat(n,D1_TB)
     &                     + 0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau
c
20      continue
      endif
c
      return
      end
C>
C> \brief Evaluate the M06-2X exchange functional
C>
C> This routine evaluates the M06-2X exchange functional [1,2]. This functional
C> is closely related to the M06, M06-HF and M06-L exchange functionals
C> except that it does not use the VS98 exchange functional, hence the
C> parameter list is defined differently.
C>
C> ### References ###
C>
C> [1] Y Zhao, DG Truhlar,
C> "A new local density functional for main-group thermochemistry,
C> transition metal bonding, thermochemical kinetics, and noncovalent
C> interactions",
C> J. Chem. Phys. <b>125</b>, 194101 (2006), DOI:
C> <a href="https://doi.org/10.1063/1.2370993">
C> 10.1063/1.2370993</a>.
C>
C> [2] Y Zhao, DG Truhlar,
C> "Density functional for spectroscopy: No long-range self-interaction
C> error, good performance for Rydberg and charge-transfer states,
C> and better performance on average than B3LYP for ground states",
C> J. Phys. Chem. A <b>110</b>, 13126-13130 (2006), DOI:
C> <a href="https://doi.org/10.1021/jp066479k">
C> 10.1021/jp066479k</a>.
C>

      Subroutine nwpwxc_x_m06_2x(param, tol_rho, ipol, nq, wght, 
     &                         rho, rgamma, tau, func, Amat, Cmat, Mmat)
c   
c$Id$
c
      implicit none
c
#include "nwpwxc_param.fh"
c
      double precision param(*) !< [Input] Parameters of the functional
                                !< - param( 1): \f$ a_0 \f$
                                !< - param( 2): \f$ a_1 \f$
                                !< - param( 3): \f$ a_2 \f$
                                !< - param( 4): \f$ a_3 \f$
                                !< - param( 5): \f$ a_4 \f$
                                !< - param( 6): \f$ a_5 \f$
                                !< - param( 7): \f$ a_6 \f$
                                !< - param( 8): \f$ a_7 \f$
                                !< - param( 9): \f$ a_8 \f$
                                !< - param(10): \f$ a_9 \f$
                                !< - param(11): \f$ a_10 \f$
                                !< - param(12): \f$ a_11 \f$
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
      double precision at, at10, at11, at0, C1, C2, fL, fNL
      double precision rrho, rho43, rho13, rhoo, rho53
      double precision Gamma2, Gamma
      double precision TauUEG, Tsig, Wsig, W1, W2, W3, W4, W5, W6
      double precision W7, W8, W9, W10, W11, Fsig
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

c     functional derivatives above FFFFFFFFFFFF


cedo       parameter( pi = 3.1415926535897932384626433832795d0 )
         
      
      parameter (F1o3=1.d0/3.d0, F1o4=1.d0/4.d0, F2o3=2.d0/3.d0, 
     &     F3o2=3.d0/2.d0,F13=1.d0/3.d0)
      parameter (F4o3=4.d0/3.d0, F4o9=4.d0/9.d0, F3o5=3.d0/5.d0)
      parameter (F83=8.d0/3.0d0, F23=2.0d0/3.d0, F53=5.d0/3.d0)
      parameter (One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0, 
     &     Five=5.0d0,Six=6.0d0, Seven=7.0d0,
     &     Eight=8.0d0, Nine=9.0d0,F10=10.d0, F11=11.d0)
      pi=acos(-1d0)      

      at0  = param(1)
      at1  = param(2)
      at2  = param(3)
      at3  = param(4)
      at4  = param(5)
      at5  = param(6)
      at6  = param(7)
      at7  = param(8)
      at8  = param(9)
      at9  = param(10)
      at10 = param(11)
      at11 = param(12)
c     if (ijzy.eq.1) then
c       at0=    3.987756D-01
c       at1=    2.548219D-01
c       at2=    3.923994D-01
c       at3=    -2.103655D+00
c       at4=    -6.302147D+00
c       at5=    1.097615D+01
c       at6=    3.097273D+01
c       at7=    -2.318489D+01
c       at8=    -5.673480D+01
c       at9=    2.160364D+01
c       at10=   3.421814D+01
c       at11=   -9.049762D+00
c      elseif (ijzy.eq.2) then
c     Parameters for M06-HF
c       at0=    1.179732D-01
c       at1=    -1.066708D+00
c       at2=    -1.462405D-01
c       at3=    7.481848D+00
c       at4=    3.776679D+00
c       at5=    -4.436118D+01
c       at6=    -1.830962D+01
c       at7=    1.003903D+02
c       at8=    3.864360D+01
c       at9=    -9.806018D+01
c       at10=   -2.557716D+01
c       at11=   3.590404D+01
c      elseif (ijzy.eq.3) then
c     Parameters for M06
c       at0=    5.877943D-01
c       at1=    -1.371776D-01
c       at2=    2.682367D-01
c       at3=    -2.515898D+00
c       at4=    -2.978892D+00
c       at5=    8.710679D+00
c       at6=    1.688195D+01
c       at7=    -4.489724D+00
c       at8=    -3.299983D+01
c       at9=    -1.449050D+01
c       at10=   2.043747D+01
c       at11=   1.256504D+01
c      elseif (ijzy.eq.4) then
c     Parameters for M06-2X
c       at0=    4.600000D-01
c       at1=    -2.206052D-01
c       at2=    -9.431788D-02
c       at3=    2.164494D+00
c       at4=    -2.556466D+00
c       at5=    -1.422133D+01
c       at6=    1.555044D+01
c       at7=    3.598078D+01
c       at8=    -2.722754D+01
c       at9=    -3.924093D+01
c       at10=   1.522808D+01
c       at11=   1.522227D+01
c     endif

      at=1.0d0

      C1 = 0.2195149727645171d0
      C2 = C1/0.804d0 
cedo      DTol=1.0D-8
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
         if(taun.lt.dtol) goto 10
            tauu=tauN 
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
            Fsig =at*(at0 + at1*W1+ at2*W2 + at3*W3
     &          + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &          + at8*W8 + at9*W9 + at10*W10 + at11*W11)

c           Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = rgamma(n,G_TT)
            Gamma = dsqrt(Gamma2)
         if(gamma.lt.dtol) goto 10
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*Fsig*qwght(n)
            func(n)=func(n)+rho43*Ax*(fL+fNL*E)*Fsig
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
            Amat(n,D1_RA)  = Amat(n,D1_RA) + dGGAdR*Fsig
     &                     + (fL+fNL*E)*Ax*rho43*dFdR
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA)
     &                     + Two*dGGAdG*Ax*rho43*Fsig 
            Mmat(n,D1_TA)  = Mmat(n,D1_TA)
     &                     + 0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau

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
             Fsig =at*(at0 + at1*W1+ at2*W2 + at3*W3
     &           + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &           + at8*W8 + at9*W9 + at10*W10 + at11*W11)


c           Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = rgamma(n,G_AA)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)

            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(Fsig)*0.5d0*qwght(n)
            func(n)=func(n)+rho43*Ax*(fL+fNL*E)*(Fsig)*.5d0
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

            Amat(n,D1_RA)  = Amat(n,D1_RA) + (dGGAdR*(Fsig)
     &                     + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA)
     &                     + dGGAdG*Ax*rho43*(Fsig)*0.5d0
            Mmat(n,D1_TA)  = Mmat(n,D1_TA)
     &                     + 0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau

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
             Fsig =at*(at0+at1*W1+ at2*W2 + at3*W3
     &           + at4*W4 + at5*W5 + at6*W6 + at7*W7
     &           + at8*W8 + at9*W9 + at10*W10 + at11*W11)


c           Gamma2 = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            Gamma2 = rgamma(n,G_BB)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
c           Ex = Ex + rho43*Ax*(fL+fNL*E)*(Fsig)*0.5d0*qwght(n)
            func(n)= func(n)+rho43*Ax*(fL+fNL*E)*(Fsig)*.5d0
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

            Amat(n,D1_RB)  = Amat(n,D1_RB) + (dGGAdR*(Fsig)
     &                     + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB)
     &                     + dGGAdG*Ax*rho43*(Fsig)*0.5d0
            Mmat(n,D1_TB)  = Mmat(n,D1_TB)
     &                     + 0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau
c
20      continue
      endif
c
      return
      end




      Subroutine nwpwxc_x_m06_d2()
      call errquit(' xm06: d2 not coded ',0,0)
      return
      end
C> @}
