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
c                              ijzy - 5 revM06-L
c                              ijzy - 6 revM06
c     References: 
c     [a]	Zhao, Y. and  Truhlar, D. G. J. Chem. Phys. 2006, 125, 194101;
c     [b]       Zhao, Y. and  Truhlar, D. G. J. Phys. Chem. A (2006),110(49),13126-13130.    
c     [c]       Grafenstein, J., Izotov, D. and Cremer, D. J. Chem. Phys.  2007, 127, 214103.


      Subroutine xc_xm06(tol_rho, fac,lfac,nlfac, rho, delrho, 
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
      double precision at, at10, at11, at0, C1, C2, fL, fNL
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

        at0= 0d0
        at1= 0d0
        at2= 0d0
        at3= 0d0
        at4= 0d0
        at5= 0d0
        at6= 0d0
        at7= 0d0
        at8= 0d0
        at9= 0d0
        at10=0d0
        at11=0d0
      if (ijzy.eq.1) then
        at0=    3.987756D-01
        at1=    2.548219D-01
        at2=    3.923994D-01
        at3=    -2.103655D+00
        at4=    -6.302147D+00
        at5=    1.097615D+01
        at6=    3.097273D+01
        at7=    -2.318489D+01
        at8=    -5.673480D+01
        at9=    2.160364D+01
        at10=   3.421814D+01
        at11=   -9.049762D+00
       elseif (ijzy.eq.2) then
C     Parameters for M06-HF
        at0=    1.179732D-01
        at1=    -1.066708D+00
        at2=    -1.462405D-01
        at3=    7.481848D+00
        at4=    3.776679D+00
        at5=    -4.436118D+01
        at6=    -1.830962D+01
        at7=    1.003903D+02
        at8=    3.864360D+01
        at9=    -9.806018D+01
        at10=   -2.557716D+01
        at11=   3.590404D+01
       elseif (ijzy.eq.3) then
C     Parameters for M06
        at0=    5.877943D-01
        at1=    -1.371776D-01
        at2=    2.682367D-01
        at3=    -2.515898D+00
        at4=    -2.978892D+00
        at5=    8.710679D+00
        at6=    1.688195D+01
        at7=    -4.489724D+00
        at8=    -3.299983D+01
        at9=    -1.449050D+01
        at10=   2.043747D+01
        at11=   1.256504D+01
       elseif (ijzy.eq.4) then
C     Parameters for M06-2X
        at0=    4.600000D-01
        at1=    -2.206052D-01
        at2=    -9.431788D-02
        at3=    2.164494D+00
        at4=    -2.556466D+00
        at5=    -1.422133D+01
        at6=    1.555044D+01
        at7=    3.598078D+01
        at8=    -2.722754D+01
        at9=    -3.924093D+01
        at10=   1.522808D+01
        at11=   1.522227D+01
        elseif (ijzy.eq.5) then
C      Parameters for revM06-L
        at0=    1.423227252D+00
        at1=    4.71820438D-01
        at2=    -1.67555701D-01
        at3=    -2.50154262D-01
        at4=    6.2487588D-02
        at5=    7.3350124D-01
        at6=    -2.359736776D+00
        at7=    -1.436594372D+00
        at8=    4.44643793D-01
        at9=    1.529925054D+00
        at10=   2.053941717D+00
        at11=   -3.6536031D-02
        elseif (ijzy.eq.6) then
C      Parameters for revM06
        at0=    6.511394014D-01
        at1=    -1.214497763D-01
        at2=    -1.367041135D-01
        at3=    3.987218551D-01
        at4=    6.056741356D-01
        at5=    -2.379738662D+00
        at6=    -1.492098351D+00
        at7=    3.031473420D+00
        at8=    5.149637108D-01
        at9=    2.633751911D+00
        at10=   9.886749252D-01
        at11=   -4.243714128D+00
      else
        call errquit("xc_xm06: illegal value of ijzy",ijzy,UERR)
      endif

      at=1.0d0
      if (ijzy.ne.4) call xc_xvs98(tol_rho, 1.0d0, lfac, 
     &           nlfac, rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,tau,Mmat,ijzy+1)


      C1 = 0.2195149727645171d0
      C2 = C1/0.804d0 
cedo      DTol=1.0D-8
      DTol=tol_rho
C
C     Scale factors for local and non-local contributions.
C
      fL  =  1.0d0
      fNL =  1.0d0
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

            Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            Gamma = dsqrt(Gamma2)
            if(gamma.lt.dtol) goto 10
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
            Ex = Ex + rho43*Ax*(fL+fNL*E)*Fsig*qwght(n)
            if(ldew) func(n)=func(n)+rho43*Ax*(fL+fNL*E)*Fsig
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
            Amat(n,1) = Amat(n,1) + dGGAdR*Fsig
     &        + (fL+fNL*E)*Ax*rho43*dFdR
            Cmat(n,1)=  Cmat(n,1) + 
     &                    Two*dGGAdG*Ax*rho43*Fsig 
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
             if(tauu.lt.dtol) goto 25
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


            Gamma2 = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if(Gamma.lt.dtol) goto 25

            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
            Ex = Ex + rho43*Ax*(fL+fNL*E)*(Fsig)*0.5d0*qwght(n)
            if(ldew) func(n)=
     =           func(n)+rho43*Ax*(fL+fNL*E)*(Fsig)*.5d0
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

            Amat(n,1) = Amat(n,1) + (dGGAdR*(Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,1)=  Cmat(n,1) + 
     &                      dGGAdG*Ax*rho43*(Fsig)*0.5d0
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
             if(tauu.lt.dtol) goto 20
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


            Gamma2 = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            Gamma2 = Four*Gamma2
            Gamma = dsqrt(Gamma2)
            if(Gamma.lt.dtol) goto 20
            s      = Cs*Gamma/rho43
            s2     = s*s
            En     = C1*s2
            Ed     = One + C2*s2
            E      = En/Ed
            Ex = Ex + rho43*Ax*(fL+fNL*E)*(Fsig)*0.5d0*qwght(n)
            if(ldew) func(n)=
     =           func(n)+rho43*Ax*(fL+fNL*E)*(Fsig)*.5d0
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

            Amat(n,2) = Amat(n,2) + (dGGAdR*(Fsig)
     &        + (fL+fNL*E)*Ax*rho43*dFdR)*0.5d0
            Cmat(n,3)=  Cmat(n,3) + 
     &                   dGGAdG*Ax*rho43*(Fsig)*0.5d0
            Mmat(n,2)=  Mmat(n,2) +
     &                  0.5d0*0.5d0*rho43*Ax*(fL+fNL*E)*dFdTau
     

c
20      continue
      endif
c
      return
      end




      Subroutine xc_xm06_d2()
      call errquit(' xm06: d2 not coded ',0,0)
      return
      end

