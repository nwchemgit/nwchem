      Subroutine xc_vwn_4(tol_rho,fac,lfac,nlfac,rho,Amat,nq,ipol,
     &           Ec,qwght,ldew,func)

      implicit none
      integer nq,ipol
      double precision fac,Ec
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]

      integer initial

      Save initial,qc,nu0c,nu1c,nu2c,qr,nu0r,nu1r,nu2r,p0,p1,p2,p3

      double precision zero,one,two,four,
     * onethird,fourthirds,onesixth,half,pi,
     * a1c,a2c,a3c,b1c,b2c,b3c,c1c,c2c,c3c,
     * d1c,d2c,d3c,
     * a1r,a2r,a3r,b1r,b2r,b3r,c1r,c2r,c3r,
     * d1r,d2r,d3r,tol_rho
      Parameter (zero=0.D0,one=1.D0,two=2.D0,four=4.D0,half=0.5D0)
      Parameter (onethird=1.D0/3.D0,fourthirds=4.D0/3.D0)
      Parameter (onesixth=1.D0/6.D0)
      Parameter (pi=3.1415926535898D0)
      Parameter (a1c=-1.D0/(3.D0*pi*pi),a2c=0.0621814D0,a3c=0.5D0*a2c)
      Parameter (b1c=1.13107D0,b2c=3.72744D0,b3c=7.06042D0)
      Parameter (c1c=13.0045D0,c2c=12.9352D0,c3c=18.0578D0)
      Parameter (d1c=-0.00475840D0,d2c=-0.10498D0,d3c=-0.32500D0)

      Parameter (a1r=-1.D0/(3.D0*pi*pi),a2r=0.0621814D0,a3r=0.5D0*a2r)
      Parameter (b1r=1.06835D0,b2r=13.0720D0,b3r=20.1231D0)
      Parameter (c1r=11.4813D0,c2r=42.7198D0,c3r=101.578D0)
      Parameter (d1r=-0.228344D0,d2r=-0.409286D0,d3r=-0.743294D0)

c--> Charge Density & Its Cubed Root

      double precision rho(ipol*(ipol+1)/2,nq)

c--> Quadrature Weights

      double precision qwght(nq)

c--> Partial Derivatives of the Correlation Energy Functional

      double precision Amat(ipol,nq)

c--> Functional Parameters (local)

      double precision qc(3),nu0c(3),nu1c(3),nu2c(3)
      double precision qr(3),nu0r(3),nu1r(3),nu2r(3)
      double precision rho3,p0,p1,p2,p3,
     * xs,vc2n,vc3n,depc,alphar,alphac,t1,t2,t3,chi_inv,
     * zeta,zeta3,zeta4,epc,efc,defc,dfzeta,
     * eta,deta,tau,dtau,f1,f2,fzeta,
     * dalphar,derpa,dderpadx,dderpadz,decdx,decdzeta

c--> auxiliary variables

      double precision aux1

      integer i,n

      Data initial/1/
c
c Correlation energy density and potential as determined via the local 
c density functional of Vosko, Wilk and Nusair (Can. J. of Phys. 58, 
c 1200 (1980)), functional IV.
c
c The interpolation formula for paramagnetic portion of the correlation 
c energy has the following form:
c
c             ep = a2*( t1 + ((2*b2)/q(2))*t3 - 
c    &                ( (b2*d2)/(d2*d2+b2*d2+c2) )*
c    &               ( t2 + 2*( (b2 + 2*d2)/q(2) )*t3 ) )
c
c             t1,t2,t3,q are defined below
c             a2,b2,c2,d2 are the fitting parameters specified for
c                         paramagnetic, ferromagnetic and spin stiffness
c
c The ferromagnetic portion and the spin stiffness are similarly interpolated.
c
c The derivative of the interpolation formula for paramagnetic portion of the 
c correlation energy has the following form:
c
c         d(ep)/dx = 2*a2*( 1/xs - (xs/chi(xs))*( 1 + b2/(xs - xso) ) )
c
c N.B.:  Converting from Rydberg to Hartree units eliminates a factor of 
c        2 from this formula.
c
c The ferromagnetic portion and the spin stiffness are similarly interpolated.
c
c***************************************************************************

c Define miscellaneous parameters.

      if( initial.eq.1 )then
       p0 = (1.D0/(fourthirds*pi))**onethird
       p1 = 0.5D0/(two**onethird - one)
       p2 = 2.D0/(3.D0*(two**onethird - one))
       p3 = 2.25D0*(two**onethird - one)
       
c Define constants used to interpolate the paramagnetic portion 
c (Monte Carlo parameters). Convert interpolation formulas from Rydberg 
c to Hartree units.

       qc(2) = sqrt(four*c2c - b2c*b2c)
       nu0c(2) = half*a2c
       nu1c(2) = -half*a2c*b2c*d2c/(d2c*d2c + b2c*d2c + c2c)
       nu2c(2) = half*a2c*(two*b2c/qc(2))*((c2c-d2c*d2c)/(d2c*d2c + 
     &           b2c*d2c + c2c))

       if(ipol.ne.1) then 
c Define constants used to interpolate the spin stiffness (RPA 
c parameters). Convert interpolation formulas from Rydberg to Hartree 
c units.

        qr(1) = sqrt(four*c1r - b1r*b1r)
        nu0r(1) = half*a1r
        nu1r(1) = -half*a1r*b1r*d1r/(d1r*d1r + b1r*d1r + c1r)
        nu2r(1) = half*a1r*(two*b1r/qr(1))*((c1r-d1r*d1r)/(d1r*d1r + 
     &            b1r*d1r + c1r))

c Define constants used to interpolate the ferromagnetic portion 
c (Monte Carlo parameters). Convert interpolation formulas from Rydberg 
c to Hartree units.

        qc(3) = sqrt(four*c3c - b3c*b3c)
        nu0c(3) = half*a3c
        nu1c(3) = -half*a3c*b3c*d3c/(d3c*d3c + b3c*d3c + c3c)
        nu2c(3) = half*a3c*(two*b3c/qc(3))*((c3c-d3c*d3c)/(d3c*d3c + 
     &            b3c*d3c + c3c))

       endif

       initial = 0
      end if

      if( ipol.eq.1 )then

c                      ======> SPIN-RESTRICTED <======

       do 100 n = 1,nq

c        Compute the interpolation variable.
        if(abs(rho(1,n)).lt.tol_rho) goto 100
        rho3=rho(1,n)**onethird
        xs = sqrt(p0/rho3)

c Compute the paramagnetic portion of the correlation energy density 
c (Monte Carlo parameters)

        chi_inv = 1.D0/(xs*xs + b2c*xs + c2c)

        t1 = log( xs*xs*chi_inv )
        t2 = log( (xs - d2c)*(xs - d2c)*chi_inv )
        t3 = atan( qc(2)/(two*xs + b2c) )

        Vc2n = nu0c(2)*t1 + nu1c(2)*t2 + nu2c(2)*t3
        Ec   = Ec + Vc2n*rho(1,n)*qwght(n)*fac
        if(ldew) func(n) = func(n) + Vc2n*rho(1,n)*fac
c Compute the paramagnetic portion of the correlation potential.
c (Note, converting from Rydberg to Hartree units eliminates a factor of 
c  2 from this formula.)

        depc = a2c*( 1.D0/xs - xs*chi_inv*( one + b2c/(xs - d2c) ) )

        Amat(1,n) = Amat(1,n) + (Vc2n - onesixth*xs*depc)*fac

  100  continue

      else

c                      ======> SPIN-UNRESTRICTED <======

       do 200 n = 1,nq
         if(abs(rho(1,n)).lt.tol_rho) goto 200

c Compute the interpolation variable.

        rho3=rho(1,n)**onethird
        xs = sqrt(p0/rho3)

c Compute the spin stiffness parameter and its derivative (RPA 
c parameterization)

        chi_inv = 1.D0/(xs*xs + b1r*xs + c1r)

        t1 = log( xs*xs*chi_inv )
        t2 = log( (xs - d1r)*(xs - d1r)*chi_inv )
        t3 = atan( qr(1)/(two*xs + b1r) )

        alphar = nu0r(1)*t1 + nu1r(1)*t2 + nu2r(1)*t3

        dalphar = a1r*( 1.D0/xs - xs*chi_inv*( one + b1r/(xs - d1r) ) )

c Compute the paramagnetic portion of the correlation energy density and
c its derivative (Monte Carlo parameterization)

        chi_inv = 1.D0/(xs*xs + b2c*xs + c2c)

        t1 = log( xs*xs*chi_inv )
        t2 = log( (xs - d2c)*(xs - d2c)*chi_inv )
        t3 = atan( qc(2)/(two*xs + b2c) )

        epc = nu0c(2)*t1 + nu1c(2)*t2 + nu2c(2)*t3

        depc = a2c*( 1.D0/xs - xs*chi_inv*( one + b2c/(xs - d2c) ) )

c Compute the ferromagnetic portion of the correlation energy density 
c and its derivative (Monte Carlo parameterization)

        chi_inv = 1.D0/(xs*xs + b3c*xs + c3c)

        t1 = log( xs*xs*chi_inv )
        t2 = log( (xs - d3c)*(xs - d3c)*chi_inv )
        t3 = atan( qc(3)/(two*xs + b3c) )

        efc = nu0c(3)*t1 + nu1c(3)*t2 + nu2c(3)*t3

        defc = a3c*( 1.D0/xs - xs*chi_inv*( one + b3c/(xs - d3c) ) )

c Compute the polarization function and its derivative.

        zeta = (rho(2,n) - rho(3,n))/rho(1,n)

        zeta3 = zeta*zeta*zeta
        zeta4 = zeta3*zeta
        if(zeta.gt.(one-1.d-10)) then
          fzeta = (two**fourthirds-two)*p1
          dfzeta = (two**onethird)*p2
        elseif(zeta.lt.(-one+1.d-10)) then
          fzeta = (two**fourthirds-two)*p1
          dfzeta = -(two**onethird)*p2
        else
          fzeta  = ( (one + zeta)**fourthirds + 
     &             (one - zeta)**fourthirds - two )*p1
          dfzeta = ((one+zeta)**onethird-(one-zeta)**onethird)*p2
        endif

c Compute Delta E_RPA(xz) and its partial derivatives

        eta  = p3*alphar
        deta = p3*dalphar

        tau  =  efc -  epc -  eta
        dtau = defc - depc - deta

        aux1 = eta + tau*zeta4
        derpa = fzeta*aux1
        dderpadx = fzeta*(deta + dtau*zeta4)
        dderpadz = dfzeta*aux1 + four*fzeta*tau*zeta3

c Compute the total correlation energy density.

        Vc3n = epc + derpa

        Ec   = Ec + Vc3n*qwght(n)*rho(1,n)*fac
        if(ldew) func(n) = func(n) + Vc3n*rho(1,n)*fac

c Compute partial derivatives of Vc3n (which is ec)
        
        decdx = depc + dderpadx 
        decdzeta = dderpadz

c Compute the total correlation potential.

        aux1 = Vc3n - onesixth*xs*decdx 
        Amat(1,n) = Amat(1,n) + (aux1 + 
     &              decdzeta*(one-zeta))*fac
        Amat(2,n) = Amat(2,n) + (aux1 -
     &              decdzeta*(one+zeta))*fac

c Compute common terms.
c
c        eta  = p3*alpha
c        deta = p3*dalpha
c
c        tau  =  ef -  ep -  eta
c        dtau = def - dep - deta
c
c Compute the total correlation energy density.
c
c        Vc3n = ep + fzeta*(eta + tau*zeta4)
c
c        Ec   = Ec + Vc3n*qwght(n)*rho(1,n)*fac
c        if(ldew) func(n) = func(n) + Vc3n*rho(1,n)*fac
c
c Compute the total correlation potential.
c
c        f1 = Vc3n - onesixth*xs*( dep + fzeta*(deta + dtau*zeta4) )
c
c        f2 = dfzeta*(eta + tau*zeta4) + fzeta*four*tau*zeta3
c
c        Amat(1,n) = Amat(1,n) + (f1 + f2*(one - zeta))*fac
c        Amat(2,n) = Amat(2,n) + (f1 - f2*(one + zeta))*fac

  200  continue

      end if

      return
      end
