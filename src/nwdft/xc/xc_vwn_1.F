      Subroutine xc_vwn_1(tol_rho,fac,lfac,nlfac,rho,Amat,nq,ipol,
     &           Ec,qwght,ldew,func)
c
      implicit none
      integer nq,ipol
      double precision fac,Ec
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]

      integer initial

      Save initial,q,nu0,nu1,nu2,p0,p1

      double precision zero,one,two,four,
     * onethird,fourthirds,onesixth,pi,
     * a1,a2,a3,b1,b2,b3,c1,c2,c3,
     * d1,d2,d3,tol_rho
      Parameter (zero=0.D0,one=1.D0,two=2.D0,four=4.D0)
      Parameter (onethird=1.D0/3.D0,fourthirds=4.D0/3.D0)
      Parameter (onesixth=1.0d0/6.D0)
      Parameter (pi=3.1415926535898D0)
      Parameter (a1=-1.0d0/(3.D0*pi*pi),a2=0.0621814D0,a3=0.5D0*a2)
      Parameter (b1=1.13107D0,b2=3.72744D0,b3=7.06042D0)
      Parameter (c1=13.0045D0,c2=12.9352D0,c3=18.0578D0)
      Parameter (d1=-0.00475840D0,d2=-0.10498D0,d3=-0.32500D0)

c--> Charge Density & Its Cubed Root

      double precision rho(ipol*(ipol+1)/2,nq)

c--> Quadrature Weights

      double precision qwght(nq)

c--> Partial Derivatives of the Correlation Energy Functional

      double precision Amat(ipol,nq)

c--> Functional Parameters (local)

      double precision q(3),nu0(3),nu1(3),nu2(3)
      double precision dalpha,rho3,p0,p1,p2,p3,
     * xs,vc2n,vc3n,decp,alpha,t1,t2,t3,chi_inv,
     * zeta,zeta3,zeta4,ecp,ecf,decf,dfzeta,fzeta,
     * decdx,decdzeta

c--> auxiliary variables

      double precision aux1

      integer i,n

      Data initial/1/
c
c Correlation energy density and potential as determined via the local 
c density functional of Vosko, Wilk and Nusair (Can. J. of Phys. 58, 
c 1200 (1980)), functional I, the Monte Carlo parameterization.
c
c The interpolation formula for paramagnetic portion of the correlation 
c energy has the following form:
c
c             ecp = a2*( t1 + ((2*b2)/q(2))*t3 - 
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
c         d(ecp)/dx = 2*a2*( 1/xs - (xs/chi(xs))*( 1 + b2/(xs - xso) ) )
c
c N.B.:  Converting from Rydberg to Hartree units eliminates a factor of 
c        2 from this formula.
c
c The ferromagnetic portion and the spin stiffness are similarly interpolated.
c
c***************************************************************************

c Define miscellaneous parameters.

       p0 = (1.0d0/(fourthirds*pi))**onethird
       p1 = 0.5D0/(two**onethird - one)
       p2 = 2.D0/(3.D0*(two**onethird - one))
       p3 = 2.25D0*(two**onethird - one)
      if( initial.eq.1 )then

       initial = 0
       
c Define constants used to interpolate the spin stiffness.

       q(1) = sqrt(four*c1 - b1*b1)
       nu0(1) = a1
       nu1(1) = -a1*b1*d1/(d1*d1 + b1*d1 + c1)
       nu2(1) = a1*(two*b1/q(1))*((c1-d1*d1)/(d1*d1 + b1*d1 + c1))

c Define constants used to interpolate the paramagnetic portion.

       q(2) = sqrt(four*c2 - b2*b2)
       nu0(2) = a2
       nu1(2) = -a2*b2*d2/(d2*d2 + b2*d2 + c2)
       nu2(2) = a2*(two*b2/q(2))*((c2-d2*d2)/(d2*d2 + b2*d2 + c2))

c Define constants used to interpolate the ferromagnetic portion.

       q(3) = sqrt(four*c3 - b3*b3)
       nu0(3) = a3
       nu1(3) = -a3*b3*d3/(d3*d3 + b3*d3 + c3)
       nu2(3) = a3*(two*b3/q(3))*((c3-d3*d3)/(d3*d3 + b3*d3 + c3))

c Convert interpolation formulas from Rydberg to Hartree units.

c
c      added factor FAC
c
       do 10 i = 1,3
        nu0(i) = 0.5D0*nu0(i)
        nu1(i) = 0.5D0*nu1(i)
        nu2(i) = 0.5D0*nu2(i)
   10  continue
c       
      end if
c
      if( ipol.eq.1 )then
c
c                      ======> SPIN-RESTRICTED <======
c
       do 100 n = 1,nq
c
c        Compute the interpolation variable.
c
         if(abs(rho(1,n)).lt.tol_rho) goto 100
        rho3=rho(1,n)**onethird
        xs = sqrt(p0/rho3)
c
c Compute the paramagnetic portion of the correlation energy density.
c
        chi_inv = 1.0d0/(xs*xs + b2*xs + c2)
c
        t1 = log( xs*xs*chi_inv )
        t2 = log( (xs - d2)*(xs - d2)*chi_inv )
        t3 = atan( q(2)/(two*xs + b2) )

        Vc2n = nu0(2)*t1 + nu1(2)*t2 + nu2(2)*t3
        Ec   = Ec + Vc2n*rho(1,n)*qwght(n)*fac
        if(ldew) func(n) = func(n) + Vc2n*rho(1,n)*fac
c
c Compute the paramagnetic portion of the correlation potential.
c (Note, converting from Rydberg to Hartree units eliminates a factor of 
c  2 from this formula.)
c
        decp = a2*( 1.0d0/xs - xs*chi_inv*( one + b2/(xs - d2) ) )
c
        Amat(1,n) = Amat(1,n) + (Vc2n - onesixth*xs*decp)*fac
c
  100  continue
c
      else
c
c                      ======> SPIN-UNRESTRICTED <======
c
       do 200 n = 1,nq
         if(abs(rho(1,n)).lt.tol_rho) goto 200
c
c Compute the interpolation variable.
c
        rho3=rho(1,n)**onethird
        xs = sqrt(p0/rho3)
c
c Compute the spin stiffness parameter and its derivative.
c
        chi_inv = 1.0d0/(xs*xs + b1*xs + c1)
c
        t1 = log( xs*xs*chi_inv )
        t2 = log( (xs - d1)*(xs - d1)*chi_inv )
        t3 = atan( q(1)/(two*xs + b1) )
c
        alpha = nu0(1)*t1 + nu1(1)*t2 + nu2(1)*t3
c
        dalpha = a1*( 1.0d0/xs - xs*chi_inv*( one + b1/(xs - d1) ) )
c
c Compute the paramagnetic portion of the correlation energy density and
c its derivative.
c
        chi_inv = 1.0d0/(xs*xs + b2*xs + c2)
c
        t1 = log( xs*xs*chi_inv )
        t2 = log( (xs - d2)*(xs - d2)*chi_inv )
        t3 = atan( q(2)/(two*xs + b2) )
c
        ecp = nu0(2)*t1 + nu1(2)*t2 + nu2(2)*t3
c
        decp = a2*( 1.0d0/xs - xs*chi_inv*( one + b2/(xs - d2) ) )
c
c Compute the ferromagnetic portion of the correlation energy density and
c its derivative.
c
        chi_inv = 1.0d0/(xs*xs + b3*xs + c3)
c
        t1 = log( xs*xs*chi_inv )
        t2 = log( (xs - d3)*(xs - d3)*chi_inv )
        t3 = atan( q(3)/(two*xs + b3) )
c
        ecf = nu0(3)*t1 + nu1(3)*t2 + nu2(3)*t3
c
        decf = a3*( 1.0d0/xs - xs*chi_inv*( one + b3/(xs - d3) ) )
c
c Compute the polarization function and its derivative.
c
        zeta = (rho(2,n) - rho(3,n))/rho(1,n)
c
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
c
c Compute the total correlation energy density.
c
        aux1=ecf-ecp
        Vc3n = ecp + aux1*fzeta
c
        Ec   = Ec + Vc3n*qwght(n)*rho(1,n)*fac
        if(ldew) func(n) = func(n) + Vc3n*rho(1,n)*fac
c
c Compute partial derivatives of Vc3n (which is ec)
c
        decdx = decp + (decf - decp)*fzeta
        decdzeta = aux1*dfzeta
c
c Compute the total correlation potential.
c
        aux1 = Vc3n - onesixth*xs*decdx
        Amat(1,n) = Amat(1,n) + (aux1 + 
     &              decdzeta*(one-zeta))*fac
        Amat(2,n) = Amat(2,n) + (aux1 -
     &              decdzeta*(one+zeta))*fac
c
c Compute common terms.
c
c        eta  = p3*alpha
c        deta = p3*dalpha
c
c        tau  =  ecf -  ecp -  eta
c        dtau = decf - decp - deta
c
c Compute the total correlation energy density.
c
c        Vc3n = ecp + fzeta*(eta + tau*zeta4)
c
c        Ec   = Ec + Vc3n*qwght(n)*rho(1,n)*fac
c        if(ldew) func(n) = func(n) + Vc3n*rho(1,n)*fac
c
c Compute the total correlation potential.
c
c        f1 = Vc3n - onesixth*xs*( decp + fzeta*(deta + dtau*zeta4) )
c
c        f2 = dfzeta*(eta + tau*zeta4) + fzeta*four*tau*zeta3
c
c        Amat(1,n) = Amat(1,n) + (f1 + f2*(one - zeta))*fac
c        Amat(2,n) = Amat(2,n) + (f1 - f2*(one + zeta))*fac
c
  200  continue
c
      end if
c
      return
      end
