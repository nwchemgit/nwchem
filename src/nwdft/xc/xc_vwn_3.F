      Subroutine xc_vwn_3(tol_rho, fac, lfac, nlfac, rho, Amat, nq, 
     &                    ipol, Ec, qwght, ldew, func)
*
* $Id: xc_vwn_3.F,v 1.5 1997-11-01 00:10:04 d3e129 Exp $
*
      implicit none
c
      integer nq, ipol
      double precision fac, Ec
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]
c
      integer initial
      Save initial, qc, nu0c, nu1c, nu2c, qr, nu0r, nu1r, nu2r, 
     &     p0, p1, p2, p3
c
      double precision zero, one, two, four, onethird, fourthirds, 
     &                 onesixth, half, pi, a1c, a2c, a3c, b1c, b2c, 
     &                 b3c, c1c, c2c, c3c, d1c, d2c, d3c, a1r, a2r, 
     &                 a3r, b1r, b2r, b3r, c1r, c2r, c3r, d1r, d2r, 
     &                 d3r, tol_rho
      Parameter (zero = 0.D0, one = 1.D0, two = 2.D0, four = 4.D0, 
     &           half = 0.5D0)
      Parameter (onethird = 1.D0/3.D0, fourthirds = 4.D0/3.D0)
      Parameter (onesixth = 1.D0/6.D0)
      Parameter (pi = 3.1415926535898D0)
      Parameter (a1c = -1.D0/(3.D0*pi*pi), a2c = 0.0621814D0, 
     &           a3c = 0.5D0*a2c)
      Parameter (b1c = 1.13107D0, b2c = 3.72744D0, b3c = 7.06042D0)
      Parameter (c1c = 13.0045D0, c2c = 12.9352D0, c3c = 18.0578D0)
      Parameter (d1c = -0.00475840D0, d2c = -0.10498D0, 
     &           d3c = -0.32500D0)
      Parameter (a1r = -1.D0/(3.D0*pi*pi), a2r = 0.0621814D0, 
     &           a3r = 0.5D0*a2r)
      Parameter (b1r = 1.06835D0, b2r = 13.0720D0, b3r = 20.1231D0)
      Parameter (c1r = 11.4813D0, c2r = 42.7198D0, c3r = 101.578D0)
      Parameter (d1r = -0.228344D0, d2r = -0.409286D0, 
     &           d3r = -0.743294D0)
c
c     Charge Density & Its Cubed Root
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Partial Derivatives of the Correlation Energy Functional
c
      double precision Amat(nq,ipol)
c
c     Functional Parameters (local)
c
      double precision qc(3), nu0c(3), nu1c(3), nu2c(3)
      double precision qr(3), nu0r(3), nu1r(3), nu2r(3)
      double precision rho3, p0, p1, p2, p3, xs, vc2n, vc3n, depr, 
     &                 depc, alphar, alphac, t1, t2, t3, chi_inv, 
     &                 zeta, zeta3, zeta4, epr, epc, efr, efc, defr, 
     &                 defc, dfzeta, eta, deta, tau, dtau, f1, f2, 
     &                 fzeta, dalphar, dalphac, derpa, dderpadx, 
     &                 dderpadz, fmpc, fmpr, dfmpc, dfmpr, coverr, 
     &                 decdx, decdzeta
c
c     auxiliary variables
c
      double precision aux1
      integer i, n
c
      Data initial/1/
c
c Correlation energy density and potential as determined via the local 
c density functional of Vosko, Wilk and Nusair (Can. J. of Phys. 58, 
c 1200 (1980)), functional III.
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
c
c     Define miscellaneous parameters.
c
      if (initial.eq.1)then
         p0 = (1.D0/(fourthirds*pi))**onethird
         p1 = 0.5D0/(two**onethird - one)
         p2 = 2.D0/(3.D0*(two**onethird - one))
         p3 = 2.25D0*(two**onethird - one)
c       
c        Define constants used to interpolate the paramagnetic portion 
c        (Monte Carlo parameters). Convert interpolation formulas from Rydberg 
c        to Hartree units.
c
         qc(2) = sqrt(four*c2c - b2c*b2c)
         nu0c(2) = half*a2c
         nu1c(2) = -half*a2c*b2c*d2c/(d2c*d2c + b2c*d2c + c2c)
         nu2c(2) = half*a2c*(two*b2c/qc(2))*((c2c-d2c*d2c)/(d2c*d2c + 
     &             b2c*d2c + c2c))
c
c        Define constants used to interpolate the paramagnetic portion 
c        (RPA parameters). Convert interpolation formulas from Rydberg 
c        to Hartree units.
c
         qr(2) = sqrt(four*c2r - b2r*b2r)
         nu0r(2) = half*a2r
         nu1r(2) = -half*a2r*b2r*d2r/(d2r*d2r + b2r*d2r + c2r)
         nu2r(2) = half*a2r*(two*b2r/qr(2))*((c2r-d2r*d2r)/(d2r*d2r + 
     &             b2r*d2r + c2r))

         if (ipol.ne.1) then 
c
c           Define constants used to interpolate the spin stiffness (Monte Carlo
c           parameters). Convert interpolation formulas from Rydberg to Hartree 
c           units.
c
            qc(1) = sqrt(four*c1c - b1c*b1c)
            nu0c(1) = half*a1c
            nu1c(1) = -half*a1c*b1c*d1c/(d1c*d1c + b1c*d1c + c1c)
            nu2c(1) = half*a1c*(two*b1c/qc(1))*((c1c-d1c*d1c)/(d1c*d1c + 
     &                b1c*d1c + c1c))
c
c           Define constants used to interpolate the spin stiffness (RPA 
c           parameters). Convert interpolation formulas from Rydberg to Hartree 
c           units.
c
            qr(1) = sqrt(four*c1r - b1r*b1r)
            nu0r(1) = half*a1r
            nu1r(1) = -half*a1r*b1r*d1r/(d1r*d1r + b1r*d1r + c1r)
            nu2r(1) = half*a1r*(two*b1r/qr(1))*((c1r-d1r*d1r)/(d1r*d1r + 
     &                b1r*d1r + c1r))
c
c           Define constants used to interpolate the ferromagnetic portion 
c           (Monte Carlo parameters). Convert interpolation formulas from Rydberg 
c           to Hartree units.
c
            qc(3) = sqrt(four*c3c - b3c*b3c)
            nu0c(3) = half*a3c
            nu1c(3) = -half*a3c*b3c*d3c/(d3c*d3c + b3c*d3c + c3c)
            nu2c(3) = half*a3c*(two*b3c/qc(3))*((c3c-d3c*d3c)/(d3c*d3c + 
     &                b3c*d3c + c3c))
c
c           Define constants used to interpolate the ferromagnetic portion 
c           (RPA parameters). Convert interpolation formulas from Rydberg 
c           to Hartree units.
c
            qr(3) = sqrt(four*c3r - b3r*b3r)
            nu0r(3) = half*a3r
            nu1r(3) = -half*a3r*b3r*d3r/(d3r*d3r + b3r*d3r + c3r)
            nu2r(3) = half*a3r*(two*b3r/qr(3))*((c3r-d3r*d3r)/(d3r*d3r + 
     &                b3r*d3r + c3r))
         endif
         initial = 0
      endif
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do 100 n = 1, nq
c
c           Compute the interpolation variable.
c
            if (abs(rho(n,1)).lt.tol_rho)goto 100
            rho3 = rho(n,1)**onethird
            xs = sqrt(p0/rho3)
c
c           Compute the paramagnetic portion of the correlation energy density 
c           (Monte Carlo parameters)
c
            chi_inv = 1.D0/(xs*xs + b2c*xs + c2c)
            t1 = log(xs*xs*chi_inv)
            t2 = log((xs - d2c)*(xs - d2c)*chi_inv)
            t3 = atan(qc(2)/(two*xs + b2c))
            Vc2n = nu0c(2)*t1 + nu1c(2)*t2 + nu2c(2)*t3
            Ec   = Ec + Vc2n*rho(n,1)*qwght(n)*fac
            if (ldew) func(n) = func(n) + Vc2n*rho(n,1)*fac
c
c           Compute the paramagnetic portion of the correlation potential.
c           (Note, converting from Rydberg to Hartree units eliminates a factor of 
c           2 from this formula.)
c
            depc = a2c*(1.D0/xs - xs*chi_inv*(one + b2c/(xs - d2c)))
            Amat(n,1) = Amat(n,1) + (Vc2n - onesixth*xs*depc)*fac
  100    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 200 n = 1, nq
            if (abs(rho(n,1)).lt.tol_rho) goto 200
c
c           Compute the interpolation variable.
c
            rho3 = rho(n,1)**onethird
            xs = sqrt(p0/rho3)
c
c           Compute the spin stiffness parameter and its derivative (RPA 
c           parameterization)
c
            chi_inv = 1.D0/(xs*xs + b1r*xs + c1r)
            t1 = log(xs*xs*chi_inv)
            t2 = log((xs - d1r)*(xs - d1r)*chi_inv)
            t3 = atan(qr(1)/(two*xs + b1r))
            alphar = nu0r(1)*t1 + nu1r(1)*t2 + nu2r(1)*t3
            dalphar = a1r*(1.D0/xs - xs*chi_inv*(one + b1r/(xs - d1r)))
c
c           Compute the paramagnetic portion of the correlation energy density and
c           its derivative (Monte Carlo parameterization)
c
            chi_inv = 1.D0/(xs*xs + b2c*xs + c2c)
            t1 = log(xs*xs*chi_inv)
            t2 = log((xs - d2c)*(xs - d2c)*chi_inv)
            t3 = atan(qc(2)/(two*xs + b2c))
            epc = nu0c(2)*t1 + nu1c(2)*t2 + nu2c(2)*t3
            depc = a2c*(1.D0/xs - xs*chi_inv*(one + b2c/(xs - d2c)))
c
c           Compute the paramagnetic portion of the correlation energy density and
c           its derivative (RPA parameterization)
c
            chi_inv = 1.D0/(xs*xs + b2r*xs + c2r)
            t1 = log(xs*xs*chi_inv)
            t2 = log((xs - d2r)*(xs - d2r)*chi_inv)
            t3 = atan(qr(2)/(two*xs + b2r))
            epr = nu0r(2)*t1 + nu1r(2)*t2 + nu2r(2)*t3
            depr = a2r*(1.D0/xs - xs*chi_inv*(one + b2r/(xs - d2r)))
c
c           Compute the ferromagnetic portion of the correlation energy density 
c           and its derivative (Monte Carlo parameterization)
c
            chi_inv = 1.D0/(xs*xs + b3c*xs + c3c)
            t1 = log(xs*xs*chi_inv)
            t2 = log((xs - d3c)*(xs - d3c)*chi_inv)
            t3 = atan(qc(3)/(two*xs + b3c))
            efc = nu0c(3)*t1 + nu1c(3)*t2 + nu2c(3)*t3
            defc = a3c*(1.D0/xs - xs*chi_inv*(one + b3c/(xs - d3c)))
c
c           Compute the ferromagnetic portion of the correlation energy density 
c           and its derivative (RPA parameterization)
c
            chi_inv = 1.D0/(xs*xs + b3r*xs + c3r)
            t1 = log(xs*xs*chi_inv)
            t2 = log((xs - d3r)*(xs - d3r)*chi_inv)
            t3 = atan(qr(3)/(two*xs + b3r))
            efr = nu0r(3)*t1 + nu1r(3)*t2 + nu2r(3)*t3
            defr = a3r*(1.D0/xs - xs*chi_inv*(one + b3r/(xs - d3r)))
c
c           Compute the polarization function and its derivative.
c
            zeta = (rho(n,2) - rho(n,3))/rho(n,1)
            zeta3 = zeta*zeta*zeta
            zeta4 = zeta3*zeta
            if (zeta.gt.(one-1.d-10)) then
               fzeta = (two**fourthirds-two)*p1
               dfzeta = (two**onethird)*p2
            elseif (zeta.lt.(-one+1.d-10)) then
               fzeta = (two**fourthirds-two)*p1
               dfzeta = -(two**onethird)*p2
            else
               fzeta  = ((one + zeta)**fourthirds + 
     &                   (one - zeta)**fourthirds - two)*p1
               dfzeta = ((one+zeta)**onethird - 
     &                   (one-zeta)**onethird)*p2
            endif
c
c           Compute Delta E_RPA(xz) and its partial derivatives
c
            eta  = p3*alphar
            deta = p3*dalphar
            tau  =  efr -  epr -  eta
            dtau = defr - depr - deta
            aux1 = eta + tau*zeta4
            derpa = fzeta*aux1
            dderpadx = fzeta*(deta + dtau*zeta4)
            dderpadz = dfzeta*aux1 + four*fzeta*tau*zeta3
c
c           Compute the total correlation energy density.
c
            fmpc = efc - epc
            fmpr = efr - epr
            coverr = fmpc/fmpr
            Vc3n = epc + coverr*derpa 
            Ec   = Ec + Vc3n*qwght(n)*rho(n,1)*fac
            if (ldew)func(n) = func(n) + Vc3n*rho(n,1)*fac
c
c           Compute partial derivatives of Vc3n (which is ec)
c        
            dfmpc = defc - depc
            dfmpr = defr - depr
            decdx = depc + dderpadx*coverr + derpa*(dfmpc/fmpr - 
     &              coverr*dfmpr/fmpr)
            decdzeta = dderpadz*coverr
c
c           Compute the total correlation potential.
c
            aux1 = Vc3n - onesixth*xs*decdx
            Amat(n,1) = Amat(n,1) + (aux1 + 
     &                  decdzeta*(one-zeta))*fac
            Amat(n,2) = Amat(n,2) + (aux1 -
     &                  decdzeta*(one+zeta))*fac
c
  200    continue
      endif
c
      return
      end
