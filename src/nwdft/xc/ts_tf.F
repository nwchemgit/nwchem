cc AJL/Begin/FDE
      Subroutine ts_tf(tol_rho, fac, lfac, nlfac, rho, Amat, nq, 
     &                    ipol, Ts, qwght, ldew, func)
c
      Implicit none
#include "errquit.fh"
c
#include "stdio.fh"
c
      integer nq, ipol
      double precision fac, Ts
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]
      double precision two_rho32, two_rho33
c
c     Charge Density
c
      double precision rho(nq,(ipol*(ipol+1))/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Partial First Derivatives of the Kinetic Energy Functional
c
      double precision Amat(nq,ipol)
c
c     Compute the partial derivatives of the Kinetic functional of Dirac.
c
      double precision P1, P2, P3, P4, tol_rho
c
c Daniel (5-23-12): P1 and P2 are for restricted calculations, P3 
c and P4 are for unrestricted calculations.
c
cc AJL: I made an error in these definitions. Updating on 23/09/2015
cc     P1 =            (1/2) *(3*PI^2)**(1/3)
cc     P2 =           (3/10) *(3*PI^2)**(1/3)
cc     P3 = 2**(2/3)*  (1/2) *(3*PI^2)**(1/3)
cc     P4 = 2**(2/3)* (3/10) *(3*PI^2)**(1/3) 
cc
cc Calculated 14/10/2014 with Wolfram Alpha
c      Parameter (P1 = 0.15468338631400680D+01)
c      Parameter (P2 = 0.92810031788404078D+00)
c      Parameter (P3 = 0.24554457015685778D+01)
c      Parameter (P4 = 0.14732674209411467D+01)

c New values
c     P1 =            (1/2) *(3*PI^2)**(2/3)
c     P2 =           (3/10) *(3*PI^2)**(2/3)
c     P3 = 2**(2/3)*  (1/2) *(3*PI^2)**(2/3)
c     P4 = 2**(2/3)* (3/10) *(3*PI^2)**(2/3) 
c
c Calculated 23/09/2015 with Wolfram Alpha
      Parameter (P1 = 0.47853900003136530D+01)
      Parameter (P2 = 0.28712340001881918D+01)
      Parameter (P3 = 0.75963331205759952D+01)
      Parameter (P4 = 0.45577998723455971D+01)
      double precision rho13, rhom23, two_third, term
      double precision rho32, rho33
      Parameter (two_third = 2.d0/3.d0)
      integer n
c
      term = 0.d0
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
c Tested against values from Jones and Gunnarsson, 
c doi: 10.1103/RevModPhys.61.689

c Reported Ts for He: 2.5605 (No basis details, HF orbitals)
c Our calculated value: 2.5600 (6-311+G, uncontracted, xc = HFexch)
c
c Reported Ts for Ne: 117.78
c Our calculated value: 117.7226 (6-311+G, uncontracted, and xc = HFexch) 
c
c For spin-polarised implementation we got identical results. All good.
c
         do n = 1, nq
c            write(luout,*) rho(n,1)
            if (rho(n,1).gt.tol_rho)then
               rho13=rho(n,1)**two_third
               Amat(n,1) = Amat(n,1) + rho13
               term = term + rho(n,1)*rho13*qwght(n)
               if(ldew) func(n) = func(n) + rho(n,1)*rho13*fac*P2
            endif
         enddo
         call dscal(nq, P1*fac, Amat(1,1), 1)
c 
         Ts = Ts + P2*term*fac
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
c The only useful initial test for this I could find was for a hydrogen atom,
c as listed in Parr & Yang's book "Density Functional Theory of Atoms and
c Molecules". On page 175 they show the analytical values for Ts for a
c hydrogen atom as being 0.2891 in spin-paired (with 0.5 e occupancy in
c each channel) and 0.4590 for the spin-unpaired (1 e in alpha channel).
c Due to the implementation of spin-paired calculations in NWChem I can't
c test the former, but for the latter the below implementation results in
c a value of 0.4540 (xc = LDA). Of course, our result is dependent on this
c functional but the similarity of the results leaves me in good
c confidence that this has been correctly implemented now.
c
c Additional note for spin-unpaired (23/09/2015):
c I have found I can also reference against the work of Iyengar,
c Ernzerhof, Maximoff and Scuseria, DOI: 10.1103/PhysRevA.63.052508 
c Here they list the energies for a whole range of gas phase atoms, with
c KS orbitals and densities obtained with B88-PW91. For H, they achieved
c Ts of 0.456, compared to 0.4631 (6-311+G uncontracted, same basis as 
c them) for our implementations and so I am confident this is correct.
c
c For an additional test, I did C(triplet) with the same basis and XC.
c In the published case, the value of Ts is 33.983, whereas for our
c implementation we got 34.0596. Again, I can only believe the
c difference is in the numerical implementation that is used in
c solving the KS equations and the accuracy of the settings for each
c type of calculation.
c
         do n = 1,nq
c            write(luout,*) rho(n,1), rho(n,2), rho(n,3)
            if (rho(n,1).gt.tol_rho)then

               rho32 = 0.d0
               rho33 = 0.d0

               if (rho(n,2).gt.tol_rho) rho32 = rho(n,2)**two_third
               if (rho(n,3).gt.tol_rho) rho33 = rho(n,3)**two_third
c Debug: Testing the functionality
c               rho32 = rho(n,2)**two_third
c               rho33 = rho(n,3)**two_third

               Amat(n,1) = Amat(n,1) + P3*rho32*fac
               Amat(n,2) = Amat(n,2) + P3*rho33*fac

               if (ldew)func(n) = func(n) + ((rho(n,2)*rho32)
     &              + (rho(n,3)*rho33))*P4*fac

               term = term+((rho(n,2)*rho32)
     &              + (rho(n,3)*rho33))*qwght(n)

c      write(LuOut,*)'from ts_tf ; rho32, rho33, qwght(n): ',
c     &                           rho32, rho33, qwght(n)
            endif
         enddo
c 
         Ts = Ts + P4*term*fac
      endif
c      write(LuOut,*)'from ts_tf ; p4, term, fac, Ex: ', 
c     &                           P4, term, fac, Ts 
      return
      end
c
cc AJL/End
