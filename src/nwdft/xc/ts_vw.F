cc AJL/Begin/FDE
c
      Subroutine ts_vw(tol_rho, fac, lfac, nlfac, rho, delrho,
     &                 Amat, Cmat, nq, ipol, Ts, qwght, ldew, func)
c
      Implicit none
#include "errquit.fh"
c
#include "stdio.fh"
c      
#include "dft2drv.fh"
c#include "dft3drv.fh"
c
      integer nq, ipol
      double precision fac, Ts
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,(ipol*(ipol+1))/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Charge Density Laplacian
c
c      double precision del2rho(nq,ipol*(ipol+1)/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Partial First Derivatives of the Kinetic Energy Functional
c
      double precision Amat(nq,ipol)
c      double precision Cmat(nq,3,ipol)
      double precision Cmat(nq,*)
c
c     Compute the partial derivatives of the Kinetic functional of Von Weiszacker
c
      double precision tol_rho
      double precision one_over_rho
c
      double precision delrho_2
      double precision temp_energy, temp_potential, term
      double precision one_eighth, one_quarter
      Parameter (one_eighth = 1.d0/8.d0)
      Parameter (one_quarter = 2.d0*one_eighth)
      integer n
c Spin-polarised variables
      double precision delrho32_2
      double precision delrho33_2
      double precision temp_energy1, temp_energy2
c      double precision temp_potential1, temp_potential2

      double precision one_over_rho32, one_over_rho33
c
      term = 0.d0
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
c Tested against values from Acharya et al. (PNAS, 1980),
c doi: 10.1073/pnas.77.12.6978

c Reported Ts for He: 2.86 (No basis details, HF orbitals)
c Our calculated value: 2.8599 (6-311+G, uncontracted, xc = HFexch)
c
c Reported Ts for Ne: 90.63
c Our calculated value: 90.6048 (6-311+G, uncontracted, and xc = HFexch) 
c
c For spin-polarised implementation we got identical results. All good.
c

         do 10 n = 1, nq
c            write(luout,*) rho(n,1), delrho(n,1,1), del2rho(n,1)
            if (rho(n,1).gt.tol_rho)then
              one_over_rho = 1.d0/rho(n,1)
              delrho_2= (delrho(n,1,1)*delrho(n,1,1))
     &                + (delrho(n,2,1)*delrho(n,2,1)) 
     &                + (delrho(n,3,1)*delrho(n,3,1))
c              delrho_2 = max(delrho_2,tol_rho_2)
c I think this is just in case of numerical instability
              temp_energy = delrho_2*one_over_rho
              term = term + temp_energy*qwght(n)
c
cc Original derivation of potential, using laplacian if density
c               temp_potential = 2*max(del2rho(n,1),tol_rho)
c I think this is just in case of numerical instability also!
c              temp_potential = 2*del2rho(n,1)
c              temp_potential = temp_energy - temp_potential
c \frac{\partial Tw[p]}{\partial \rho}
c = \frac{\partial tw}{\partial \rho} -
c   \triangledown . \frac{\partial tw}{\partial \triangledown \rho}
c              Amat(n,1) =Amat(n,1) + 
c     &               one_eighth*temp_potential*one_over_rho*fac

c
ccc Updated version, using method of Pople and Johnson (1992, CPL)
cc \frac{\partial tw}{\partial \rho} = 
              Amat(n,1) = Amat(n,1) - 
     &               one_eighth*temp_energy*one_over_rho*fac
cc \frac{\partial tw}{\partial |\triangledown \rho|^{2}} = \frac{\partial tw}{\partial\gamma^{\alpha\alpha}}
cc                                                       = \frac{1}{8 \rho} 
              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) +
cc A factor of two is missing when converted from being \partial 
cc \gamma^{\alpha\alpha} to \partial \triangledown \rho.
cc Therefore we account this here in to our numbers, and it seems to
cc produce the correct results.
c     &        one_eighth*one_over_rho*fac
     &        one_quarter*one_over_rho*fac
            endif
   10    continue
c 
         Ts = Ts + one_eighth*term*fac

      else
c
c        ======> SPIN-UNRESTRICTED <======
c
c Tested against values from Iyengar,
c Ernzerhof, Maximoff and Scuseria, DOI: 10.1103/PhysRevA.63.052508,

c Reported Ts for H: 0.500 (6-311+G, uncontracted, xc = B88-PW91)
c Our calculated value: 0.5088 (6-311+G, uncontracted, xc = B88-PW91)
c
c Reported Ts for C: 32.197 (triplet, same setup as above)
c Our calculated value: 34.0596 (triplet, same setup as above). 
c
c As with the TF functional, I put these differences down to numerical
c implementation seeing as the cited reference is ~15 years old. 
c
c For MPI implementation we got identical results. All good.
c
         do 20 n = 1, nq
c            write(luout,*) rho(n,1), delrho(n,1,1)+delrho(n,1,2),
c     &                     del2rho(n,1), del2rho(n,2), del2rho(n,3)
            if (rho(n,1).gt.tol_rho) then

                temp_energy1 = 0.0d0
                temp_energy2 = 0.0d0

                one_over_rho32 = 1.d0/rho(n,2)
                one_over_rho33 = 1.d0/rho(n,3)

                delrho32_2 = (delrho(n,1,1)*delrho(n,1,1))
     &              + (delrho(n,2,1)*delrho(n,2,1))
     &              + (delrho(n,3,1)*delrho(n,3,1))          

                delrho33_2 = (delrho(n,1,2)*delrho(n,1,2))
     &              + (delrho(n,2,2)*delrho(n,2,2))
     &              + (delrho(n,3,2)*delrho(n,3,2))

               if (rho(n,2).gt.tol_rho) 
     &           temp_energy1 = delrho32_2*one_over_rho32
               if (rho(n,3).gt.tol_rho)
     &           temp_energy2 = delrho33_2*one_over_rho33
               temp_energy = temp_energy1 + temp_energy2
               term = term + temp_energy*qwght(n)

c               temp_potential1 = 0.5d0*temp_energy1
c               temp_potential2 = 0.5d0*temp_energy2
c
c Needs correcting to use Cmat               
c               temp_potential1 = temp_potential1-del2rho(n,1)     
c               temp_potential2 = temp_potential2-del2rho(n,2)

               if (rho(n,2).gt.tol_rho) then
                 Amat(n,1) =Amat(n,1) - 
c     &             one_quarter*temp_potential1*one_over_rho32*fac
     &             one_eighth*temp_energy1*one_over_rho32*fac
                 Cmat(n,D1_GAA) = Cmat(n,D1_GAA) +
     &             one_eighth*one_over_rho32*fac
               endif
c
               if (rho(n,3).gt.tol_rho) then
                 Amat(n,2) =Amat(n,2) - 
c     &             one_quarter*temp_potential2*one_over_rho33*fac
     &             one_eighth*temp_energy2*one_over_rho33*fac
                 Cmat(n,D1_GBB) = Cmat(n,D1_GBB) +
     &             one_eighth*one_over_rho33*fac
               endif
c
            endif
   20    continue
c
         Ts = Ts + one_eighth*term*fac

      endif
c
      return
      end
c
cc AJL/End
