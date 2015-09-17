cc AJL/Begin/FDE
c
      Subroutine ts_vw(tol_rho, fac, lfac, nlfac, rho, delrho, del2rho, 
     &                 Amat, Cmat, nq, ipol, Ts, qwght, ldew, func)
c
      Implicit none
#include "errquit.fh"
c
#include "stdio.fh"
c      
c#include "dft2drv.fh"
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
c Need to remove this by altering derivative
c
      double precision del2rho(nq,ipol*(ipol+1)/2)
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
      double precision temp_potential1, temp_potential2

      double precision one_over_rho32, one_over_rho33
c
      term = 0.d0
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c AJL/Begin - only this case applies as kinetic energy is 1e
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
c               temp_potential = 2*max(del2rho(n,1),tol_rho)
c I think this is just in case of numerical instability also!
              temp_potential = 2*del2rho(n,1)
              temp_potential = temp_energy - temp_potential
c \frac{\partial Tw[p]}{\partial \rho}
c = \frac{\partial tw}{\partial \rho} -
c   \triangledown . \frac{\partial tw}{\partial \triangledown \rho}
              Amat(n,1) =Amat(n,1) + 
     &               one_eighth*temp_potential*one_over_rho*fac
ccc Updated version:
cc \frac{\partial tw}{\partial \rho} = 
c              Amat(n,1) = Amat(n,1) - 
c     &               one_eighth*temp_energy*one_over_rho*fac
cc \frac{\partial tw}{\partial |\triangledown \rho|^{2}} = \frac{\partial tw}{\partial\gamma^{\alpha\alpha}}
cc                                                       = \frac{1}{8 \rho} 
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) +
c     &        one_eighth*one_over_rho*fac
cc Something is wrong with this as it won't give the same answers
cc as the "brute-force" implementation depending on del2rho....
cc I get quite close results if I swap 1/8 for 1/4 above, but not good
cc enough!
            endif
   10    continue
c 
         Ts = Ts + one_eighth*term*fac

      else
c
c        ======> SPIN-UNRESTRICTED <======
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

               temp_potential1 = 0.5d0*temp_energy1
               temp_potential2 = 0.5d0*temp_energy2
c               
               temp_potential1 = temp_potential1-del2rho(n,1)     
               temp_potential2 = temp_potential2-del2rho(n,2)

               if (rho(n,2).gt.tol_rho) 
     & Amat(n,1) =Amat(n,1) + 
     &           one_quarter*temp_potential1*one_over_rho32*fac
               if (rho(n,3).gt.tol_rho) 
     & Amat(n,2) =Amat(n,2) + 
     &           one_quarter*temp_potential2*one_over_rho33*fac
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
