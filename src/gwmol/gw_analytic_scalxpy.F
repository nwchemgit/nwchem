      subroutine gw_analytic_scalxpy(XPY,omega,mynpoles,npoles)
c
c     Scale the eigenvector T with the inverse square root of the
c     neutral charge excitations (Equation (33)). 
c
c     Equation (33) of the manuscript has the 1/\Omega^{1/2} missing!
c     (see Equations (42a) and (42b) from Comp. Phys. Comm. 208, 149
c     (2016)).
c

      implicit none
      integer,intent(in) :: mynpoles, npoles
      double precision,intent(in) :: omega(npoles)
      double precision,intent(out) :: XPY(mynpoles,npoles)

      integer ipole,jpole

!$omp parallel do simd
      do ipole=1,npoles
        XPY(:,ipole) = XPY(:,ipole)/dsqrt(omega(ipole))
      enddo
!$omp end parallel do simd     

      end subroutine
