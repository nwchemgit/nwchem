      subroutine gw_analytic_scalxpy(XPY,omega,mynpoles,npoles)
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
