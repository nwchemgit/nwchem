      subroutine gw_analytic_scalxpy(XPY,omega,mynpoles,npoles)
      implicit none
      integer,intent(in) :: mynpoles, npoles
      double precision,intent(in) :: omega(npoles)
      double precision,intent(out) :: XPY(mynpoles,npoles)

      integer ipole,jpole

!$omp parallel do collapse(2)
      do ipole=1,npoles
        do jpole=1,mynpoles
          XPY(jpole,ipole) = XPY(jpole,ipole)/dsqrt(omega(ipole))
        end do
      enddo
!$omp end parallel do      

      end subroutine
