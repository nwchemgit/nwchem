#include "gw_simd.fh"
      subroutine gw_analytic_addrpa(rpa,eia,mynpoles,lpole)
c
c     Adds squared eigenvalue differences to diagonal of
c     RPA matrix. Each rank has all columns and only a
c     subset of rows. 
c
c     This is the first term of (A-B)^{1/2}(A+B)(A-B)^{1/2} in
c     Equation (32), (34).

      implicit none
      integer,intent(in) :: mynpoles,lpole
      double precision,intent(in) :: eia(mynpoles)
      double precision,intent(out) :: rpa(mynpoles,*)

      integer ipole
!$omp parallel do MYSIMD
      do ipole=lpole,lpole+mynpoles-1
        rpa(ipole-lpole+1,ipole) = rpa(ipole-lpole+1,ipole) +
     $                             eia(ipole-lpole+1)**2   
      enddo
!$omp end parallel do MYSIMD

      end subroutine
