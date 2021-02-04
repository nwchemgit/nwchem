      subroutine gw_analytic_addrpa(rpa,eia,mynpoles,lpole)
      implicit none
      integer,intent(in) :: mynpoles,lpole
      double precision,intent(in) :: eia(mynpoles)
      double precision,intent(out) :: rpa(mynpoles,*)

      integer ipole

      !$omp parallel do
      do ipole=lpole,lpole+mynpoles-1
        rpa(ipole-lpole+1,ipole) = rpa(ipole-lpole+1,ipole) +
     $                             eia(ipole-lpole+1)**2   
      enddo
      !$omp end parallel do

      end subroutine
