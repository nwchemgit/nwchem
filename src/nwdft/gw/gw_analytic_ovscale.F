      subroutine gw_analytic_ovscale(id,OV,eia,nri,npoles)
      implicit none
      integer,intent(in)            :: id, nri, npoles
      double precision, intent(in)  :: eia(npoles)
      double precision, intent(out) :: OV(nri,npoles)

      integer ipole

      if (id.eq.1) then
        do ipole=1,npoles
          OV(1:nri,ipole) = OV(1:nri,ipole)*dsqrt(eia(ipole))
        enddo
      else
        do ipole=1,npoles
          OV(1:nri,ipole) = OV(1:nri,ipole)/dsqrt(eia(ipole))
        enddo
      endif

      end subroutine
