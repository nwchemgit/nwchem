      subroutine gw_analytic_ovscale(id,OV,eia,nri,npoles)
      implicit none
      integer,intent(in)            :: id, nri, npoles
      double precision, intent(in)  :: eia(npoles)
      double precision, intent(out) :: OV(nri,npoles)

      integer ipole,iMaxThreads,iri

      if (id.eq.1) then

!$omp   parallel do collapse(2)
        do ipole=1,npoles
          do iri=1,nri
            OV(iri,ipole) = OV(iri,ipole)*dsqrt(eia(ipole))
          enddo
        enddo
!$omp   end parallel do

      else

!$omp   parallel do collapse(2)
        do ipole=1,npoles
          do iri=1,nri
            OV(iri,ipole) = OV(iri,ipole)/dsqrt(eia(ipole))
          enddo
        enddo
!$omp   end parallel do

      endif

      end subroutine
