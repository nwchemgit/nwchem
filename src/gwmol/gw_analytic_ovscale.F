      subroutine gw_analytic_ovscale(id,OV,eia,nri,npoles)
c
c     Scale Occ-Vir ERI block with the root of the eigenvalue
c     differences. This scaling will produce the right 
c     matrix structure for the (A-B)^{1/2}(A+B)(A-B)^{1/2}
c     matrix.
c
      implicit none
      integer,intent(in)            :: id, nri, npoles
      double precision, intent(in)  :: eia(npoles)
      double precision, intent(out) :: OV(nri,npoles)

      integer ipole,iMaxThreads,iri

      !
      ! Scale the elements
      !
      if (id.eq.1) then

!$omp   parallel do simd collapse(2)
        do ipole=1,npoles
          do iri=1,nri
            OV(iri,ipole) = OV(iri,ipole)*dsqrt(eia(ipole))
          enddo
        enddo
!$omp   end parallel do simd

      !
      ! Get the original ERIs
      !
      else

!$omp   parallel do simd collapse(2)
        do ipole=1,npoles
          do iri=1,nri
            OV(iri,ipole) = OV(iri,ipole)/dsqrt(eia(ipole))
          enddo
        enddo
!$omp   end parallel do simd

      endif

      end subroutine
