      subroutine tddft_grad_fix_exc(nq,exc,wght)
      implicit none
c
c     $Id$
c
c     This subroutine tweaks the energy expression in exc for the 
c     derivatives of the weights contribution. The situation is that
c     for the derivatives of the excitation energy there will be a 
c     contribution involving the derivatives of the weights. However
c     as the weights themselves get multiplied into the partial 
c     derivatives of the functional they are incorporated also in the
c     excitation energy expression. This subroutine devides them out
c     again so that subroutine dftg_wderiv can generate the proper 
c     derivatives of the weights contribution.
c
c     Written by Huub van Dam, December 2006.
c
c     Input:
c    
      integer nq  ! the number of grid points in the current batch
c
      double precision wght(nq) ! the quadrature weights
c
c     Input/Output:
c
      double precision exc(nq)  ! the excitation energy 
      ! on input : including the quadrature weight
      ! on output: excluding the quadrature weight
c
c     Local:
c
      integer iq ! counter over grid points
c
c     Code:
c
      do iq = 1, nq
        if (dabs(wght(iq)).gt.0.0d0) then
          exc(iq) = exc(iq)/wght(iq)
        endif
      enddo
c
      end
