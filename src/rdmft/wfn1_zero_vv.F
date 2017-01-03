c-----------------------------------------------------------------------
c
      subroutine wfn1_zero_vv(wfn1_param,wfn1_wave,wfn1_deriv)
      implicit none
C>
C> \brief Zero the virtual-virtual block of the correlation function
C> matrices
C>
C> Rotations amongst the virtual-correlation functions are not going
C> to affect the total energy. Hence the energy based optimization
C> scheme will not provide an effective way to reduce the 
C> virtual-virtual block to zero. This in turn generates a spurious 
C> impression that the calculation is not converging. 
C> To avoid this confusion this routine explicitly zeroes the 
C> virtual-virtual blocks.
C>
#include "global.fh"
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
c
      type(wfn1_prm), intent(in)    :: wfn1_param
      type(wfn1_wfn), intent(in)    :: wfn1_wave
      type(wfn1_drv), intent(inout) :: wfn1_deriv
c
      integer :: nmo   ! local copy of nmo
c
      integer :: nproc ! the number of processors
      integer :: iproc ! the current processor rank
c
      integer :: ii    ! counter
      integer :: jj    ! counter
      integer :: ilo   ! patch limit
      integer :: ihi   ! patch limit
      integer :: jlo   ! patch limit
      integer :: jhi   ! patch limit
      integer :: ni    ! size of patch dimension i
c
      logical, allocatable :: ovirt(:)
c
      double precision, allocatable :: buf(:,:)
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      nmo = wfn1_param%nmo
      allocate(ovirt(1:nmo))
      call ga_sync()
c
c     Do the alpha-correlation-function matrix
c
      ovirt = .true.
      do ii = 1, wfn1_param%nea
        ovirt(wfn1_wave%icnfa(ii)) = .false.
      enddo
      call ga_distribution(wfn1_deriv%g_dca,iproc,ilo,ihi,jlo,jhi)
      ni = ihi-ilo+1
      allocate(buf(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_deriv%g_dca,ilo,ihi,jlo,jhi,buf,ni)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if (ovirt(ii).and.ovirt(jj)) buf(ii,jj) = 0.0d0
        enddo
      enddo
      call ga_put(wfn1_deriv%g_dca,ilo,ihi,jlo,jhi,buf,ni)
      deallocate(buf)
c
c     Do the beta-correlation-function matrix
c
      ovirt = .true.
      do ii = 1, wfn1_param%neb
        ovirt(wfn1_wave%icnfb(ii)) = .false.
      enddo
      call ga_distribution(wfn1_deriv%g_dcb,iproc,ilo,ihi,jlo,jhi)
      ni = ihi-ilo+1
      allocate(buf(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_deriv%g_dcb,ilo,ihi,jlo,jhi,buf,ni)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if (ovirt(ii).and.ovirt(jj)) buf(ii,jj) = 0.0d0
        enddo
      enddo
      call ga_put(wfn1_deriv%g_dcb,ilo,ihi,jlo,jhi,buf,ni)
      deallocate(buf)
c
      deallocate(ovirt)
      call ga_sync()
c
      end
c
c-----------------------------------------------------------------------
