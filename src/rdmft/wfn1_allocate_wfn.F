c-----------------------------------------------------------------------
c
      subroutine wfn1_allocate_wfn(wfn1_param,wfn1_inst,wfn1_wave)
      implicit none
#include "wfn1_param.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_wfn.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
!>
!> \brief Allocate the wavefunction components
!>
!> Allocate global arrays for the natural orbitals, the correlation
!> functions, as well as arrays to store the electron configuration.
!>
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_wfn),     intent(inout) :: wfn1_wave
c
      integer nblock(2) ! the number of blocks in each dimension
      integer nproc     ! the number of processor ranks
      integer iproc     ! the current processor rank
c
      integer ii        ! counter
      integer nfnc      ! number of functions
c
      integer, allocatable :: map(:,:)
c
      character*18 pname
      pname = 'wfn1_allocate_wfn:'
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      allocate(map(nproc,2))
c
      call util_atom_blocked(nproc,wfn1_inst%geom,wfn1_inst%basis,
     &                       map,nblock)
c
c     dimension 2 does not have to be on atom boundaries so remap the
c     blocking of it and then allocate the Natural orbital arrays
c
      nfnc = (wfn1_param%nmo+nblock(2)-1)/nblock(2)
      do ii = 1, nblock(2)
        map(ii,2) = 1+(ii-1)*nfnc
      enddo
c
      if (.not.ga_create_irreg(MT_DBL,wfn1_param%nbf,wfn1_param%nmo,
     &                         'Na',map(:,1),nblock(1),
     &                         map(:,2),nblock(2),wfn1_wave%g_na))
     &  call errquit(pname//" could not allocate g_na",0,MA_ERR)
      if (.not.ga_create_irreg(MT_DBL,wfn1_param%nbf,wfn1_param%nmo,
     &                         'Nb',map(:,1),nblock(1),
     &                         map(:,2),nblock(2),wfn1_wave%g_nb))
     &  call errquit(pname//" could not allocate g_nb",0,MA_ERR)
c
c     Now remap dimension 1 and allocate the Correlation function arrays
c
      nfnc = (wfn1_param%nmo+nblock(1)-1)/nblock(1)
      do ii = 1, nblock(1)
        map(ii,1) = 1+(ii-1)*nfnc
      enddo
c
      if (.not.ga_create_irreg(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,
     &                         'Ca',map(:,2),nblock(2),
     &                         map(:,1),nblock(1),wfn1_wave%g_ca))
     &  call errquit(pname//" could not allocate g_na",0,MA_ERR)
      if (.not.ga_create_irreg(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,
     &                         'Cb',map(:,2),nblock(2),
     &                         map(:,1),nblock(1),wfn1_wave%g_cb))
     &  call errquit(pname//" could not allocate g_nb",0,MA_ERR)
c
      deallocate(map)
c
      allocate(wfn1_wave%dcnta(1:wfn1_param%nea))
      allocate(wfn1_wave%dcntb(1:wfn1_param%neb))
      allocate(wfn1_wave%icnfa(1:wfn1_param%nea))
      allocate(wfn1_wave%icnfb(1:wfn1_param%neb))
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_deallocate_wfn(wfn1_wave)
      implicit none
#include "wfn1_wfn.fh"
#include "global.fh"
#include "errquit.fh"
!>
!> \brief Deallocate wavefunction components
!>
      type(wfn1_wfn), intent(inout) :: wfn1_wave
      character*20 pname
      pname = "wfn1_deallocate_wfn:"
!>
      deallocate(wfn1_wave%icnfa)
      deallocate(wfn1_wave%icnfb)
      deallocate(wfn1_wave%dcnta)
      deallocate(wfn1_wave%dcntb)
      if (.not.ga_destroy(wfn1_wave%g_na)) then
        call errquit(pname//" failed to dealloc g_na",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_wave%g_nb)) then
        call errquit(pname//" failed to dealloc g_nb",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_wave%g_ca)) then
        call errquit(pname//" failed to dealloc g_ca",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_wave%g_cb)) then
        call errquit(pname//" failed to dealloc g_cb",0,GA_ERR)
      endif
!>
      end
c
c-----------------------------------------------------------------------
