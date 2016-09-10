c-----------------------------------------------------------------------
c
      subroutine wfn1_allocate_drv(wfn1_param,wfn1_deriv)
      implicit none
#include "wfn1_param.fh"
#include "wfn1_drv.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
!>
!> \brief Allocate the derivative components
!>
!> Allocate global arrays for the natural orbitals, and the correlation
!> functions.
!>
      type(wfn1_prm),     intent(in)  :: wfn1_param
      type(wfn1_drv),     intent(out) :: wfn1_deriv
c
      integer nproc     ! the number of processor ranks
      integer iproc     ! the current processor rank
c
      character*18 pname
      pname = 'wfn1_allocate_drv:'
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
c     The natural orbital components
c
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,
     &                         'dNa',-1,-1,wfn1_deriv%g_dna))
     &  call errquit(pname//" could not allocate g_dna",0,MA_ERR)
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,
     &                         'dNb',-1,-1,wfn1_deriv%g_dnb))
     &  call errquit(pname//" could not allocate g_dnb",0,MA_ERR)
c
c     The correlation function components
c
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,
     &                         'dCa',-1,-1,wfn1_deriv%g_dca))
     &  call errquit(pname//" could not allocate g_dna",0,MA_ERR)
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,
     &                         'dCb',-1,-1,wfn1_deriv%g_dcb))
     &  call errquit(pname//" could not allocate g_dnb",0,MA_ERR)
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_deallocate_drv(wfn1_deriv)
      implicit none
#include "wfn1_drv.fh"
#include "global.fh"
#include "errquit.fh"
!>
!> \brief Deallocate wavefunction components
!>
      type(wfn1_drv), intent(inout) :: wfn1_deriv
      character*20 pname
      pname = "wfn1_deallocate_drv:"
!>
      if (.not.ga_destroy(wfn1_deriv%g_dna)) then
        call errquit(pname//" failed to dealloc g_dna",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_deriv%g_dnb)) then
        call errquit(pname//" failed to dealloc g_dnb",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_deriv%g_dca)) then
        call errquit(pname//" failed to dealloc g_dca",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_deriv%g_dcb)) then
        call errquit(pname//" failed to dealloc g_dcb",0,GA_ERR)
      endif
!>
      end
c
c-----------------------------------------------------------------------
