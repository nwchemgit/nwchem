c-----------------------------------------------------------------------
c
      subroutine wfn1_print_analysis(wfn1_param,wfn1_inst,wfn1_wave)
      implicit none
C>
C> \brief Print the wavefunction analysis
C>
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
c
      type(wfn1_prm),     intent(in) :: wfn1_param
      type(wfn1_prminst), intent(in) :: wfn1_inst
      type(wfn1_wfn),     intent(in) :: wfn1_wave
c
      double precision, allocatable :: e1a(:), e1b(:)
      double precision, allocatable :: occa(:), occb(:)
c
      double precision :: tanalyze
c
      integer, allocatable :: irrep(:)
c
      integer :: ilo, ihi
c
      character*255 blob
c
      tanalyze = 0.025d0
      ilo = 1
      ihi = wfn1_param%nmo
c
      allocate(e1a(1:wfn1_param%nmo),e1b(1:wfn1_param%nmo))
      allocate(occa(1:wfn1_param%nmo),occb(1:wfn1_param%nmo))
      allocate(irrep(1:wfn1_param%nmo))
      e1a   = 0.0d0
      e1b   = 0.0d0
      irrep = 0
c
      call wfn1_calc_occ(wfn1_wave%g_ca,wfn1_param%nmo,wfn1_param%nea,
     &                   wfn1_wave%icnfa,occa)
      call wfn1_calc_occ(wfn1_wave%g_cb,wfn1_param%nmo,wfn1_param%neb,
     &                   wfn1_wave%icnfb,occb)
c
      blob = "WFN1 Final Alpha Natural Orbital Analysis"
      call movecs_print_anal(wfn1_inst%basis,ilo,ihi,tanalyze,
     &                       wfn1_wave%g_na,blob,.true.,e1a,.false.,
     &                       irrep,.true.,occa)
c
      blob = "WFN1 Final Beta Natural Orbital Analysis"
      call movecs_print_anal(wfn1_inst%basis,ilo,ihi,tanalyze,
     &                       wfn1_wave%g_nb,blob,.true.,e1b,.false.,
     &                       irrep,.true.,occb)
c
      deallocate(occa,occb)
      deallocate(e1a,e1b)
c
      end
c
c-----------------------------------------------------------------------
