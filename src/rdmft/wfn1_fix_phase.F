c-----------------------------------------------------------------------
c
      subroutine wfn1_fix_phase(wfn1_wave)
      implicit none
C>
C> \brief Fix the phase of the orbitals
C>
#include "wfn1_wfn.fh"
c
      type(wfn1_wfn), intent(inout) :: wfn1_wave
c
      call movecs_fix_phase(wfn1_wave%g_na)
      call movecs_fix_phase(wfn1_wave%g_nb)
      call movecs_fix_phase(wfn1_wave%g_ca)
      call movecs_fix_phase(wfn1_wave%g_cb)
c
      end
c
c-----------------------------------------------------------------------
