c-----------------------------------------------------------------------
c
      subroutine wfn1_configuration(wfn1_param,wfn1_wave)
      implicit none
C>
C> \brief Set the electron configuration up
C>
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
      type(wfn1_prm), intent(in)    :: wfn1_param
      type(wfn1_wfn), intent(inout) :: wfn1_wave
c
      integer ii ! counter
c
      do ii = 1, wfn1_param%nea
        wfn1_wave%icnfa(ii) = ii
      enddo
c
      do ii = 1, wfn1_param%neb
        wfn1_wave%icnfb(ii) = ii
      enddo
c
      end
