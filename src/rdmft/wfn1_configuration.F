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
        wfn1_wave%dcnta(ii) = 1.0d0
      enddo
      if (wfn1_param%nea.gt.0) then
        wfn1_wave%dcnta(wfn1_param%nea) = wfn1_param%fraca
      endif
c
      do ii = 1, wfn1_param%neb
        wfn1_wave%icnfb(ii) = ii
        wfn1_wave%dcntb(ii) = 1.0d0
      enddo
      if (wfn1_param%neb.gt.0) then
        wfn1_wave%dcntb(wfn1_param%neb) = wfn1_param%fracb
      endif
c
      end
c
c-----------------------------------------------------------------------
