c-----------------------------------------------------------------------
c
      double precision function wfn1_dnrm2(wfn1_deriv)
      implicit none
C>
C> \brief Compute the norm of the gradients
C>
#include "wfn1_drv.fh"
c
      type(wfn1_drv), intent(in) :: wfn1_deriv
c
      double precision :: dnrm2
      double precision :: term2
c
      dnrm2 = 0.0d0
      call ga_normf(wfn1_deriv%g_dna,term2)
      dnrm2 = dnrm2 + term2**2
      call ga_normf(wfn1_deriv%g_dnb,term2)
      dnrm2 = dnrm2 + term2**2
      call ga_normf(wfn1_deriv%g_dca,term2)
      dnrm2 = dnrm2 + term2**2
      call ga_normf(wfn1_deriv%g_dcb,term2)
      dnrm2 = dnrm2 + term2**2
      dnrm2 = sqrt(dnrm2)
      wfn1_dnrm2 = dnrm2
c
      return
      end function wfn1_dnrm2
c
c-----------------------------------------------------------------------
