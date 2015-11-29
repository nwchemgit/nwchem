C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Convert a gradient to a search direction
C>
C> The gradient of the energy points into the direction along which
C> the energy increases. However, the energy expression of interest
C> is strictly speaking unbound. The direction we should be searching
C> in is down-hill in the electronic coordinates but up-hill in the
C> Lagrangian multipliers. Hence we need to apply a trivial
C> transformation to the gradient to generate the required search 
C> direction.
C>
      subroutine wfn1_nxt_gradient_to_search(nbf,grad)
      implicit none
c
#include "wfn1_nxt.fh"
c
      integer nbf !< [Input] The number of basis functions
      double precision grad(8*nbf*nbf) !< [In/Output] The gradient /
      !< search direction
c
      call dscal(nbf*nbf,-1.0d0,grad(wfn1_nxt_aoccl()),1)
      call dscal(nbf*nbf,-1.0d0,grad(wfn1_nxt_boccl()),1)
      call dscal(nbf*nbf,-1.0d0,grad(wfn1_nxt_anatl()),1)
      call dscal(nbf*nbf,-1.0d0,grad(wfn1_nxt_bnatl()),1)
      end
C>
C> @}
