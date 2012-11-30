C> \ingroup nwxc
C> @{
C>
C> \brief Retrieve the Coulomb Attenuation parameters
C>
C> Using the Coulomb Attenuation or range separation the Coulomb
C> operator is written as:
C> \f{eqnarray*}{
C>   \frac{1}{r_{12}} &=& \frac{1-[\alpha+\beta\mathrm{erf}(\gamma r_{12})]}{r_{12}} 
C>                     +  \frac{\alpha+\beta\mathrm{erf}(\gamma r_{12})}{r_{12}}
C> \f}
C> This routine retrieves the parameters \f$ \alpha \f$, \f$ \beta \f$,
C> and \f$ \gamma \f$. To enable returning multiple values at once this
C> routine is implemented as a subroutine.
C>
      subroutine nwxc_get_cam(alpha,beta,gamma)
c     implicit none
#include "nwxcP.fh"
      double precision alpha !< [Output] The alpha parameter
      double precision beta  !< [Output] The beta parameter
      double precision gamma !< [Output] The gamma parameter
c     alpha = nwxc_cam_alpha
c     beta  = nwxc_cam_beta
c     gamma = nwxc_cam_gamma
      end
C> @}
