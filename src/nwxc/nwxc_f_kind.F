C> \ingroup nwxc
C> @{
C>
C> \brief Does this functional have any exchange-correlation terms?
C>
C> A function to query the functional for the presence of any 
C> exchange-correlation terms. With the current hybrid and double hybrid
C> functionals it has become possible to specify functionals that 
C> produce just pure Hartree-Fock or MP2 without any DFT 
C> exchange-correlation at all. Although unusual the code has to be able
C> to handle these cases correctly. Hence there is a need to be able to
C> query whether any exchange-correlation terms need to be evaluated.
C>
C> \return Return .true. if the functional has an exchange-correlation
C> term of any kind, and .false. otherwise.
C>
      logical function nwxc_has_xc()
      implicit none
      logical  nwxc_is_lda, nwxc_is_gga, nwxc_is_mgga
      external nwxc_is_lda, nwxc_is_gga, nwxc_is_mgga
      nwxc_has_xc = nwxc_is_lda() .or. nwxc_is_gga() .or. nwxc_is_mgga()
      return
      end
C>
C> \brief Is this functional an LDA functional?
C>
C> An LDA functional is a functional that depends on the local density
C> only. 
C>
C> The recommended way to query the nature of a functional (e.g. for
C> memory allocation purposes) is:
C> \code
C>    if (nwxc_is_mgga()) then
C>      ... do something ...
C>    else if (nwxc_is_gga()) then
C>      ... do something ...
C>    else if (nwxc_is_lda()) then
C>      ... do something ...
C>    else
C>      ... bomb ...
C>    endif
C> \endcode
C>
C> \return Return .true. if the functional depends on the density only,
C> and .false. otherwise.
C>
      logical function nwxc_is_lda()
      implicit none
#include "nwxcP.fh"
      nwxc_is_lda = nwxc_depends_on_rho .and. 
     +              (.not. nwxc_depends_on_gamma) .and.
     +              (.not. nwxc_depends_on_tau)
      return
      end
C>
C> \brief Is this functional a GGA functional?
C>
C> A GGA functional is a functional that depends on the local density
C> and the density gradient, but not on the kinetic energy density.
C>
C> The recommended way to query the nature of a functional (e.g. for
C> memory allocation purposes) is:
C> \code
C>    if (nwxc_is_mgga()) then
C>      ... do something ...
C>    else if (nwxc_is_gga()) then
C>      ... do something ...
C>    else if (nwxc_is_lda()) then
C>      ... do something ...
C>    else
C>      ... bomb ...
C>    endif
C> \endcode
C>
C> \return Return .true. if the functional depends on the density and
C> the density gradient, and .false. otherwise.
C>
      logical function nwxc_is_gga()
      implicit none
#include "nwxcP.fh"
      nwxc_is_gga = nwxc_depends_on_rho .and. 
     +              nwxc_depends_on_gamma .and.
     +              (.not. nwxc_depends_on_tau)
      return
      end
C>
C> \brief Is this functional a meta-GGA functional?
C>
C> A meta-GGA functional is a functional that depends on the local
C> density, the density gradient, and on the kinetic energy density.
C>
C> The recommended way to query the nature of a functional (e.g. for
C> memory allocation purposes) is:
C> \code
C>    if (nwxc_is_mgga()) then
C>      ... do something ...
C>    else if (nwxc_is_gga()) then
C>      ... do something ...
C>    else if (nwxc_is_lda()) then
C>      ... do something ...
C>    else
C>      ... bomb ...
C>    endif
C> \endcode
C>
C> \return Return .true. if the functional depends on the density,
C> the density gradient, and the kinetic energy density, and .false.
C> otherwise.
C>
      logical function nwxc_is_mgga()
      implicit none
#include "nwxcP.fh"
      nwxc_is_mgga = nwxc_depends_on_rho .and. 
     +               nwxc_depends_on_gamma .and.
     +               nwxc_depends_on_tau
      return
      end
C> @}
