C> \ingroup nwxc
C> @{
C>
C> \brief Retrieve the Hartree-Fock exchange weight
C>
C> Simply look up and return the weight for the Hartree-Fock exchange.
C>
C> \return Returns the weight of the Hartree-Fock exchange term
C>
      double precision function nwxc_wght_hfx()
      implicit none
#include "nwxcP.fh"
      nwxc_wght_hfx = nwxc_wght_hf
      return
      end
C>
C> \brief Retrieve the MP2 correlation weight
C>
C> Simply look up and return the weight for the MP2 correlation term.
C>
C> \return Returns the weight of the MP2 correlation term
C>
      double precision function nwxc_wght_mp2c()
      implicit none
#include "nwxcP.fh"
      nwxc_wght_mp2c = nwxc_wght_mp2
      return
      end
C> @}
