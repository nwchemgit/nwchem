C> \ingroup nwxc
C> @{
C>
C> \brief Does this functional have Hartree-Fock exchange?
C>
C> A function to query the functional for the presence of a
C> Hartree-Fock exchange term. This function just queries the 
C> functional definition. In practice the Hartree-Fock exchange term
C> must be implemented in the application code as this library cannot
C> provide that.
C>
C> \return Returns .true. if the functional contains a Hartree-Fock 
C> exchange term (i.e. is a hybrid functional), and .false. otherwise.
C>
      logical function nwxc_has_hfx()
      implicit none
#include "nwxcP.fh"
      nwxc_has_hfx = nwxc_wght_hf .gt. 0.0d0
      return
      end
C>
C> \brief Does this functional have MP2?
C>
C> A function to query the functional for the presence of an MP2
C> correlation term. This function just queries the 
C> functional definition. In practice the MP2 correlation term
C> must be implemented in the application code as this library cannot
C> provide that.
C>
C> \return Returns .true. if the functional contains an MP2 correlation
C> term (i.e. is a double hybrid functional), and .false. otherwise.
C>
      logical function nwxc_has_mp2c()
      implicit none
#include "nwxcP.fh"
      nwxc_has_mp2c = nwxc_wght_mp2 .gt. 0.0d0
      return
      end
C>
C> \brief Does this functional have Coulomb Attenuation?
C>
C> A function to query the functional for the use Coulomb Attentuated
C> operators. With Coulomb Attentuation part of the expression is 
C> incorporate in the Density Functional component and part needs to
C> be done by the application code.
C>
C> \return Returns .true. if the functional use Coulomb Attenuation,
C> and .false. otherwise.
C>
      logical function nwxc_has_cam()
      implicit none
#include "nwxcP.fh"
      nwxc_has_cam = nwxc_cam_gamma .ne. 0.0d0 .or.
     +               nwxc_cam_beta  .gt. 0.0d0 .or.
     +               nwxc_cam_alpha .gt. 0.0d0
      return
      end
C>
C> \brief Does this functional have dispersion corrections?
C>
C> A function to query the functional for the presence of a
C> dispersion correction term. This function just queries the 
C> functional definition. In practice the dispersion correction term
C> must be implemented in the application code as this library cannot
C> provide that.
C>
C> \return Returns .true. if the functional contains a dispersion
C> correction term (i.e. is a hybrid functional), and .false. otherwise.
C>
      logical function nwxc_has_disp()
      implicit none
#include "nwxcP.fh"
      call errquit("not sure how to do this yet",0,0)
      nwxc_has_disp = nwxc_wght_hf .gt. 0.0d0
      return
      end
C> @}
