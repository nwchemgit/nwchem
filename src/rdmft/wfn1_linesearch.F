c-----------------------------------------------------------------------
c
      subroutine wfn1_linesearch(wfn1_param,wfn1_inst,wfn1_mat,
     &                           wfn1_wave,wfn1_deriv,wfn1_energy,
     &                           wfn1_trial_wv,wfn1_trial_drv,
     &                           g_bfmo,g_momo)
      implicit none
C>
C> \brief Perform a linesearch for the minimum total energy
C>
C> Given the gradient of the energy wrt. the wave function parameters
C> this gradient is used as the search direction to find the minimum
C> of the total energy. The basic approach is to start from the
C> energy and gradient at the initial, and the energy after taking a
C> step, and fit a parabola through these three data points. Looking
C> for the minimum of that parabola a suggestion for the minimum is
C> obtained. In subsequent steps the energy at three points might be
C> used to fit a parabola. 
C>
C> When enough progress has been made the subroutine returns the 
C> current energies, the current wavefunction, and the current gradients
C> for the next linesearch (or the final solution if the calculation
C> is converged).
C>
#include "global.fh"
#include "wfn1_trace.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "wfn1_en.fh"
c
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
      type(wfn1_drv),     intent(inout) :: wfn1_deriv
      type(wfn1_wfn),     intent(inout) :: wfn1_trial_wv
      type(wfn1_drv),     intent(inout) :: wfn1_trial_drv
      type(wfn1_en),      intent(inout) :: wfn1_energy
c
      integer, intent(inout) :: g_bfmo ! nbf x nmo scratch GA
      integer, intent(inout) :: g_momo ! nmo x nmo scratch GA
c
      integer             npts        !< The number of points
      integer             mxpts       !< The maximum number of points
      parameter          (mxpts = 10) !< Set the number points
      double precision xx(mxpts)      !< The x-coordinates
      double precision fx(mxpts)      !< The function values
      double precision df(mxpts)      !< The function gradient
c
      npts = 1
      xx(npts) = 0.0d0
      fx(npts) = wfn1_energy%E_tot
      df(npts) = wfn1_trace(wfn1_deriv,wfn1_deriv)
c
      end subroutine wfn1_linesearch
c
c-----------------------------------------------------------------------
