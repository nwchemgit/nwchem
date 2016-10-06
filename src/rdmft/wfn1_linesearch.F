c-----------------------------------------------------------------------
c
      subroutine wfn1_linesearch(wfn1_param,wfn1_inst,wfn1_mat,
     &                           wfn1_wave,wfn1_deriv,wfn1_energy,
     &                           wfn1_trial_wv,wfn1_trial_drv,
     &                           wfn1_trial_en,
     &                           g_bfmo,g_momo1,g_momo2)
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
      type(wfn1_en),      intent(inout) :: wfn1_energy
      type(wfn1_wfn),     intent(inout) :: wfn1_trial_wv
      type(wfn1_drv),     intent(inout) :: wfn1_trial_drv
      type(wfn1_en),      intent(inout) :: wfn1_trial_en
c
      integer, intent(inout) :: g_bfmo  ! nbf x nmo scratch GA
      integer, intent(inout) :: g_momo1 ! nmo x nmo scratch GA
      integer, intent(inout) :: g_momo2 ! nmo x nmo scratch GA
c
      integer, external :: wfn1_idmin
c
      integer             npts        !< The number of points
      integer             mxpts       !< The maximum number of points
      parameter          (mxpts = 10) !< Set the number points
      double precision xx(mxpts)      !< The x-coordinates
      double precision fx(mxpts)      !< The function values
      double precision df(mxpts)      !< The function gradient
      double precision tr             !< The trust radius
      double precision aa, bb, cc     !< The factor of the parabola
      double precision xm, fxm, dfm
c
      integer indxmin
      integer ix1, ix2, ix3
c
      tr = 2.0d0
c
c     Store the initial point
c
      npts = 1
      xx(npts) = 0.0d0
      fx(npts) = wfn1_energy%E_tot
      df(npts) = wfn1_trace(wfn1_deriv,wfn1_deriv)
c
c     Take one step
c
      npts = npts + 1
      xx(npts) = 1.0d0
      call wfn1_gen_trial_wfn(xx(npts),wfn1_param,wfn1_wave,wfn1_deriv,
     &                        wfn1_trial_wv,g_momo1,g_momo2)
      call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &              wfn1_trial_drv,wfn1_trial_en,g_bfmo,g_momo1)
      fx(npts) = wfn1_trial_en%E_tot
      df(npts) = wfn1_trace(wfn1_deriv,wfn1_trial_drv)
      call wfn1_f0df0f1(xx(2),fx(1),df(1),fx(2),tr,aa,bb,cc,xm,fxm)
c
c     Always select the next step (hoping that using the quadratic
c     approximation will pay off)
c
      call wfn1_gen_trial_wfn(xm,wfn1_param,wfn1_wave,wfn1_deriv,
     &                        wfn1_trial_wv,g_momo1,g_momo2)
      call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &              wfn1_trial_drv,wfn1_trial_en,g_bfmo,g_momo1)
      dfm = wfn1_trace(wfn1_deriv,wfn1_trial_drv)
      call wfn1_insert(mxpts,npts,xx,fx,df,xm,wfn1_trial_en%E_tot,dfm)
c
      do while (wfn1_trial_en%E_tot > wfn1_energy%E_tot)
        indxmin = wfn1_idmin(npts,fx,1)
        call wfn1_closest_three(npts,xx,indxmin,ix1,ix2,ix3)
        call wfn1_f1f2f3(xx(ix1),xx(ix2),xx(ix3),
     &                   fx(ix1),fx(ix2),fx(ix3),
     &                   tr,aa,bb,cc,xm,fxm)
        call wfn1_gen_trial_wfn(xm,wfn1_param,wfn1_wave,wfn1_deriv,
     &                          wfn1_trial_wv,g_momo1,g_momo2)
        call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &                wfn1_trial_drv,wfn1_trial_en,g_bfmo,g_momo1)
        dfm = wfn1_trace(wfn1_deriv,wfn1_trial_drv)
        call wfn1_insert(mxpts,npts,xx,fx,df,xm,wfn1_trial_en%E_tot,dfm)
      enddo
c
      end subroutine wfn1_linesearch
c
c-----------------------------------------------------------------------
