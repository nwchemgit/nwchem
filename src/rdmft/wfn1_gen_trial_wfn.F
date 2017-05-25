c-----------------------------------------------------------------------
c
      subroutine wfn1_gen_trial_wfn(alpha,wfn1_param,wfn1_wave,
     &                              wfn1_deriv,wfn1_trial_wv,g_t,g_exp)
C>
C> \brief Generate a trial wave function
C>
C> Given the current wave function, the derivatives of the energy
C> wrt. the wave function parameters, and a step size along the
C> gradient, construct the wave function obtained by applying the step.
C>
      implicit none
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh" 
#include "wfn1_param.fh" 
c
      double precision, intent(in)  :: alpha ! the step size
c
      type(wfn1_wfn), intent(in)    :: wfn1_wave     ! the current wave
                                                     ! function
      type(wfn1_wfn), intent(inout) :: wfn1_trial_wv ! the trial wave
                                                     ! function
      type(wfn1_drv), intent(in)    :: wfn1_deriv    ! the derivative
      type(wfn1_prm), intent(in)    :: wfn1_param    ! the calculation 
                                                     ! parameters
c
      integer, intent(inout) :: g_t   ! nmo x nmo scratch global array
      integer, intent(inout) :: g_exp ! nmo x nmo scratch global array
c
      integer :: nmo
c
      nmo = wfn1_param%nmo
c     call wfn1_tst_exp(alpha,wfn1_deriv%g_dna,g_t,g_exp)
      call wfn1_exp(alpha,wfn1_deriv%g_dna,g_t,g_exp)
      call ga_dgemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_na,g_exp,
     &              0.0d0,wfn1_trial_wv%g_na)
c
c     call wfn1_tst_exp(alpha,wfn1_deriv%g_dnb,g_t,g_exp)
      call wfn1_exp(alpha,wfn1_deriv%g_dnb,g_t,g_exp)
      call ga_dgemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_nb,g_exp,
     &              0.0d0,wfn1_trial_wv%g_nb)
c
c     call wfn1_tst_exp(alpha,wfn1_deriv%g_dca,g_t,g_exp)
      call wfn1_exp(alpha,wfn1_deriv%g_dca,g_t,g_exp)
      call ga_dgemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_ca,g_exp,
     &              0.0d0,wfn1_trial_wv%g_ca)
c
c     call wfn1_tst_exp(alpha,wfn1_deriv%g_dcb,g_t,g_exp)
      call wfn1_exp(alpha,wfn1_deriv%g_dcb,g_t,g_exp)
      call ga_dgemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_cb,g_exp,
     &              0.0d0,wfn1_trial_wv%g_cb)

c
      end
c
c-----------------------------------------------------------------------
