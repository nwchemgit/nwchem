c-----------------------------------------------------------------------
c
      subroutine wfn1_e_g_2el(wfn1_param,wfn1_wave,wfn1_deriv,
     &                        v_dca,v_dcb,g_2aa,g_2bb,
     &                        g_waa,g_wab,g_wba,g_wbb,
     &                        g_saa,g_sab,g_sba,g_sbb,
     &                        g_bfmo,g_momo,E_2el)
      implicit none
C>
C> \brief Calculate the 2-electron terms
C>
C> In the algorithm we use the fact that the 2-electron integrals 
C> are the most expensive part of the calculation. Hence the factor
C> \f$W\f$ and its derivatives should be calculated only once.
C> As the less expensive factor the entropy and its derivatives can be
C> evaluated more than once. 
C>
C> The above argument suggests the following algorithm:
C> 1. Calculate \f$S_{rs}\f$
C> 2. Calculate the \f$W_{rs}(1-S_{rs})\f$ terms including the 
C>    derivatives of \f$W_{rs}\f$ wrt. the wave function parameters.
C> 3. Calculate the derivatives of \f$S_{rs}\f$ wrt. the wave
C>    function parameters.
C>
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      type(wfn1_prm), intent(in)    :: wfn1_param  ! the job parameters
      type(wfn1_wfn), intent(in)    :: wfn1_wave   ! the wavefunction
      type(wfn1_drv), intent(inout) :: wfn1_deriv  ! the derivative
c
      double precision, intent(inout) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(inout) :: v_dcb(1:wfn1_param%nmo)
c
      integer, intent(inout) :: g_2aa  ! nmo x nmo matrix for aa-block
                                       ! of the 2-electron density
                                       ! matrix diagonal
      integer, intent(inout) :: g_2bb  ! nmo x nmo matrix for bb-block
                                       ! of the 2-electron density
                                       ! matrix diagonal
      integer, intent(inout) :: g_waa  ! nea x nea matrix for aa-block
                                       ! of the 2-electron interaction
      integer, intent(inout) :: g_wab  ! nea x neb matrix for ab-block
                                       ! of the 2-electron interaction
      integer, intent(inout) :: g_wba  ! neb x nea matrix for ba-block
                                       ! of the 2-electron interaction
      integer, intent(inout) :: g_wbb  ! neb x nea matrix for ba-block
                                       ! of the 2-electron interaction
      integer, intent(inout) :: g_saa  ! nea x nea matrix for aa-block
                                       ! of the pair entropy
      integer, intent(inout) :: g_sab  ! nea x neb matrix for ab-block
                                       ! of the pair entropy
      integer, intent(inout) :: g_sba  ! neb x nea matrix for ba-block
                                       ! of the pair entropy
      integer, intent(inout) :: g_sbb  ! neb x nea matrix for ba-block
                                       ! of the pair entropy
      integer, intent(inout) :: g_bfmo ! nbf x nmo matrix
      integer, intent(inout) :: g_momo ! nmo x nmo matrix
c
      double precision, intent(out) :: E_2el ! the two-electron energy
c
c     Local
c
      double precision, allocatable :: occa(:) ! a-occupation numbers
      double precision, allocatable :: occb(:) ! b-occupation numbers
c
      integer :: nproc
      integer :: iproc
c
      integer :: nmo            ! number of MOs
      integer :: nea            ! number of alpha electrons
      integer :: neb            ! number of beta electrons
      integer :: ir, is, it     ! correlation function labels
      integer :: ii, jj, kk, ll ! natural orbital labels
      integer :: ia, ib, ic, id ! atomic orbital labels
c
      integer :: iclo,ichi,jclo,jchi
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      nmo = wfn1_param%nmo
      nea = wfn1_param%nea
      neb = wfn1_param%neb
      allocate(occa(1:nmo),occb(1:nmo))
      occa = 0.0d0
      occb = 0.0d0
c
c     Calculate occupation numbers
c     ----------------------------
c
c     Calculate alpha-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,occa)
c
c     Calculate beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,occb)
c
c     Calculate alpha-alpha-occupation numbers
c
      call wfn1_calc_occ2(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                    g_2aa)
c
c     Calculate beta-beta-occupation numbers
c
      call wfn1_calc_occ2(wfn1_wave%g_cb,nmo,nea,wfn1_wave%icnfb,
     &                    g_2bb)
c
c     Calculate entropies
c     -------------------
c
      end subroutine wfn1_e_g_2el
c
c-----------------------------------------------------------------------
