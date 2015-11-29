C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Compute the energy and gradient of the current state vector
C>
C> Computes the energy and gradient of the state vector given the
C> various sets of integrals. This routine can generate only the 
C> Hartree-Fock energy expression (which has less degenerate states).
C>
      subroutine wfn1_nxt2_hf_fock(rtdb,geom,nbf,nea,neb,h1,eri,
     &           erix,ov,ov12,ovm12,state,fock,grad,Etot,levelshift)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "wfn1_nxt.fh"
#include "geom.fh"
c
      integer rtdb !< [Input] The RTDB handle (needed for the DFT terms)
      integer geom !< [Input] The geometry handle (needed for the
                   !< nuclear repulsion energy)
      integer nbf  !< [Input] The number of basis functions
      integer nea  !< [Input] The number of alpha electrons
      integer neb  !< [Input] The number of beta electrons
c
      double precision h1(nbf,nbf)           !< [Input] The 1-electron
                                             !< integrals
      double precision eri(nbf,nbf,nbf,nbf)  !< [Input] The Coulomb
                                             !< integrals
      double precision erix(nbf,nbf,nbf,nbf) !< [Input] The exchange
                                             !< integrals
      double precision ov(nbf,nbf)           !< [Input] The overlap
                                             !< integrals
      double precision ov12(nbf,nbf)         !< [Input] The sqrt(overlap)
      double precision ovm12(nbf,nbf)        !< [Input] The 1/sqrt(overlap)
c
      double precision state(8*nbf*nbf+2*nbf) !< [Input] The state vector
      double precision fock(8*nbf*nbf+2*nbf)  !< [Output] The Fock matrices
      double precision grad(8*nbf*nbf+2*nbf)  !< [Output] The gradient
c
      double precision Etot !< [Output] The total energy
      double precision levelshift !< [In/Output] The levelshifter
c
      integer l_occa, k_occa !< Memory for alpha electron occupation
                             !< numbers
      integer l_occb, k_occb !< Memory for beta electron occupation
                             !< numbers
c
      integer l_dfa, k_dfa   !< Memory for alpha electron occupation
                             !< numbers
      integer l_dfb, k_dfb   !< Memory for beta electron occupation
                             !< numbers
c
      integer l_da, k_da !< Memory for the alpha electron density
                         !< matrix
      integer l_db, k_db !< Memory for the beta electron density
                         !< matrix
c
      integer l_fa, k_fa !< Memory for the alpha electron Fock
                         !< matrix
      integer l_fb, k_fb !< Memory for the beta electron Fock
                         !< matrix
c
      integer l_ta, k_ta !< Scratch matrix
      integer l_tb, k_tb !< Scratch matrix
      integer l_t1, k_t1 !< Scratch matrix
      integer l_t2, k_t2 !< Scratch matrix
      integer l_t4, k_t4 !< Scratch matrix
c
      integer k_tv !< Scratch vector
c
      double precision Enuc   !< The nuclear repulsion energy
      double precision E1el   !< The 1-electron energy
      double precision E2elC  !< The 2-electron Coulomb energy
      double precision E2elX  !< The 2-electron eXchange energy
      double precision Ef1el  !< The 1-electron energy
      double precision Ef2elC !< The 2-electron Coulomb energy
      double precision Ef2elX !< The 2-electron eXchange energy
      double precision ELo    !< The occupation function Lagrangian
      double precision ELn    !< The natural orbital Lagrangian
      double precision Ef     !< The Fock energy
c
      double precision t_power !< temporary copy of the exchange power
c
      t_power = wfn1_nxt_x_pow
      wfn1_nxt_x_pow = 1.0d0
      call wfn1_nxt2_energy_gradient(rtdb,geom,nbf,nea,neb,h1,eri,
     &          erix,ov,ov12,ovm12,state,fock,grad,Etot,levelshift)
      wfn1_nxt_x_pow = t_power
c
      end
C>
C> @}
