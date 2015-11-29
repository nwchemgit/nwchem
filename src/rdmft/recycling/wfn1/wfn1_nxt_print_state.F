C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Print the next step format wavefunction 
C>
C> Print the components of the wavefunction
C>
      subroutine wfn1_nxt_print_state(state,nbf)
      implicit none
C
#include "wfn1_nxt.fh"
C
      integer nbf !< [Input] The number of basis functions
      double precision state(8*nbf*nbf) !< [Input] The wavefunction
c
      write(*,*)
      write(*,*)'alpha occupation functions'
      write(*,*)'--------------------------'
      call hess_hssout(state(wfn1_nxt_aocc()),nbf,nbf,nbf)
      write(*,*)
      write(*,*)'beta occupation functions'
      write(*,*)'-------------------------'
      call hess_hssout(state(wfn1_nxt_bocc()),nbf,nbf,nbf)
      write(*,*)
      write(*,*)'alpha natural orbitals'
      write(*,*)'----------------------'
      call hess_hssout(state(wfn1_nxt_anat()),nbf,nbf,nbf)
      write(*,*)
      write(*,*)'beta natural orbitals'
      write(*,*)'---------------------'
      call hess_hssout(state(wfn1_nxt_bnat()),nbf,nbf,nbf)
      return
c
      write(*,*)
      write(*,*)'alpha occupation function lagrangian'
      write(*,*)'------------------------------------'
      call hess_hssout(state(wfn1_nxt_aoccl()),nbf,nbf,nbf)
      write(*,*)
      write(*,*)'beta occupation function lagrangian'
      write(*,*)'-----------------------------------'
      call hess_hssout(state(wfn1_nxt_boccl()),nbf,nbf,nbf)
      write(*,*)
      write(*,*)'alpha natural orbital lagrangian'
      write(*,*)'--------------------------------'
      call hess_hssout(state(wfn1_nxt_anatl()),nbf,nbf,nbf)
      write(*,*)
      write(*,*)'beta natural orbital lagrangian'
      write(*,*)'-------------------------------'
      call hess_hssout(state(wfn1_nxt_bnatl()),nbf,nbf,nbf)
C
      end
C> @}
