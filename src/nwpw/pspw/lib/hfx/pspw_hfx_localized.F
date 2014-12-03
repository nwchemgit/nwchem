
*     ************************************************
*     *                                              *
*     *        pspw_potential_hfx_localize           *
*     *                                              *
*     ************************************************

      subroutine pspw_potential_hfx_localize(psi_r,Hpsi_r)
      implicit none
      real*8 psi_r(*)
      real*8 Hpsi_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

      !*** variables to use from pspw_hfx.fh ***
      !* ispin - 1-restricted, 2-unrestricted
      !* norbs(2) - number of up orbitals and down orbitals
      !* n2ft3d   - number of real space grid points ~ (nfft1+2)*nfft2*nfft3/nprocs
      !*******************************************

      !*** the exchange correction has to be added (or appended) to Hpsi_r not just set


      !**** need to set these variables ****
      ehfx = 0.0
      phfx = 0.0  !* phfx = 2*ehfx


      write(*,*) "Hello from pspw_potential_localized_hfx",Hpsi_r(11)


      return
      end
