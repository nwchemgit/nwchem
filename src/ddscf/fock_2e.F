      subroutine fock_2e( geom1, ao_basis, nfock, jfac, kfac,
     $                    tol2e1, oskel1, vg_dens, vg_fock )
c$Id: fock_2e.F,v 1.7 1995-05-24 23:51:24 vg038 Exp $
c
c     wrapper routine for ao_fock_2e and riscf_fock_2e
c
c
      implicit none
#include "cscf.fh"
#include "cfock.fh"

      integer geom1, ao_basis        ! [input] parameter handles
      integer nfock                  ! [input] number of Fock matrices
      double precision jfac(nfock)   ! [input] Coulomb prefactor
      double precision kfac(nfock)   ! [input] exchange prefactor
      double precision tol2e1        ! [input] integral selection threshold
      logical oskel1                 ! [input] toggle skeleton Fock matrix
      integer vg_dens(nfock)         ! [input] array of handles to densities
      integer vg_fock(nfock)         ! [input] array of handles to Fock matrices

      if (nriscf.gt.1) then
        call riscf_fock_2e( geom1, ao_basis, riscf_basis, nfock, 
     $       jfac, kfac, tol2e1, oskel1, vg_dens, vg_fock )
      else if (nriscf.eq.1) then
        if (rifock) then
          call riscf_fock_2e( geom1, ao_basis, riscf_basis, nfock, 
     $         jfac, kfac, tol2e1, oskel1, vg_dens, vg_fock )
        else
          call ao_fock_2e( geom1, ao_basis, nfock, jfac, kfac,
     $         tol2e1, oskel1, vg_dens, vg_fock )
        endif
      else
          call ao_fock_2e( geom1, ao_basis, nfock, jfac, kfac,
     $                    tol2e1, oskel1, vg_dens, vg_fock )
      endif

      return
      end
