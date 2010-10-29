      subroutine rhf_fock_2e( geom, basis, g_dens, g_fock, tol2e,
     $     ocoul, oexch, oskel )
*
* $Id$
*
      implicit none
c     
c     Convenience wrapper routine that computes the standard
c     closed shell fock matrix from the input density matrix
c     no questions asked.  It is only called from the initial
c     guess routine with a superposition of atomic densities.
c     Now also called from moints_full to generate an effective
c     frozen-core hamiltonian
c     
      integer geom, basis
      integer g_dens, g_fock    ! [input] handles
      double precision tol2e    ! [input] integral selection threshold
      logical ocoul, oexch      ! [input] compute coulomb/exchange 
      logical oskel             ! [input] compute skeleton fock matrix
c     
c     
      integer vg_dens(2)
      integer vg_fock(2)
      integer nfock
      double precision jfac(2), kfac(2)

      if (ocoul) then
         jfac(1) = 1.d0
      else
         jfac(1) = 0.d0
      endif
      if (oexch) then
         kfac(1) = -0.5d0
      else
         kfac(1) = 0.d0
      endif
      nfock = 1
      vg_dens(1) = g_dens
      vg_fock(1) = g_fock

      call fock_2e( geom, basis, nfock, jfac, kfac, tol2e, oskel,
     $     vg_dens, vg_fock, .false.)
      
      return
      end
