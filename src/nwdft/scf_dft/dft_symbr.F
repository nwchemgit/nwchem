      subroutine dft_symbr(geom,basis,ipol,
     .     g_a_dens,g_b_dens, g_tmpa)
c     
c     $Id$
c     
      implicit none
#include "global.fh"
#include "stdio.fh"
      integer geom,basis
      integer  g_a_dens,g_b_dens, g_tmpa
      integer ipol
c
c     This routine is stolen from uhf.F
c     Since UHF can break spatial symmetry by localizing the orbitals
c     the densities may not be totally symmetric, but since the Hamiltonian
c     is symmetric contraction with the integrals projects out the totally 
c     symmetric component ... hence we can symmetrize the densities and
c     exploit symmetry.  Compute the max change in any element due to
c     symmetrizing and print a warning if it is big.
c
      double precision errmaxa,errmaxb
c
      call ga_copy(g_a_dens,g_tmpa)
      call sym_symmetrize(geom, basis, .true., g_a_dens)
      call ga_dadd(1d0, g_a_dens, -1d0, g_tmpa, g_tmpa)
      call ga_maxelt(g_tmpa, errmaxa)
      if (errmaxa.gt.1d-4) then
         if (ga_nodeid().eq.0) then
            write(luout,77) 'alpha',errmaxa
            call util_flush(luout)
         endif
      endif
      if(ipol.eq.1) return
      call ga_copy(g_b_dens,g_tmpa)
      call sym_symmetrize(geom, basis, .true., g_b_dens)
      call ga_dadd(1d0, g_b_dens, -1d0, g_tmpa, g_tmpa)
      call ga_maxelt(g_tmpa, errmaxb)
      if (errmaxb.gt.1d-4) then
         if (ga_nodeid().eq.0) then
            write(luout,77) 'beta',errmaxb
            call util_flush(luout)
         endif
      endif
 77         format(' Warning: spatial symmetry breaking: ',
     $           1p,A5,' = ',1d9.2)
      return
      end

