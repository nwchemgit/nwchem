      subroutine fock_2e_task(
     $     geom, basis,
     $     iat, jat, kat, lat,
     $     ijk_prev,
     $     nfock, vg_dens, vg_fock,
     $     jfac, kfac, 
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp,
     $     sijkl, q4_atom, tol2e, den_tol)
      implicit none
#include "bas.fh"
#include "cfock.fh"
c
c     $Id: fock_2e_task.F,v 1.5 1995-11-03 20:41:55 d3g681 Exp $
c
c     Given an atomic quartet, fetch the necessary blocks
c     of the density matrices, call fock_2e_b to add in
c     in the fock matrix contribution and then accumulate
c     the contributions.
c
c     This routine will change soon to accomodate more general
c     blocking strategies
c
      integer geom
      integer basis
      integer iat, jat, kat, lat
      integer ijk_prev(3,2)
      integer nfock
      integer vg_dens(nfock), vg_fock(nfock)
      double precision jfac(nfock), kfac(nfock)
      double precision dij(nfock,*), dik(nfock,*), dli(nfock,*),
     $     djk(nfock,*), dlj(nfock,*), dlk(nfock,*)
      double precision fij(nfock,*), fik(nfock,*), fli(nfock,*),
     $     fjk(nfock,*), flj(nfock,*), flk(nfock,*)
      double precision tmp(*)
      double precision sijkl
      double precision q4_atom
      double precision tol2e, den_tol
c
c     Given a task which comprises computing the contributions
c     arising from a quartet of atoms, do everything necessary
c
      double precision denmax, dtol2e
      logical otest
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
c
c     With caching get density and flush necessary fock blocks
c     and return magnitude of density blocks for screening
c
      if (.not. bas_ce2bfr(basis, iat, ilo, ihi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
      if (.not. bas_ce2bfr(basis, jat, jlo, jhi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
      if (.not. bas_ce2bfr(basis, kat, klo, khi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
      if (.not. bas_ce2bfr(basis, lat, llo, lhi))
     $     call errquit('fock_2e_task: ce2bfr failed', 0)
      call fock_2e_cache_dens_fock(
     $     ilo, jlo, klo, llo,
     $     ihi, jhi, khi, lhi,
     $     ijk_prev,
     $     nfock, vg_dens, vg_fock,
     $     jfac, kfac, 
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     denmax, tmp)
c
c     Apply screening on the density and schwarz inequality
c     (do not use density screening if writing integrals to file)
c
      if (owritefile) then
         otest  = sijkl .gt. den_tol
      else
         otest  = (sijkl*denmax) .gt. den_tol
      endif
c
      if (otest) then
c
         if (owritefile) then
            dtol2e = den_tol
         else
            dtol2e = den_tol / denmax
         endif
c     
c     Construct these Fock atom-blocks
c     
         call fock_2e_b(basis, nfock, dtol2e,
     $        q4_atom,
     $        iat, jat, kat, lat, 
     $        dij, dik, dli, djk, dlj, dlk,
     $        fij, fik, fli, fjk, flj, flk)
c     
c     Update F blocks
c     
         call fock_upd_blk(nfock, vg_fock,
     $        llo, lhi, ilo, ihi, kfac, fli, tmp)
         call fock_upd_blk(nfock, vg_fock,
     $        llo, lhi, jlo, jhi, kfac, flj, tmp)
         call fock_upd_blk(nfock, vg_fock,
     $        llo, lhi, klo, khi, jfac, flk, tmp)
      end if
c     
      end
