      subroutine fock_2e_cache_dens_fock(
     $     ilo, jlo, klo, llo,
     $     ihi, jhi, khi, lhi,
     $     ijk_prev,
     $     nfock, vg_dens, vg_fock,
     $     jfac, kfac, 
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp)
      implicit none
#include "cscfps.fh"
c
c     $Id$
c
c     For the given ranges of i,j,k,l fetch the six blocks
c     of the density matrices (ij,ik,il,jk,jl,kl) and accumulate
c     the fock matrix blocks.  Scale the fock matrix blocks
c     by the necessary factors for Coulomb/Exchange.
c     Take advantage of caching by storing in ijk_prev() the
c     previous i,j,k ranges (l assumed to be inner loop and
c     thus not to benefit from caching).  If any of ij,ik,jk
c     have not changed then nothing need be done for those
c     blocks.
c
c     !! ijk_prev() must be preserved between calls and set
c     to invalid basis function labels (-1) on the first call.
c     
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi
      integer ijk_prev(3,2)     ! old values (i/j/k, lo/hi)
      integer nfock             ! No. of matrices
      integer vg_dens(nfock), vg_fock(nfock) ! GA handles
      double precision jfac(nfock), kfac(nfock)
      double precision
     $     dij(nfock,*), dik(nfock,*), dli(nfock,*),
     $     djk(nfock,*), dlj(nfock,*), dlk(nfock,*)
      double precision
     $     fij(nfock,*), fik(nfock,*), fli(nfock,*),
     $     fjk(nfock,*), flj(nfock,*), flk(nfock,*)
      double precision tmp(*)
c
      integer idim, jdim, kdim, ldim
c
      if (oscfps) call pstat_on(ps_gaops)
c
c     Assume that blocking is fixed and therefore only
c     need to compare ilo with ilo_prev, etc.
c
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
      kdim = khi - klo + 1
      ldim = lhi - llo + 1
c     
c     Blocks of D/F without l label ... caching is useful
c     
      if (ijk_prev(1,1).ne.ilo .or. ijk_prev(2,1).ne.jlo) then
         if (ijk_prev(1,1) .ne.-1)
     $        call fock_upd_blk( nfock, vg_fock,
     $        ijk_prev(1,1), ijk_prev(1,2),
     $        ijk_prev(2,1), ijk_prev(2,2),
     $        jfac, fij, tmp )
         call fock_get_blk(nfock, vg_dens,
     $        ilo, ihi, jlo, jhi, dij, tmp)
         call dfill((idim*jdim*nfock), 0.0d0, fij, 1)
      end if
      if (ijk_prev(1,1).ne.ilo .or. ijk_prev(3,1).ne.klo) then
         if (ijk_prev(1,1).ne.-1)
     $        call fock_upd_blk(nfock, vg_fock,
     $        ijk_prev(1,1), ijk_prev(1,2),
     $        ijk_prev(3,1), ijk_prev(3,2), kfac, fik, tmp)
         call fock_get_blk(nfock, vg_dens,
     $        ilo, ihi, klo, khi, dik, tmp)
         call dfill((idim*kdim*nfock), 0.0d0, fik, 1)
      end if
      if (ijk_prev(2,1).ne.jlo .or. ijk_prev(3,1).ne.klo) then
         if (ijk_prev(2,1) .ne. -1)
     $        call fock_upd_blk(nfock, vg_fock,
     $        ijk_prev(2,1), ijk_prev(2,2),
     $        ijk_prev(3,1), ijk_prev(3,2), kfac, fjk, tmp)
         call fock_get_blk(nfock, vg_dens,
     $        jlo, jhi, klo, khi, djk, tmp)
         call dfill((jdim*kdim*nfock), 0.0d0, fjk, 1)
      end if
c
c     Blocks with l label always change
c
      call fock_get_blk(nfock, vg_dens,
     $     llo, lhi, klo, khi, dlk, tmp)
      call fock_get_blk(nfock, vg_dens,
     $     llo, lhi, ilo, ihi, dli, tmp)
      call fock_get_blk(nfock, vg_dens,
     $     llo, lhi, jlo, jhi, dlj, tmp)
      call dfill((kdim*ldim*nfock), 0.0d0, flk, 1)
      call dfill((idim*ldim*nfock), 0.0d0, fli, 1)
      call dfill((jdim*ldim*nfock), 0.0d0, flj, 1)
c     
      ijk_prev(1,1) = ilo
      ijk_prev(2,1) = jlo
      ijk_prev(3,1) = klo
      ijk_prev(1,2) = ihi
      ijk_prev(2,2) = jhi
      ijk_prev(3,2) = khi
c     
      if (oscfps) call pstat_off(ps_gaops)
c
      end

