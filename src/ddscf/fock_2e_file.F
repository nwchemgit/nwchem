      subroutine fock_2e_from_file(geom, basis, nfock, ablklen,
     $     jfac, kfac, tol2e, oskel,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, vg_dens, vg_fock)
      implicit none
#include "util.fh"
#include "cscfps.fh"
c
c     $Id: fock_2e_file.F,v 1.1 1995-12-13 01:24:41 d3g681 Exp $
c
c     Accumulate the contribution to the fock matrices from
c     integrals store in the integral file.  Simply read thru
c     the file getting a range of indices, fetch the corresponding
c     density matrix blocks and then read the integrals in that
c     block.
c
      integer geom, basis
      integer nfock
      integer ablklen
      double precision jfac(nfock)
      double precision kfac(nfock)
      double precision tol2e
      logical oskel
      double precision dij(nfock*ablklen),dik(nfock*ablklen)
      double precision dli(nfock*ablklen),djk(nfock*ablklen)
      double precision dlj(nfock*ablklen),dlk(nfock*ablklen)
      double precision fij(nfock*ablklen),fik(nfock*ablklen)
      double precision fli(nfock*ablklen),fjk(nfock*ablklen)
      double precision flj(nfock*ablklen),flk(nfock*ablklen)
      double precision tmp(ablklen)
      integer vg_dens(nfock)
      integer vg_fock(nfock)
c
      integer maxints
      parameter (maxints = 2048)
      integer i(maxints), j(maxints), k(maxints), l(maxints)
      double precision g(maxints)
      double precision den_tol, denmax, dtol2e
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi
      integer ijk_prev(3,2)
      integer nints, blklen
      logical more
c
      logical int2e_get_bf_range, int2e_file_read
      external int2e_get_bf_range, int2e_file_read
c
      den_tol = max(tol2e*0.01d0,1.0d-300) ! To avoid a hard zero
c
      ijk_prev(1,1) = -1
      ijk_prev(2,1) = -1
      ijk_prev(3,1) = -1
      ijk_prev(1,2) = -1
      ijk_prev(2,2) = -1
      ijk_prev(3,2) = -1
c
      blklen = nfock*ablklen
      call dfill(blklen, 0.0d0, fij, 1)
      call dfill(blklen, 0.0d0, fik, 1)
      call dfill(blklen, 0.0d0, fli, 1)
      call dfill(blklen, 0.0d0, fjk, 1)
      call dfill(blklen, 0.0d0, flj, 1)
      call dfill(blklen, 0.0d0, flk, 1)
c
c     Loop over blocks of integral labels
c
 10   if (int2e_get_bf_range(ilo,ihi,jlo,jhi,klo,khi,llo,lhi)) then
c
c     Get matrices for this block of labels
c
         call fock_2e_cache_dens_fock(
     $        ilo, jlo, klo, llo,
     $        ihi, jhi, khi, lhi,
     $        ijk_prev,
     $        nfock, vg_dens, vg_fock,
     $        jfac, kfac, 
     $        dij, dik, dli, djk, dlj, dlk,
     $        fij, fik, fli, fjk, flj, flk,
     $        tmp)
c
         call fock_density_screen(nfock,
     $        ilo, jlo, klo, llo,
     $        ihi, jhi, khi, lhi,
     $        ilo, jlo, klo, llo,
     $        ihi, jhi, khi, lhi,
     $        dij, dik, dli, djk, dlj, dlk, denmax)
c
         dtol2e = min(1d-7, den_tol/max(1d-10,denmax))
c     
c     Loop over batches of integrals in this block
c     
 20      if (oscfps) call pstat_on(ps_intio)
         more = int2e_file_read(maxints, nints, i, j, k, l, g)
         if (oscfps) call pstat_off(ps_intio)
c
         if (more) then
c
            if (util_print('int2efile read', print_never))
     $           call print_integ_list(nints, i, j, k, l, g)
c
            if (oscfps) call pstat_on(ps_fock_add)
            call fock_2e_label(nfock, dtol2e,
     $           nints, i, j, k, l, g,
     $           ilo, ihi, jlo, jhi, klo, khi, llo, lhi,
     $           dij, dik, dli, djk, dlj, dlk,
     $           fij, fik, fli, fjk, flj, flk)
            if (oscfps) call pstat_off(ps_fock_add)
c
            goto 20
c
         end if
c     
c     Update F blocks
c     
         call fock_upd_blk(nfock, vg_fock,
     $        llo, lhi, ilo, ihi, kfac, fli, tmp)
         call fock_upd_blk(nfock, vg_fock,
     $        llo, lhi, jlo, jhi, kfac, flj, tmp)
         call fock_upd_blk(nfock, vg_fock,
     $        llo, lhi, klo, khi, jfac, flk, tmp)
c     
         goto 10
c     
      end if
c
      if (ijk_prev(1,1) .ne. -1) then
         call fock_upd_blk(nfock, vg_fock, 
     $        ijk_prev(1,1), ijk_prev(1,2),
     $        ijk_prev(2,1), ijk_prev(2,2),
     $        jfac, fij, tmp)
         call fock_upd_blk(nfock, vg_fock,
     $        ijk_prev(2,1), ijk_prev(2,2),
     $        ijk_prev(3,1), ijk_prev(3,2),
     $        kfac, fjk, tmp )
         call fock_upd_blk( nfock, vg_fock,
     $        ijk_prev(1,1), ijk_prev(1,2),
     $        ijk_prev(3,1), ijk_prev(3,2),
     $        kfac, fik, tmp )
      end if
c     
      end
