      subroutine fock_2e_from_file(geom, basis, nfock, ablklen,
     $     jfac, kfac, tol2e, oskel,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, vg_dens, vg_fock)
      implicit none
#include "util.fh"
#include "cscfps.fh"
#include "cfock.fh"
c
c     $Id$
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
      double precision den_tol, denmax, dtol2e
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi
      integer ijk_prev(3,2)
      integer blklen
c
      logical int2e_get_bf_range, int2e_file_read
      external int2e_get_bf_range, int2e_file_read
c
      if (oscfps) call pstat_on(ps_fock_io)
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
         call fock_init_cmul(ihi-ilo+1,jhi-jlo+1,lhi-llo+1)
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
         dtol2e = min(dentolmax, den_tol/max(1d-10,denmax), 
     $        den_tol/max(1d-10,denmax**2))
c     
         call int2e_file_fock_block(nfock, dtol2e,
     $           dij, dik, dli, djk, dlj, dlk,
     $           fij, fik, fli, fjk, flj, flk)
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
      if (oscfps) call pstat_off(ps_fock_io)
c
      end
      subroutine fock_2e_rep_from_file(geom, basis, nfock, nbf,
     $     jfac, kfac, tol2e, oskel, dens, fock)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "cscfps.fh"
#include "cfock.fh"
#include "global.fh"
      integer geom, basis
      integer nfock, nbf
      double precision jfac(nfock), kfac(nfock)
      double precision tol2e
      logical oskel
      double precision dens(nfock,nbf*nbf), fock(nfock,nbf*nbf)
c
      double precision den_tol, denmax
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi, i, j
c
      logical int2e_get_bf_range, int2e_file_read
      external int2e_get_bf_range, int2e_file_read
      integer idamax
      external idamax
c
      if (oscfps) call pstat_on(ps_fock_io)
c
      denmax = 0.0d0
      do i = 1, nfock
         j = idamax(nbf*nbf, dens(i,1), nfock)
         denmax = max(denmax, abs(dens(i,j)))
      end do
      den_tol = min(dentolmax,tol2e/denmax,tol2e/denmax**2) ! Threshold to screen integs only
c
      if (ga_nodeid().eq.0 .and. util_print('fockfile',print_debug))
     $     write(6,*) 'fockfile: tols ',tol2e, dentolmax, denmax,den_tol
c
      call fock_init_cmul(nbf,nbf,nbf) ! lookup table for f build
c
c     Loop over blocks of integral labels
c
 10   if (int2e_get_bf_range(ilo,ihi,jlo,jhi,klo,khi,llo,lhi)) then
         if (ilo.ne.1 .or. jlo.ne.1 .or. klo.ne.1 .or. llo.ne.1)
     $        call errquit('not a replicated file?',0, UNKNOWN_ERR)
         if (ihi.ne.nbf .or. jhi.ne.nbf .or. khi.ne.nbf.or. lhi.ne.nbf)
     $        call errquit('not a replicated file?',0, UNKNOWN_ERR)
c
         call int2e_file_rep_fock(nfock, nbf,
     $        jfac, kfac, den_tol, dens, fock)
c
         goto 10
c     
      end if
c
      if (oscfps) call pstat_off(ps_fock_io)
c
      end


