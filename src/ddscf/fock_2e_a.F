      block data cfock_data
c
c$Id: fock_2e_a.F,v 1.20 1996-01-10 16:44:51 d3g681 Exp $
c
#include "cfock.fh"
      data task_bf/-1/          ! Chunking of funstions per task
      data filesize /0/         ! Max size of integral file on disk
      data memsize  /0/         ! Max size of integral file in memory
      data oreadfile /.false./  ! If reading integrals from file
      data owritefile/.false./  ! If writing integrals to file

c     see scf_init for input of task_chunks etc.

      end
      subroutine fock_2e_a( geom, basis, nfock, ablklen,
     $     jfac, kfac, tol2e, oskel,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, vg_dens, vg_fock,
     $     blocks, nblock)
      implicit none
c     
c     $Id: fock_2e_a.F,v 1.20 1996-01-10 16:44:51 d3g681 Exp $
c     
#include "tcgmsg.fh"
#include "global.fh"
#include "schwarz.fh"
#include "geom.fh"
#include "cfock.fh"
#include "sym.fh"
#include "util.fh"
#include "cscfps.fh"
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
      integer blocks(2,*)
      integer nblock
c     
c
      double precision ncomputed
      common/crap/ncomputed
c
      integer natoms, nproc, ijkl, next
      integer ib, jb, kb, lb, lbhi
      integer iatlo, jatlo, katlo, latlo
      integer iathi, jathi, kathi, lathi
      integer ijk_prev(3,2)     ! (i/j/k, lo/hi)
      integer blklen
      double precision sij, smax, skl ! schwarz info
      double precision den_tol  ! Density screen tolerance ( tol2e *  factor )
      logical otest, otestij
      double precision schwarz_atom_block
      external schwarz_atom_block
      logical sym_atom_block_pair, sym_atom_block_quartet
      external sym_atom_block_pair, sym_atom_block_quartet
c     
      integer nxtask
      external nxtask
      external cfock_data ! For T3D linker
c
c     How many atoms are there?
c     
      ncomputed = 0
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('fock_2e: geom_ncent failed', 0)
      if (natoms .gt. 1000) call errquit
     $     ('fock_2e: hard dimension failed', natoms)
c     
c     How many processors are there?
c     
      nproc = ga_nnodes()
c     
c     Empirically observe that to compute the fock matrix/energy
c     to an accuracy of t need to screen on density*integrals
c     to a precision of about 10^-2*t
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
      ijkl = 0
      next = nxtask(nproc, 1)
c     
c     Loop thru blocked atomic quartets
c
      smax = schwarz_max()
      do ib = nblock, 1, -1
         iatlo = blocks(1,ib)
         iathi = blocks(2,ib)
         if (util_print('report loop', print_debug)) 
     $        call fock_2e_report(iathi, natoms)
         do jb = 1, ib
            jatlo = blocks(1,jb)
            jathi = blocks(2,jb)
            sij = schwarz_atom_block(iatlo,iathi,jatlo,jathi)
            otestij = sij*smax .ge. tol2e
*
*     Cannot use symmetry of pairs since the blocking of atoms means
*     that we may have (kl) > (ij).
*
*            if (otestij .and. oskel) otestij = sym_atom_block_pair
*     $           (geom, iatlo, iathi, jatlo, jathi)
            if (otestij) then
               do kb = ib, 1, -1
                  katlo = blocks(1,kb)
                  kathi = blocks(2,kb)
                  lbhi = kb
                  if (ib .eq. kb) lbhi = jb
                  do lb = 1, lbhi
                     latlo = blocks(1,lb)
                     lathi = blocks(2,lb)
                     skl = schwarz_atom_block(katlo,kathi,latlo,lathi)
                     otest = (sij*skl) .gt. tol2e
                     if (otest .and. oskel) otest=sym_atom_block_quartet
     $                    (geom, iatlo, iathi, jatlo, jathi, 
     $                       katlo, kathi, latlo, lathi)
c     
c     Load balance over non-zero interactions
c     
                     if (otest .and. (ijkl .eq. next)) then
                        call fock_2e_task(
     $                       geom, basis, oskel,
     $                       iatlo, jatlo, katlo, latlo,
     $                       iathi, jathi, kathi, lathi,
     $                       ijk_prev,
     $                       nfock, vg_dens, vg_fock,
     $                       jfac, kfac, 
     $                       dij, dik, dli, djk, dlj, dlk,
     $                       fij, fik, fli, fjk, flj, flk,
     $                       tmp, tol2e, den_tol)
c
                        next = nxtask(nproc, 1)
                     end if
                     if (otest) ijkl = ijkl + 1
                  end do
               enddo
            endif
         end do
      end do
c     
      if (ijk_prev(1,1) .ne. -1) then
         if (oscfps) call pstat_on(ps_gaops)
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
         if (oscfps) call pstat_off(ps_gaops)
      endif
c
      if (util_print('ncomputed',print_debug))
     $   write(6,*) ' ncomputed ', ncomputed
c     
      next = nxtask(-nproc, 1)
      call ga_sync()
c     
      end
      logical function sym_atom_block_pair(geom, 
     $     iatlo, iathi, jatlo, jathi)
      implicit none
#include "sym.fh"
      integer geom
      integer iatlo, iathi, jatlo, jathi
c
      double precision q2
      integer iat, jat, jhi
      logical status, oij
c
      status = .false.
      oij = iatlo .eq. jatlo
      do iat = iatlo, iathi
         jhi = jathi
         if (oij) jhi = iat
         do jat = jatlo, jhi
            status = sym_atom_pair(geom, iat, jat, q2)
            if (status) goto 10
         enddo
      enddo
c
 10   sym_atom_block_pair = status
c
      end
      logical function sym_atom_block_quartet(geom, 
     $     iatlo, iathi, jatlo, jathi, katlo, kathi, latlo, lathi)
      implicit none
#include "sym.fh"
      integer geom
      integer iatlo, iathi, jatlo, jathi, katlo, kathi, latlo, lathi
c
      integer iat, jat, kat, lat, lattop, jattop, kattop
      double precision q4
      logical otest, oij, okl, oikjl
c
      oij   = iatlo.eq.jatlo
      okl   = katlo.eq.latlo
      oikjl = (iatlo.eq.katlo) .and. (jatlo.eq.latlo)
c
      otest = .false.
      do iat = iatlo, iathi
         jattop = jathi
         if (oij) jattop = iat
         do jat = jatlo, jattop
            kattop = kathi
            if (oikjl) kattop = iat
            do kat = katlo, kattop
               lattop = lathi
               if (okl) lattop = kat
               if (oikjl .and. kat.eq.iat) lattop = jat
               do lat = latlo,lattop
                  otest = sym_atom_quartet(geom, 
     $                 iat, jat, kat, lat, q4)
                  if (otest) goto 10
               enddo
            enddo
         enddo
      enddo
c
 10   sym_atom_block_quartet = otest
c
      end
      double precision function schwarz_atom_block(
     $     iatlo,iathi,jatlo,jathi)
      implicit none
#include "schwarz.fh"
      integer iatlo, iathi, jatlo, jathi
c
      integer iat, jat, jhi
      double precision s
      logical oij
c
      s = 0.0d0
      oij = iatlo .eq. jatlo
      do iat = iatlo, iathi
         jhi = jathi
         if (oij) jhi = iat
         do jat = jatlo, jhi
            s = max(s, schwarz_atom(iat,jat))
         enddo
      enddo
c
      schwarz_atom_block = s
c
      end
