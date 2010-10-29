      block data cfock_data
c
c$Id$
c
#include "cfock.fh"
      data task_bf/-1/          ! Chunking of funstions per task
      data filesize /0/         ! Max size of integral file on disk
      data memsize  /0/         ! Max size of integral file in memory
      data oreadfile /.false./  ! If reading integrals from file
      data owritefile/.false./  ! If writing integrals to file
      data dentolmax/1d-6/
      data odensityscreen/.true./

c     see scf_init for input of task_chunks etc.

      end
      subroutine fock_2e_a( geom, basis, nfock, ablklen,
     $     jfac, kfac, tol2e, oskel,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, vg_dens, vg_fock,
     $     blocks, nblock)
      implicit none
#include "errquit.fh"
c     
c     $Id$
c     
#include "tcgmsg.fh"
#include "global.fh"
#include "schwarz.fh"
#include "geom.fh"
#include "cfock.fh"
#include "sym.fh"
#include "util.fh"
#include "cscfps.fh"
#include "nwc_const.fh"
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
      double precision ncomputed
      common/ctemp/ncomputed
c
      integer natoms, nproc, ijkl, next
      integer ib, jb, kb, lb, lbhi, i
      integer iatlo, jatlo, katlo, latlo
      integer iathi, jathi, kathi, lathi
      integer ijk_prev(3,2)     ! (i/j/k, lo/hi)
      integer blklen
      double precision sij, smax, skl ! schwarz info
      double precision denmax, dentol, temp, integ_acc
      logical otest
      double precision schwarz_atom_block
      external schwarz_atom_block
      logical sym_atom_block_quartet
      external sym_atom_block_quartet
c     
      integer nxtask
      external nxtask
      external cfock_data ! For T3D linker
c
c     Figure out the threshold on which to screen just integrals
c     by dividing by maximum density matrix element.  10x for safety
c     since really wanted max over atomic blocks
c
      denmax = 1d-10            ! To avoid a hard zero
      do i = 1, nfock
         call ga_maxelt(vg_dens(i), temp)
         denmax = max(temp, denmax)
      enddo
      dentol = min(dentolmax,tol2e/(10.0d0*denmax)) ! For screen integs only
c
c     SEE AO_REPLICATED.F to make the next parameters consistent
c     
c     NOTE THAT ONLY USER SPECIFIED INTEGRAL ACCURACY IS NOW USED
c
      if (intacc .eq. 0.0d0) then
         integ_acc = min(1d-10,max(0.01d0*dentol,1d-20)) ! Variable
*         call int_acc_set(integ_acc)
      else
         integ_acc = intacc     ! User controlled
         call int_acc_set(integ_acc)
      endif
c
c     How many atoms are there?
c     
      ncomputed = 0
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('fock_2e: geom_ncent failed', 0, GEOM_ERR)
      if (natoms .gt. nw_max_atom) call errquit
     $     ('fock_2e: hard dimension failed', natoms, UNKNOWN_ERR)
c     
c     How many processors are there?
c     
      nproc = ga_nnodes()
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
*
*     Cannot use symmetry of pairs since the blocking of atoms means
*     that we may have (kl) > (ij).
*
            if (sij*smax .ge. dentol) then
               do kb = ib, 1, -1
                  katlo = blocks(1,kb)
                  kathi = blocks(2,kb)
                  lbhi = kb
                  if (ib .eq. kb) lbhi = jb
                  do lb = 1, lbhi
                     latlo = blocks(1,lb)
                     lathi = blocks(2,lb)
                     skl = schwarz_atom_block(katlo,kathi,latlo,lathi)
                     otest = (sij*skl) .gt. dentol
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
     $                       tmp, tol2e, dentol)
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
c     Given a quartet of ranges of atoms, determin if there are any
c     interactions that are symmetry unique.
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
c     Given a pair of atom ranges return the max interaction
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
