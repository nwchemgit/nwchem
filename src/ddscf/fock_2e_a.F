      block data cfock_data
c
c$Id: fock_2e_a.F,v 1.14 1995-11-13 06:41:14 d3g681 Exp $
c
#include "cfock.fh"
      data task_chunks/-1/      ! Defines size of parallel tasks
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
     $     tmp, vg_dens, vg_fock )
      implicit none
c     
c     $Id: fock_2e_a.F,v 1.14 1995-11-13 06:41:14 d3g681 Exp $
c     
#include "tcgmsg.fh"
#include "global.fh"
#include "schwarz.fh"
#include "geom.fh"
#include "cfock.fh"
#include "sym.fh"
#include "util.fh"
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
      integer natoms, nproc, ijkl, next
      integer iat, jat, kat, lat, kl, lathi
      integer ijk_prev(3,2)     ! (i/j/k, lo/hi)
      integer blklen
      double precision sij, smax, skl(1000), sijkl ! schwarz info
      integer lat_list(1000), nkl ! list of lat interacting with kat
      double precision den_tol  ! Density screen tolerance ( tol2e *  factor )
      double precision q4_atom  ! Symmetry scaling factor
      double precision q4_junk  ! Scratch
      logical otest, otestij
c     
      integer nxtask
      external nxtask
c
c     How many atoms are there?
c     
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
c     Determine appropriate task chunking
c
      call fock_2e_est_task(geom, oskel, tol2e)
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
      q4_atom = 1.0d0
c     
      ijkl = 0
      next = nxtask(nproc, task_chunks)
c     
c     Loop thru atomic triplets (i,j,k)
c     
      smax = schwarz_max()
      do kat = natoms, 1, -1
         if (util_print('report loop', print_debug)) 
     $        call fock_2e_report(kat, natoms)
c
c     Build list of atoms lat that overlap with atom kat
c
         nkl = 0
         do lat = 1, kat
            skl(lat) =  schwarz_atom(kat,lat)
            if (smax*skl(lat) .ge. tol2e) then
               nkl = nkl + 1
               lat_list(nkl) = lat
            end if
         end do
c
         do iat = natoms, kat, -1
            do jat = iat, 1, -1
               sij = schwarz_atom(iat,jat)
               otestij = sij*smax .ge. tol2e
               if (otestij .and. oskel) then
                  otestij = sym_atom_pair(geom, iat, jat, q4_junk)
               end if
               if (otestij) then
                  lathi = kat
                  if (iat .eq. kat) lathi = jat
c     
                  do kl = 1, nkl
                     lat = lat_list(kl)
c     
c     Load balance over non-zero interactions
c     
                     q4_atom = 1.0d0
                     sijkl = sij*skl(lat)
                     otest = (lat .le. lathi) .and.
     $                    (sijkl .ge. tol2e)
                     if (otest .and. oskel) then
                        otest = sym_atom_quartet(geom, 
     $                       iat, jat, kat, lat, q4_atom)
                     end if
c     
                     if (otest .and. (ijkl .eq. next)) then
                        iscreen(3) = iscreen(3) + 1
c
                        call fock_2e_task(
     $                       geom, basis,
     $                       iat, jat, kat, lat,
     $                       ijk_prev,
     $                       nfock, vg_dens, vg_fock,
     $                       jfac, kfac, 
     $                       dij, dik, dli, djk, dlj, dlk,
     $                       fij, fik, fli, fjk, flj, flk,
     $                       tmp,
     $                       sijkl, q4_atom, tol2e, den_tol)
c
                        next = nxtask(nproc, task_chunks)
                     end if
                     if (otest) ijkl = ijkl + 1
                  end do
               else
                  if (ga_nodeid().eq.0 .and. kat.eq.1) 
     $                 iscreen(1) = iscreen(1) + 1
               end if
            end do
         end do
      end do
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
      next = nxtask(-nproc, task_chunks)
      call ga_sync()
c     
      end
