      block data cfock_data
c
c$Id: fock_2e_a.F,v 1.6 1995-04-05 19:51:33 d3g681 Exp $
c
#include "cfock.fh"
      data task_chunks/-1/

c     see scf_init for input of task_chunks

      end
      subroutine fock_2e_a( geom, basis, nfock, ablklen,
     $     jfac, kfac, tol2e, oskel,
     $     dij, dik, dli, djk, dlj, dlk,
     $     fij, fik, fli, fjk, flj, flk,
     $     tmp, vg_dens, vg_fock )
      implicit none
c     
c     $Id: fock_2e_a.F,v 1.6 1995-04-05 19:51:33 d3g681 Exp $
c     
#include "tcgmsg.fh"
#include "global.fh"
#include "schwarz.fh"
#include "geom.fh"
#include "cfock.fh"
#include "sym.fh"
#include "util.fh"
c     
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
c     
c     
      integer natoms, nproc, ijkl, next
      integer iat, jat, kat, lat, kl, lathi
      integer idim, jdim, kdim, ldim
      integer i_prev, j_prev, k_prev
      integer blklen
      double precision denmax   ! Max Magnitude of atomic density blocks
      double precision sij, smax, skl(1000), sijkl ! schwarz info
      integer lat_list(1000), nkl ! list of lat interacting with kat
      double precision den_tol  ! Density screen tolerance ( tol2e *  factor )
      double precision q4_atom  ! Symmetry scaling factor
      double precision q4_junk  ! Scratch
      double precision dtol2e
      integer ntask
      logical otest, otestij
c     
c     
      integer nxtask
      double precision dabsmax, dabssum
      external dabsmax, dabssum, nxtask
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
c     For load imbalance to be about 1% we need about 100 tasks per
c     processor.  Many more than this is unecessary.
c     Task_chunks controls parallelization of the nested ijk.  
c     The available parallelism is of order natoms**3.
c     The default values of -1 is used to check for user input
c     
      ntask = max(natoms**3 / nproc, 1)
      if (task_chunks .eq. -1) then
         task_chunks = max(ntask/100,1)
      end if
      if (ga_nodeid().eq.0 .and. util_print('ntask',print_high)) then
         write(6,*) ' ntask, task_chunks, ',ntask, task_chunks
         call util_flush(6)
      end if
c     
c     Empirically observe that to compute the fock matrix/energy
c     to an accuracy of t need to screen on density*integrals
c     to a precision of about 10^-2*t
c     
      den_tol = max(tol2e*0.01d0,1.0d-300) ! To avoid a hard zero
c     
      i_prev = -1               ! For caching of F/D blocks
      j_prev = -1
      k_prev = -1
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
         nkl = 0
         do lat = 1, kat
            skl(lat) =  schwarz_atom(kat,lat)
            if (smax*skl(lat) .ge. tol2e) then
               nkl = nkl + 1
               lat_list(nkl) = lat
            end if
         end do
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
c     Load balance now only over non-zero interactions
c     
                     q4_atom = 1.0d0
                     otest = (lat .le. lathi) .and.
     $                    (sij*skl(lat) .ge. tol2e)
                     if (otest .and. oskel) then
                        otest = sym_atom_quartet(geom, 
     $                       iat, jat, kat, lat, q4_atom)
                     end if
c     
                     if (otest .and. (ijkl .eq. next)) then
                        iscreen(3) = iscreen(3) + 1
c     
c     Get blocks of D/F for coulomb interaction
c     
                        if (i_prev.ne.iat .or. j_prev.ne.jat) then
                           if (i_prev.ne.-1)
     $                          call fock_upd_atom_blk( nfock, vg_fock,
     $                          basis, i_prev, j_prev, jfac, fij, tmp )
                           call fock_get_atom_blk( nfock, vg_dens,
     $                          basis, iat, jat, dij, tmp, idim, jdim )
                           call dfill((idim*jdim*nfock), 0.0d0, fij, 1)
                        end if
                        call fock_get_atom_blk( nfock, vg_dens, basis,
     $                       lat, kat, dlk, tmp, ldim, kdim )
                        call dfill((kdim*ldim*nfock), 0.0d0, flk, 1)
c     
c     Get blocks of D/F for exchange interaction
c     
                        if (i_prev.ne.iat .or. k_prev.ne.kat) then
                           if (i_prev.ne.-1)
     $                          call fock_upd_atom_blk( nfock, vg_fock,
     $                          basis, i_prev, k_prev, kfac, fik, tmp )
                           call fock_get_atom_blk( nfock, vg_dens,
     $                          basis, iat, kat, dik, tmp, idim, kdim )
                           call dfill((idim*kdim*nfock), 0.0d0, fik, 1)
                        end if
                        if (j_prev.ne.jat .or. k_prev.ne.kat) then
                           if (j_prev .ne. -1)
     $                          call fock_upd_atom_blk( nfock, vg_fock,
     $                          basis, j_prev, k_prev, kfac, fjk, tmp )
                           call fock_get_atom_blk( nfock, vg_dens,
     $                          basis, jat, kat, djk, tmp, jdim, kdim )
                           call dfill((jdim*kdim*nfock), 0.0d0, fjk, 1)
                        end if
                        call fock_get_atom_blk( nfock, vg_dens, basis,
     $                       lat, iat, dli, tmp, ldim, idim )
                        call fock_get_atom_blk( nfock, vg_dens, basis,
     $                       lat, jat, dlj, tmp, ldim, jdim )
                        call dfill((idim*ldim*nfock), 0.0d0, fli, 1)
                        call dfill((jdim*ldim*nfock), 0.0d0, flj, 1)
c     
                        i_prev = iat
                        j_prev = jat
                        k_prev = kat
c     
c     Apply screening on the density and schwarz inequality
c     
                        call fock_density_screen(nfock,
     $                       idim, jdim, kdim, ldim,
     $                       dij, dik, dli, djk, dlj, dlk, denmax)
c
                        sijkl = schwarz_atom(iat,jat) *
     $                       schwarz_atom(kat,lat)
c
                        if (denmax*sijkl .ge. den_tol) then
c     
c     Construct these Fock atom-blocks
c     
                           dtol2e = den_tol/denmax
                           call fock_2e_b( basis, nfock, dtol2e,
     $                          q4_atom,
     $                          iat, jat, kat, lat, 
     $                          idim, jdim, kdim, ldim,
     $                          dij, dik, dli, djk, dlj, dlk,
     $                          fij, fik, fli, fjk, flj, flk )
c     
c     Update F blocks
c     
                           call fock_upd_atom_blk( nfock, vg_fock,
     $                          basis, lat, iat, kfac, fli, tmp)
                           call fock_upd_atom_blk( nfock, vg_fock,
     $                          basis, lat, jat, kfac, flj, tmp)
                           call fock_upd_atom_blk( nfock, vg_fock,
     $                          basis, lat, kat, jfac, flk, tmp)
                        end if
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
      if (i_prev.ne.-1) then
         call fock_upd_atom_blk( nfock, vg_fock, basis,
     $        i_prev, j_prev, jfac, fij, tmp )
         call fock_upd_atom_blk( nfock, vg_fock, basis,
     $        j_prev, k_prev, kfac, fjk, tmp )
         call fock_upd_atom_blk( nfock, vg_fock, basis,
     $        i_prev, k_prev, kfac, fik, tmp )
      end if
c     
      next = nxtask(-nproc, task_chunks)
      call ga_sync()
c     
      end
