      subroutine ao_fock_2e( geom, basis, nfock, jfac, kfac,
     $     tol2e, oskel, vg_dens, vg_fock, asym)
c
c$Id$
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "schwarz.fh"
#include "cscfps.fh"
#include "cfock.fh"
#include "util.fh"
c
c  Distributed-data AO 2e-Fock construction routine
c
c      F(i) = jfac(i)*J[D(i)] + kfac(i)*K[D(i)]
c
c     arguments
c     
      integer geom, basis            ! [input] parameter handles
      integer nfock                  ! [input] number of Fock matrices
      double precision jfac(nfock)   ! [input] Coulomb prefactor
      double precision kfac(nfock)   ! [input] exchange prefactor
      double precision tol2e         ! [input] integral selection threshold
      logical oskel                  ! [input] toggle skeleton Fock matrix
      integer vg_dens(nfock)         ! [input] array of handles to densities
      integer vg_fock(nfock)         ! [output] array of handles to Fock matrices
      logical asym                   ! [input] flag to antisymmetrize the fock matrices
c     
c     local variables
c     
      integer natoms, maxd, ablklen, nshb, nbf, maxg, lenscr, nsh
      integer i, mem_rep
      integer l_dij, l_dik, l_dli, l_djk, l_dlj, l_dlk
      integer l_fij, l_fik, l_fli, l_fjk, l_flj, l_flk
      integer k_dij, k_dik, k_dli, k_djk, k_dlj, k_dlk
      integer k_fij, k_fik, k_fli, k_fjk, k_flj, k_flk
      integer l_atmp, k_atmp
      integer l_block, k_block, nblock, maxblock
      integer tdim
      parameter (tdim = 100000)
      logical status
      logical int2e_test_mem, int2e_buf_write, int2e_file_rewind
      external int2e_test_mem, int2e_buf_write, int2e_file_rewind
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('ao_fock_2e: bas_numbf failed', 0, BASIS_ERR)
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('ao_fock_2e: bas_numcont failed', 0, BASIS_ERR)
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('ao_fock_2e: geom_ncent failed', 0, GEOM_ERR)
c
c     Default is replicated data.  This can be overridden by setting
c     'fock:replicated' as false.  It is also overriden here when the
c     replicated data algorithm would exceed memory requirements
c
      call intb_mem_2e4c(maxg, lenscr)
      mem_rep = nbf*nbf*2*nfock + nsh*nsh + maxquartet
     $        + ma_sizeof(MT_INT,4*maxquartet+nbf*(nbf+1),MT_DBL)
     $        + nbf*(nbf+1)/2
     $        + ma_sizeof(MT_INT,2*tdim,MT_DBL) 
     $        + 8*ma_sizeof_overhead(MT_DBL)
      if (owritefile) mem_rep = mem_rep + max(0,memsize)
      if (owritefile .or. (.not. (owritefile .or. oreadfile)))
     $      mem_rep = mem_rep + maxeri*5 + lenscr
      if (oreplicated) then
         oreplicated = (mem_rep .lt. ma_inquire_avail(mt_dbl))
      endif
c
      if (oreplicated) then
         call ao_replicated( geom, basis, nfock, jfac, kfac,
     $        tol2e, oskel, vg_dens, vg_fock, asym )
         return
      endif
c     
c     allocate necessary local temporary arrays on the stack
c     
c     l_scr ... workspace for integral routines
c     l_d** ... ** block of density matrix
c     l_f** ... ** block of fock matrix
c     
c     k_* are the offsets corrsponding to the l_* handles
c
c     If reading/writing to/from cache/file then rewind
c
      if (oreadfile .or. owritefile) then
         if (.not. int2e_file_rewind()) call errquit
     $        ('ao_fock_2e: failed rewinding integral file',0,
     &       UNKNOWN_ERR)
      endif
c
c     Halve exchange factor to conform to internal definition
c
      call dscal(nfock, 0.5d0, kfac, 1)
c
c     Determine appropriate task chunking and max no. of bf in a
c     block of the density/fock matrix
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('ao_fock_2e: geom_ncent failed', 0, GEOM_ERR)
      if (.not. ma_push_get(MT_INT, 2*natoms, 'fock2e:block',
     $     l_block, k_block))call errquit('fock_2e: ma failed',2*natoms,
     &       MA_ERR)
      call fock_2e_block_atoms(basis, oskel, tol2e,
     $     int_mb(k_block), nblock, maxblock)
c
      ablklen = maxblock**2
      maxd = nfock*ablklen
c     
      status = .true.
      status = status .and. ma_push_get(MT_DBL, maxd, 'dij',
     $     l_dij, k_dij)
      status = status .and. ma_push_get(MT_DBL, maxd, 'dik',
     $     l_dik, k_dik)
      status = status .and. ma_push_get(MT_DBL, maxd, 'dli',
     $     l_dli, k_dli)
      status = status .and. ma_push_get(MT_DBL, maxd, 'djk',
     $     l_djk, k_djk)
      status = status .and. ma_push_get(MT_DBL, maxd, 'dlj',
     $     l_dlj, k_dlj)
      status = status .and. ma_push_get(MT_DBL, maxd, 'dlk',
     $     l_dlk, k_dlk)
      status = status .and. ma_push_get(MT_DBL, maxd, 'fij',
     $     l_fij, k_fij)
      status = status .and. ma_push_get(MT_DBL, maxd, 'fik',
     $     l_fik, k_fik)
      status = status .and. ma_push_get(MT_DBL, maxd, 'fli',
     $     l_fli, k_fli)
      status = status .and. ma_push_get(MT_DBL, maxd, 'fjk',
     $     l_fjk, k_fjk)
      status = status .and. ma_push_get(MT_DBL, maxd, 'flj',
     $     l_flj, k_flj)
      status = status .and. ma_push_get(MT_DBL, maxd, 'flk',
     $     l_flk, k_flk)
      status = status .and. ma_push_get(MT_DBL, ablklen, 'atmp',
     $     l_atmp, k_atmp)
      if (.not. status) call errquit('fock_2e: d/f** failed', maxd,
     &       MA_ERR)
c
c     If we can hold all of the integrals in core then we should
c     do this ... i.e., the default should be semi-direct, though
c     right now we must store all of the integrals.
c
      if (.not. (oreadfile .or. owritefile))
     $     owritefile = int2e_test_mem(geom, basis, tol2e*0.01d0)
c
      if (oreadfile) then
         call fock_2e_from_file( geom, basis, nfock, ablklen,
     $        jfac, kfac, tol2e, oskel,
     $        dbl_mb(k_dij), dbl_mb(k_dik), dbl_mb(k_dli), 
     $        dbl_mb(k_djk), dbl_mb(k_dlj), dbl_mb(k_dlk), 
     $        dbl_mb(k_fij), dbl_mb(k_fik), dbl_mb(k_fli), 
     $        dbl_mb(k_fjk), dbl_mb(k_flj), dbl_mb(k_flk),
     $        dbl_mb(k_atmp), vg_dens, vg_fock)
      else
c
c        tell the screening statistics we are starting another fock build
c
         iscreen(6) = iscreen(6) + 1
c
         call fock_2e_a( geom, basis, nfock, ablklen,
     $        jfac, kfac, tol2e, oskel,
     $        dbl_mb(k_dij), dbl_mb(k_dik), dbl_mb(k_dli), 
     $        dbl_mb(k_djk), dbl_mb(k_dlj), dbl_mb(k_dlk), 
     $        dbl_mb(k_fij), dbl_mb(k_fik), dbl_mb(k_fli), 
     $        dbl_mb(k_fjk), dbl_mb(k_flj), dbl_mb(k_flk),
     $        dbl_mb(k_atmp), vg_dens, vg_fock,
     $        int_mb(k_block), nblock)
c
c        print the integral screening statistics for this fock build
c        schwarz_print must be called by all nodes.
c
         if (util_print('screening statistics', print_debug)) then
            if ( .not. bas_numcont(basis, nshb) ) call errquit(
     $         'ao_fock_2e: problem with call to bas_numcont', basis,
     &       BASIS_ERR)
            call schwarz_print(natoms, nshb)
         endif
C
      endif
c     
      status = .true.
      status = status .and. ma_pop_stack(l_atmp)
      status = status .and. ma_pop_stack(l_flk)
      status = status .and. ma_pop_stack(l_flj)
      status = status .and. ma_pop_stack(l_fjk)
      status = status .and. ma_pop_stack(l_fli)
      status = status .and. ma_pop_stack(l_fik)
      status = status .and. ma_pop_stack(l_fij)
      status = status .and. ma_pop_stack(l_dlk)
      status = status .and. ma_pop_stack(l_dlj)
      status = status .and. ma_pop_stack(l_djk)
      status = status .and. ma_pop_stack(l_dli)
      status = status .and. ma_pop_stack(l_dik)
      status = status .and. ma_pop_stack(l_dij)
      status = status .and. ma_pop_stack(l_block)
      if (.not. status) call errquit('fock_2e: ma_pop?', 0, MA_ERR)
c     
      call ga_sync()
c     
c     fock_2e_a stuffs contributions into both the lower and upper
c     triangles ... need to symmetrize
c
      do i=1,nfock
        call ga_dscal(vg_fock(i), 4.0d0)
        if (.not. oskel .and. .not. asym) 
     &       call ga_symmetrize(vg_fock(i))
      end do
c
c     Reset exchange factors
c
      call dscal(nfock, 2.d0, kfac, 1)
c
c     Disable writing integrals to cache/file and enable reading for
c     next call if we were just writing.  Also flush trailing block.
c
      if (owritefile) then
         oreadfile = .true.
         owritefile = .false.
         if (.not. int2e_buf_write(.true.))
     $        call errquit('ao_fock_2e: write of last record failed',0,
     &       UNKNOWN_ERR)
      endif
c     
      end
