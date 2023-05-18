      subroutine int_2e_sifs( geom, basis, tol2e, oskel, aoint2, ninfo, 
     &    info, ifmt, ibvtyp, ibitv, clabs)
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
c     driver routine to write 2e integrals to sifs file based on
c     ao_fock_2e
c
c     arguments
c     
      integer geom, basis            ! [input] parameter handles
      integer aoint2 
      integer ninfo, ifmt, ibvtyp
      integer ibitv
      integer info(ninfo)
      double precision tol2e         ! [input] integral selection threshold
      logical oskel                  ! [input] toggle skeleton Fock matrix
c     
c     local variables
c     
      integer natoms, ablklen, nshb, nbf, maxg, lenscr, nsh
      integer i, mem_rep
      integer l_dij, l_dik, l_dli, l_djk, l_dlj, l_dlk
      integer l_fij, l_fik, l_fli, l_fjk, l_flj, l_flk
      integer k_dij, k_dik, k_dli, k_djk, k_dlj, k_dlk
      integer k_fij, k_fik, k_fli, k_fjk, k_flj, k_flk
      integer l_atmp, k_atmp
      integer l_block, k_block, nblock, maxblock
      logical status
      logical int2e_test_mem, int2e_buf_write, int2e_file_rewind
      external int2e_test_mem, int2e_buf_write, int2e_file_rewind
      integer ncomputed, nwritten, num
      common/sif2e/ncomputed, nwritten, num
      integer max2e, l2rec, n2max, ierr
      integer clabs(4,*)
c
cgk debug
*     write(*,*)'gk: entered int_2e_sifs'
cgk end
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('ao_replicated: bas_numbf failed', 0, BASIS_ERR)
cgk debug
*     write(*,*)'gk: nbf= ', nbf
cgk end
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('ao_replicated: bas_numcont failed', 0, BASIS_ERR)
cgk debug
*     write(*,*)'gk: nsh= ', nsh
cgk end
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('ao_fock_2e: geom_ncent failed', 0, GEOM_ERR)
cgk debug
*     write(*,*)'gk: natoms= ', natoms
cgk end
c
c     Default is replicated data.  This can be overridden by setting
c     'fock:replicated' as false.  It is also overriden here when the
c     replicated data algorithm would exceed memory requirements
c
*     call intb_mem_2e4c(maxg, lenscr)
*     mem_rep = nbf*nbf*2*nfock + nsh*nsh
*     if (owritefile) mem_rep = mem_rep + max(0,memsize)
*     if (owritefile .or. (.not. (owritefile .or. oreadfile)))
*    $      mem_rep = mem_rep + maxeri*5 + lenscr
*     if (oreplicated) then
*        oreplicated = (mem_rep .lt. ma_inquire_avail(mt_dbl))
*     endif
c
*     if (oreplicated) then
*        call ao_replicated( geom, basis, nfock, jfac, kfac,
*    $        tol2e, oskel, vg_dens, vg_fock, asym )
*        return
*     endif
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
*     if (oreadfile .or. owritefile) then
*        if (.not. int2e_file_rewind()) call errquit
*    $        ('ao_fock_2e: failed rewinding integral file',0,
*    &       UNKNOWN_ERR)
*     endif
c
c     Halve exchange factor to conform to internal definition
c
*     call dscal(nfock, 0.5d0, kfac, 1)
c
c     Determine appropriate task chunking and max no. of bf in a
c     block of the density/fock matrix
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('ao_fock_2e: geom_ncent failed', 0, GEOM_ERR)
      if (.not. ma_push_get(MT_INT, 2*natoms, 'fock2e:block',
     $     l_block, k_block))call errquit('fock_2e: ma failed',2*natoms,
     &       MA_ERR)
cgk debug
*     write(*,*)'gk: calling fock_2e_block_atoms'
cgk end
      call fock_2e_block_atoms(basis, oskel, tol2e,
     $     int_mb(k_block), nblock, maxblock)
cgk debug
*     write(*,*)'gk: back from fock_2e_block_atoms'
cgk end
c
      ablklen = maxblock**2
*     maxd = nfock*ablklen
c     
      status = .true.
*     status = status .and. ma_push_get(MT_DBL, maxd, 'dij',
*    $     l_dij, k_dij)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'dik',
*    $     l_dik, k_dik)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'dli',
*    $     l_dli, k_dli)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'djk',
*    $     l_djk, k_djk)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'dlj',
*    $     l_dlj, k_dlj)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'dlk',
*    $     l_dlk, k_dlk)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'fij',
*    $     l_fij, k_fij)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'fik',
*    $     l_fik, k_fik)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'fli',
*    $     l_fli, k_fli)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'fjk',
*    $     l_fjk, k_fjk)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'flj',
*    $     l_flj, k_flj)
*     status = status .and. ma_push_get(MT_DBL, maxd, 'flk',
*    $     l_flk, k_flk)
      status = ma_push_get(MT_DBL, ablklen, 'atmp',
     $     l_atmp, k_atmp)
      if (.not. status) call errquit('fock_2e: ablklen failed', ablklen,
     &       MA_ERR)
c
c     If we can hold all of the integrals in core then we should
c     do this ... i.e., the default should be semi-direct, though
c     right now we must store all of the integrals.
c
cgk debug
*     write(*,*)'gk: before int2e_test_mem  oreadfile = ', oreadfile
*     write(*,*)'gk: before int2e_test_mem  owritefile = ', owritefile
cgk end
*     if (.not. (oreadfile .or. owritefile))
*    $     owritefile = int2e_test_mem(geom, basis, tol2e*0.01d0)
cgk debug
*      write(*,*)'gk: after int2e_test_mem  owritefile = ', owritefile
cgk end
c
*     if (oreadfile) then
*        call fock_2e_from_file( geom, basis, nfock, ablklen,
*    $        jfac, kfac, tol2e, oskel,
*    $        dbl_mb(k_dij), dbl_mb(k_dik), dbl_mb(k_dli), 
*    $        dbl_mb(k_djk), dbl_mb(k_dlj), dbl_mb(k_dlk), 
*    $        dbl_mb(k_fij), dbl_mb(k_fik), dbl_mb(k_fli), 
*    $        dbl_mb(k_fjk), dbl_mb(k_flj), dbl_mb(k_flk),
*    $        dbl_mb(k_atmp), vg_dens, vg_fock)
*     else
c
c        tell the screening statistics we're starting another fock build
c
*        iscreen(6) = iscreen(6) + 1
c
cgk debug
*     write(*,*)'gk: calling int_2e_sifs_a'
cgk end

         call int_2e_sifs_a( geom, basis, ablklen, tol2e, oskel,
     $        dbl_mb(k_atmp), int_mb(k_block), nblock,  aoint2, ninfo, 
     &        info, ifmt, ibvtyp, ibitv, clabs)

cgk debug
*     write(*,*)'gk: back from int_2e_sifs_a'
cgk end
c
c        print the integral screening statistics for this fock build
c        schwarz_print must be called by all nodes.
c
         if (util_print('screening statistics', print_debug)) then
            if ( .not. bas_numcont(basis, nshb) ) call errquit(
     $         'ao_fock_2e: problem with call to bas_numcont', basis,
     &       BASIS_ERR)
cgk debug
*     write(*,*)'gk: calling schwarz_print'
cgk end
            call schwarz_print(natoms, nshb)
cgk debug
*     write(*,*)'gk: back from schwarz_print'
cgk end
         endif
C
*     endif
c     
      status = .true.
      status = status .and. ma_pop_stack(l_atmp)
      status = status .and. ma_pop_stack(l_block)
      if (.not. status) call errquit('fock_2e: ma_pop?', 0, MA_ERR)
c     
      call ga_sync()
      write(*,*) 'number of 2e integrals computed = ', ncomputed
      write(*,*) 'number of 2e integrals written to aoints = ', nwritten
      
      end
