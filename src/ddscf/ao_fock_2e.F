      subroutine ao_fock_2e( geom, basis, nfock, jfac, kfac,
     $                    tol2e, oskel, vg_dens, vg_fock )
c
c$Id: ao_fock_2e.F,v 1.3 1995-11-02 21:27:39 gg502 Exp $
c
      implicit none
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
c     
c     local variables
c     
      integer natoms, maxd, max_at_bf, ablklen, nshb
      integer iatom, i, ilo, ihi
      integer l_dij, l_dik, l_dli, l_djk, l_dlj, l_dlk
      integer l_fij, l_fik, l_fli, l_fjk, l_flj, l_flk
      integer k_dij, k_dik, k_dli, k_djk, k_dlj, k_dlk
      integer k_fij, k_fik, k_fli, k_fjk, k_flj, k_flk
      integer l_atmp, k_atmp
      logical status
c     
c     allocate necessary local temporary arrays on the stack
c     
c     l_scr ... workspace for integral routines
c     l_d** ... ** block of density matrix
c     l_f** ... ** block of fock matrix
c     
c     k_* are the offsets corrsponding to the l_* handles
c
      if (oscfps) call pstat_on(ps_fock_2e)
c
c     If reading/writing to/from cache/file then rewind
c
      if (oreadfile .or. owritefile) call int2e_file_rewind
c
c     Halve exchange factor to conform to internal definition
c
      call dscal(nfock, 0.5d0, kfac, 1)
c     
c     Determine maximum no. of basis functions on any atom
c     
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('ao_fock_2e: geom_ncent failed', 0)
      max_at_bf = 0
      do iatom = 1, natoms
         if (.not. bas_ce2bfr(basis, iatom, ilo, ihi))
     $        call errquit('ao_fock_2e: bas_ce2bfr failed', iatom)
         max_at_bf = max(max_at_bf, ihi-ilo+1)
      end do
      ablklen = max_at_bf**2
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
      if (.not. status) call errquit('fock_2e: d/f** failed', maxd)
c$$$c     
c$$$c     Zero screening statistics
c$$$c     
c$$$      call ifill(nscreen, 0, iscreen, 1)
c     
      if (oreadfile) then
         call fock_2e_from_file( geom, basis, nfock, ablklen,
     $        jfac, kfac, tol2e, oskel,
     $        dbl_mb(k_dij), dbl_mb(k_dik), dbl_mb(k_dli), 
     $        dbl_mb(k_djk), dbl_mb(k_dlj), dbl_mb(k_dlk), 
     $        dbl_mb(k_fij), dbl_mb(k_fik), dbl_mb(k_fli), 
     $        dbl_mb(k_fjk), dbl_mb(k_flj), dbl_mb(k_flk),
     $        dbl_mb(k_atmp), vg_dens, vg_fock )
      else
C
C        Tell the screening statistics we're starting another Fock build
C
         iscreen(6) = iscreen(6) + 1
C
         call fock_2e_a( geom, basis, nfock, ablklen,
     $        jfac, kfac, tol2e, oskel,
     $        dbl_mb(k_dij), dbl_mb(k_dik), dbl_mb(k_dli), 
     $        dbl_mb(k_djk), dbl_mb(k_dlj), dbl_mb(k_dlk), 
     $        dbl_mb(k_fij), dbl_mb(k_fik), dbl_mb(k_fli), 
     $        dbl_mb(k_fjk), dbl_mb(k_flj), dbl_mb(k_flk),
     $        dbl_mb(k_atmp), vg_dens, vg_fock )
C
C        Print the integral screening statistics for this Fock build
C
         if (ga_nodeid().eq.0 .and. util_print('screening statistics',
     $      print_debug)) then
            if ( .NOT. bas_numcont(basis, nshb) ) Call ErrQuit(
     $         'ao_fock_2e: problem with call to bas_numcont', basis)
            call schwarz_print(natoms, nshb)
         Endif
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
      if (.not. status) call errquit('fock_2e: ma_pop?', 0)
c     
      call ga_sync()
c     
c     fock_2e_a stuffs contributions into both the lower and upper
c     triangles ... need to symmetrize
c
      do i=1,nfock
        call ga_dscal(vg_fock(i), 4.0d0)
        if (.not. oskel) call ga_symmetrize(vg_fock(i))
      end do
c
c     Reset exchange factors
c
      call dscal(nfock, 2.d0, kfac, 1)
c
      if (oscfps) call pstat_off(ps_fock_2e)
c
c     Disable writing integrals to cache/file and enable reading for
c     next call if we were just writing.  Also flush trailing block.
c
      if (owritefile) then
         oreadfile = .true.
         owritefile = .false.
         call int2e_buf_write
c     insert munge of incomplete blocks here.
      endif
c     
      end
