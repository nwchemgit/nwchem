      subroutine rhf_fock_2e(geom, basis, g_dens, g_fock, tol2e,
     &       ocoul, oexch)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "schwarz.fh"
#include "cscfps.fh"
c     
c     arguments
c     
      integer geom, basis, g_dens, g_fock ! [input] handles
      double precision tol2e    ! [input] integral selection threshold
      logical ocoul, oexch      ! [input] compute coulomb/exchange 
c                               !         contribution to fock matrix
c     
c     local variables
c     
      integer natoms, maxd, iatom, max_at_bf, ilo, ihi
      integer l_dij, l_dik, l_dli, l_djk, l_dlj, l_dlk
      integer l_fij, l_fik, l_fli, l_fjk, l_flj, l_flk
      integer k_dij, k_dik, k_dli, k_djk, k_dlj, k_dlk
      integer k_fij, k_fik, k_fli, k_fjk, k_flj, k_flk
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
c     Determine maximum no. of basis functions on any atom
c     
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('rhf_fock_2e: geom_ncent failed', 0)
      max_at_bf = 0
      do iatom = 1, natoms
         if (.not. bas_ce2bfr(basis, iatom, ilo, ihi))
     $        call errquit('rhf_fock_2e: bas_ce2bfr failed', iatom)
         max_at_bf = max(max_at_bf, ihi-ilo+1)
      enddo
      maxd = max_at_bf**2
c     
      status = ma_push_get(MT_DBL, maxd, 'dij', l_dij, k_dij)
      status = ma_push_get(MT_DBL, maxd, 'dik', l_dik, k_dik)
      status = ma_push_get(MT_DBL, maxd, 'dli', l_dli, k_dli)
      status = ma_push_get(MT_DBL, maxd, 'djk', l_djk, k_djk)
      status = ma_push_get(MT_DBL, maxd, 'dlj', l_dlj, k_dlj)
      status = ma_push_get(MT_DBL, maxd, 'dlk', l_dlk, k_dlk)
      status = ma_push_get(MT_DBL, maxd, 'fij', l_fij, k_fij)
      status = ma_push_get(MT_DBL, maxd, 'fik', l_fik, k_fik)
      status = ma_push_get(MT_DBL, maxd, 'fli', l_fli, k_fli)
      status = ma_push_get(MT_DBL, maxd, 'fjk', l_fjk, k_fjk)
      status = ma_push_get(MT_DBL, maxd, 'flj', l_flj, k_flj)
      status = ma_push_get(MT_DBL, maxd, 'flk', l_flk, k_flk)
      if (.not. status) call errquit('fock_2e: d/f** failed', maxd)
c     
c     Zero screening statistics
c     
      call ifill(nscreen, 0, iscreen, 1)
c     
      call rhf_fock_2e_a(geom, basis, g_dens, g_fock, 
     $     dbl_mb(k_dij), dbl_mb(k_dik), dbl_mb(k_dli), 
     $     dbl_mb(k_djk), dbl_mb(k_dlj), dbl_mb(k_dlk), 
     $     dbl_mb(k_fij), dbl_mb(k_fik), dbl_mb(k_fli), 
     $     dbl_mb(k_fjk), dbl_mb(k_flj), dbl_mb(k_flk), maxd,
     $     tol2e, ocoul, oexch)
c     call ga_print(g_fock)
c     
      status = ma_pop_stack(l_flk)
      status = ma_pop_stack(l_flj)
      status = ma_pop_stack(l_fjk)
      status = ma_pop_stack(l_fli)
      status = ma_pop_stack(l_fik)
      status = ma_pop_stack(l_fij)
      status = ma_pop_stack(l_dlk)
      status = ma_pop_stack(l_dlj)
      status = ma_pop_stack(l_djk)
      status = ma_pop_stack(l_dli)
      status = ma_pop_stack(l_dik)
      status = ma_pop_stack(l_dij)
c     
      call ga_sync()
c     
c     fock_2e_a stuffs contributions into both the lower and upper
c     triangles ... need to symmetrize and multiply by 2
c     
      call ga_symmetrize(g_fock)
      call ga_dscal(g_fock, 2.0d0)
c
      if (oscfps) call pstat_off(ps_fock_2e)
c     
      end
