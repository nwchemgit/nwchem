      subroutine rhf_fock_2e(geom, basis, g_dens, g_fock, tol2e,
     &       ocoul, oexch)
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "schwarz.fh"
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
      integer l_dij, l_dik, l_dil, l_djk, l_djl, l_dkl
      integer l_fij, l_fik, l_fil, l_fjk, l_fjl, l_fkl
      integer k_dij, k_dik, k_dil, k_djk, k_djl, k_dkl
      integer k_fij, k_fik, k_fil, k_fjk, k_fjl, k_fkl
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
      status = ma_push_get(MT_DBL, maxd, 'dil', l_dil, k_dil)
      status = ma_push_get(MT_DBL, maxd, 'djk', l_djk, k_djk)
      status = ma_push_get(MT_DBL, maxd, 'djl', l_djl, k_djl)
      status = ma_push_get(MT_DBL, maxd, 'dkl', l_dkl, k_dkl)
      status = ma_push_get(MT_DBL, maxd, 'fij', l_fij, k_fij)
      status = ma_push_get(MT_DBL, maxd, 'fik', l_fik, k_fik)
      status = ma_push_get(MT_DBL, maxd, 'fil', l_fil, k_fil)
      status = ma_push_get(MT_DBL, maxd, 'fjk', l_fjk, k_fjk)
      status = ma_push_get(MT_DBL, maxd, 'fjl', l_fjl, k_fjl)
      status = ma_push_get(MT_DBL, maxd, 'fkl', l_fkl, k_fkl)
      if (.not. status) call errquit('fock_2e: d/f** failed', maxd)
c     
c     Zero screening statistics
c     
      call ifill(nscreen, 0, iscreen, 1)
c     
      call rhf_fock_2e_a(geom, basis, g_dens, g_fock, 
     $     dbl_mb(k_dij), dbl_mb(k_dik), dbl_mb(k_dil), 
     $     dbl_mb(k_djk), dbl_mb(k_djl), dbl_mb(k_dkl), 
     $     dbl_mb(k_fij), dbl_mb(k_fik), dbl_mb(k_fil), 
     $     dbl_mb(k_fjk), dbl_mb(k_fjl), dbl_mb(k_fkl), maxd,
     $     tol2e, ocoul, oexch)
c     call ga_print(g_fock)
c     
      status = ma_pop_stack(l_fkl)
      status = ma_pop_stack(l_fjl)
      status = ma_pop_stack(l_fjk)
      status = ma_pop_stack(l_fil)
      status = ma_pop_stack(l_fik)
      status = ma_pop_stack(l_fij)
      status = ma_pop_stack(l_dkl)
      status = ma_pop_stack(l_djl)
      status = ma_pop_stack(l_djk)
      status = ma_pop_stack(l_dil)
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
      end
