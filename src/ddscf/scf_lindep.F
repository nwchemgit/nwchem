      subroutine scf_lindep(rtdb, geom, basis, g_u, n_indep)
c     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "util.fh"
      integer rtdb, geom, basis ! [input] Handles
      integer g_u              ! [input] GA handle for output vectors
      integer n_indep          ! [output] No. of linearly indep vectors
c
c     Return in g_u the linearly independent eigenvectors of the overlap
c     matrix in the given basis.  Note that these vectors are orthonormal
c     over the standard inner product, but are NOT orthonormal over the
c     overlap matrix.
c
c     n_indep returns the number of vectors
c
c     To project a set of MOs C into the linearly independent space 
c     use U * UT * C.
c
      integer n_dep            ! No. of dependent vectors
      double precision tol_default, tol
      parameter (tol_default = 1d-5)
c
      logical oprint
      integer g_over
      integer l_ev, k_ev
      integer i, j, nbf
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      oprint = util_print('lindep', print_default) .and. 
     $      ga_nodeid().eq.0
c
c     Make and diag the overlap matrix
c
       if (.not. bas_numbf(basis, nbf)) 
     $      call errquit('scf_lin: bad basis?',0, BASIS_ERR)
       if (.not. ma_push_get(MT_DBL, nbf, 'Evals', l_ev, k_ev))
     $      call errquit('scf_lindep: push of evals failed', nbf,
     &       MA_ERR)
       g_over  = ga_create_atom_blocked(geom, basis,'Temp Over')
c
       call ga_zero(g_over)
       call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
       if (util_print('ao overlap', print_debug))
     $      call ga_print(g_over)
c
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
       call ga_pdsyev(g_over, g_u, dbl_mb(k_ev), 0)
#else
       call ga_diag_std(g_over, g_u, dbl_mb(k_ev))
#endif
#else
       call ga_diag_std_seq(g_over, g_u, dbl_mb(k_ev))
#endif
c
c     Get user input parameters ... only one allowed
c
c     lindep:n_dep  ... force number of linearly dep vectors vectors
c     lindep:tol    ... tolerance 
c
       if (rtdb_get(rtdb, 'lindep:n_dep', mt_int, 1, n_dep)) then
          if (oprint) then
             write(6, 2) n_dep
 2           format(/' Forcing the number of linearly dependent ',
     $            'vectors to be', i3)
             call util_flush(6)
          endif
       else
          n_dep = -1
       endif
       if (rtdb_get(rtdb, 'lindep:tol', mt_dbl, 1, tol)) then
          if (oprint) then
             write(6, 3) tol
 3           format(/' Setting tolerance for linear dependence to',
     $            1p,d9.2)
             call util_flush(6)
          endif
       else
          tol = -1.0d0
       endif
c
       if ((tol.ne.-1.0d0) .and. (n_dep.ne.-1)) call errquit
     $      ('scf_lindep: only one of n_dep and tol permitted',0,
     &       UNKNOWN_ERR)
c
c     Determine linearly independent set
c
       if (n_dep .eq. -1) then
          if (tol .eq. -1.0d0) tol = tol_default
          do i = 0, nbf-1
             if (dbl_mb(k_ev+i).gt.tol) goto 333
          enddo
 333      n_dep = i
       endif
c
       if (n_dep.gt.0 .and. oprint .and. ga_nodeid().eq.0) then
          write(6,334) n_dep
 334      format(/' !! The overlap matrix has ',i3,' vectors deemed ',
     $         'linearly dependent with'/'    eigenvalues:')
          write(6,33) (dbl_mb(k_ev+j),j=0,n_dep-1)
 33       format(1p,8d9.2)
          write(6,*)
          call util_flush(6)
       endif
c
c     Move the linearly independent vectors down to the beginning
c     of the vector set
c
       if (n_dep.gt.0 .and. ga_nodeid().eq.0) then
          do i = 1, n_dep
             j = nbf - i + 1
             call ga_get(g_u, 1, nbf, j, j, dbl_mb(k_ev), nbf)
             call ga_put(g_u, 1, nbf, i, i, dbl_mb(k_ev), nbf)
          enddo
       endif
       call ga_sync
c
       n_indep = nbf - n_dep
c
       if (.not. ga_destroy(g_over)) call errquit
     $      ('scf_lindep: destroy of over?', 0, GA_ERR)
       if (.not. ma_pop_stack(l_ev)) call errquit
     $      ('scf_lindep: destroy of ev?', 0, MA_ERR)
c
       end
