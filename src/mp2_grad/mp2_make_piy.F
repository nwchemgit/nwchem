      subroutine mp2_make_piy(no_lo, no_hi, g_pij, g_liy, evals)
*
* $Id$
*
      implicit none
#include "global.fh"
      integer no_lo, no_hi      ! [input] Range of active occupied orbitals
      integer g_pij             ! [input/output] Core+active density
      integer g_liy             ! [input] Core-Active Lagrangian
      double precision evals(*) ! [input] Eigenvalues
c
c     Piy = 1/2 L1iy (ei - ey)^-1
c
      integer i, y
      double precision liy
c
      call ga_sync()
      do y = 1+ga_nodeid(), no_lo-1, ga_nnodes() ! Should use locality of y
         do i = no_lo, no_hi
            call ga_get(g_liy, i-no_lo+1, i-no_lo+1, y, y, liy, 1)
            liy = -liy  / (evals(i)-evals(y))
            call ga_put(g_pij, i, i, y, y, liy, 1)
            call ga_put(g_pij, y, y, i, i, liy, 1)
         enddo
      enddo
c
      end
      subroutine mp2_make_wiy(no_lo, no_hi, g_wij, g_liy, evals)
      implicit none
#include "global.fh"
      integer no_lo, no_hi      ! [input] Range of active occupied orbitals
      integer g_wij             ! [input/output] Core+active weighted density
      integer g_liy             ! [input] Core-Active Lagrangian
      double precision evals(*) ! [input] Eigenvalues
c
c     Wiy = 1/2 L1iy 
c
      integer i, y
      double precision liy
c
      call ga_sync()
      do y = 1+ga_nodeid(), no_lo-1, ga_nnodes() ! Should use locality of y
         do i = no_lo, no_hi
            call ga_get(g_liy, i-no_lo+1, i-no_lo+1, y, y, liy, 1)
            liy = 0.5d0*liy
            call ga_put(g_wij, i, i, y, y, liy, 1)
            call ga_put(g_wij, y, y, i, i, liy, 1)
         enddo
      enddo
c
      end
