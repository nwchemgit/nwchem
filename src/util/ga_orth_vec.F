      subroutine ga_orthog_vec(n, nvec, g_m, g_x, j)
C$Id$
      implicit none
#include "global.fh"
      integer n                 ! vector length
      integer nvec              ! no. of vectors
      integer g_m               ! GA handle for matrix
      integer g_x               ! GA handle for vector
      integer j                 ! Column for vector
c
c     orthogonalize the vector x(1:n,j) to the vectors g(1:n, 1:nvec)
c
c     Note that x is NOT normalized.
c
      integer i, iter
      double precision xm, mm
c
      do iter = 1, 2
         do i = 1, nvec
            xm = ga_ddot_patch(
     $           g_x, 'n', 1, n, j, j,
     $           g_m, 'n', 1, n, i, i)
            mm = ga_ddot_patch(
     $           g_m, 'n', 1, n, i, i, 
     $           g_m, 'n', 1, n, i, i)
            call ga_dadd_patch(
     $           1.0d0, g_x, 1, n, j, j,
     $           -xm/mm, g_m, 1, n, i, i,
     $           g_x, 1, n, 1, 1)
         enddo
      enddo
c
      end

