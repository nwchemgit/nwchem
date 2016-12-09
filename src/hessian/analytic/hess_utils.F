c
c     Wrapper to nga matrices
c
      subroutine hess_copymat(ndens,ndim1,ndim2,ndim3,g_a,g_b)
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "util.fh"
c
      integer ndens, ndim1, ndim2, ndim3, g_a, g_b
c
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer i
c
      do i = 1,ndens
        alo(1) = 1 + (i-1)*ndim1 ! n3xyz + 1
        ahi(1) = i * ndim1     ! 2 * n3xyz
        alo(2) = 1
        ahi(2) = ndim2 !nbf
        alo(3) = 1
        ahi(3) = ndim3 !nbf
        blo(1) = 1 + (i-1)*ndim1 ! n3xyz + 1
        bhi(1) = i * ndim1     ! 2 * n3xyz
        blo(2) = 1
        bhi(2) = ndim2 !nbf
        blo(3) = 1
        bhi(3) = ndim3 !nbf
        call nga_copy_patch('N', g_a, alo, ahi, g_b, blo, bhi)
      end do
c
      return
      end
