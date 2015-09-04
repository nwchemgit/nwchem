      subroutine pphilo(id, tot, lo, up, num)
c
c     Subroutine to determine the load balancing for manually
c     parallaleized loops.
c      
      implicit none
#include "global.fh"
      integer id     ! Node ID
      integer tot    ! Total number of elements
      integer lo     ! Lower bound
      integer up     ! Upper bound
      integer num    ! Number of elements on this node
      integer nnodes ! Number of nodes
c
      nnodes = ga_nnodes()
      if (nnodes > 1) then
        lo = 1 + (id * tot ) / nnodes
        up = ((id +1) * tot) / nnodes
        num = up - lo + 1
      else
        lo = 1
        up = tot
        num = tot
      end if


      end subroutine pphilo
