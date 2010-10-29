c $Id$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : 1e- derivative memory routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_d1emem(maxang)
      implicit none
#include "apiP.fh"
#include "hnd_pointers.fh"
c
      integer scrmem, maxbuf, maxang, ll
c
c     set memory requirements
c
      maxbuf = 0
      isz_1e = max(isz_1e,maxbuf)
c
c     Largest block of scratch from first derivative comes from 
c     hnd_tvd_ij.F:
c     3(L+3)(L+5)+3(L+3)^2+12(L+1)^2   for T
c     (L+1)(3(L+3)^2+6(L+1)^2)         for V
c
c     We are also setting up the pointers into the scratch block
c
      tvp (1) = 1                                 ! T (L+3)(L+5)
      tvp (2) = (maxang+3)*(maxang+5)+1           ! T (L+3)(L+5)
      tvp (3) = 2*tvp (2)+1                       ! T (L+3)(L+5)
      tvp (4) = 3*tvp (2)+1                       ! T (L+3)(L+3)
      tvp (5) = tvp (4)+(maxang+3)**2             ! T (L+3)(L+3)
      tvp (6) = tvp (5)+(maxang+3)**2             ! T (L+3)(L+3)
      tvp (7) = tvp (5)+(maxang+3)**2             ! T (L+1)(L+3)
      do ll = 8, 19
         tvp (ll) = tvp (ll-1)+(maxang+1)**2      ! T 11*(L+1)(L+3) + V (L+2)(L+3)^2
      enddo
      tvp (20) = tvp(19)+(maxang+3)**2*(maxang+2) ! V (L+2)(L+3)^2
      tvp (21) = tvp(20)+(maxang+3)**2*(maxang+2) ! V (L+2)(L+3)^2
      tvp (22) = tvp(21)+(maxang+3)**2*(maxang+2) ! V (L+2)(L+1)^2
      do ll = 23, 27                              ! V 5*(L+2)(L+1)^2 + T 2nd (L+1)^2
         tvp(ll) = tvp(ll-1)+(maxang+1)**2*(maxang+2)
      enddo
      do ll = 28, 54
         tvp (ll) = 0                             ! Pointers for 2nd der
      enddo 
      scrmem = tvp(27)+(maxang+1)**2*(maxang+2)
      mem_1e = max(mem_1e,scrmem)
c
      return
      end
