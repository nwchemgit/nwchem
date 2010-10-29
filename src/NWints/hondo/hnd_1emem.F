c $Id$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : 1e- memory routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_1emem(maxang)
      implicit none
#include "apiP.fh"
#include "hnd_pointers.fh"
c
      integer scrmem, maxbuf, maxang
c
c     set memory requirements
c
      maxbuf = 0
      isz_1e = max(isz_1e,maxbuf)
c
c     Scratch needed in hnd_stv_int.F:
c     3(L+1)(L+3)            for S
c     scr for S + 3(L+1)^2   for T
c     3(L+1)^3               for V
c
      stp(1) = 1                               ! S
      stp(2) = (maxang+1)*(maxang+3)+1         ! S
      stp(3) = 2*stp (2)+1                     ! S
      stp(4) = 3*stp (2)+1                     ! T
      stp(5) = stp (4)+(maxang+1)**2           ! T
      stp(6) = stp (5)+(maxang+1)**2           ! T
      stp(7) = stp (6)+(maxang+1)**2           ! V
      stp(8) = stp (7)+(maxang+1)**3           ! V
      stp(9) = stp (8)+(maxang+1)**3           ! V
      scrmem = stp (9)+(maxang+1)**3  
      mem_1e = max(mem_1e,scrmem)
c
      return
      end
