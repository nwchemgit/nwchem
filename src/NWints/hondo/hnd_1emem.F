cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : 1e- memory routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_1emem
      implicit none
#include "apiP.fh"
      integer scrmem, maxbuf
c
c     set memory requirements
c     The HONDO integrals assume that the nw code sets up the
c     buffer size appropriately.
c     The scrmem size is the size of a g-g shell of integrals
c     (the maximum that the code needs since it can only do up to g).
c
      scrmem = 255
      maxbuf = 0
c
      isz_1e = max(isz_1e,maxbuf)
      mem_1e = max(mem_1e,scrmem) 
c
      return
      end
