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
c
      scrmem = 0
      maxbuf = 0
c
      isz_1e = max(isz_1e,scrmem)
      mem_1e = max(mem_1e,maxbuf) 
c
      return
      end
