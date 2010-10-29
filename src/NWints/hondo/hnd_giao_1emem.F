c $Id$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : 1e- giao memory routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_giao_1emem(maxang)
      implicit none
#include "apiP.fh"
#include "hnd_pointers.fh"
c
      integer scrmem, maxbuf, maxang, i
c
c     set memory requirements
c
      maxbuf = 0
      isz_1e = max(isz_1e,maxbuf)
c
c     Scratch needed in various GIAO routines
c     hnd_giaoh01     9(L+1)(L+2)(L+2)
c     hnd_giaoh11     11(3)(L+2)(L+2)((2L+3)/2+1)
c     hnd_giaol10     9(L+1)(L+2)
c     hnd_giaos10     6(L+2)(L+1)
c     hnd_giaotv10    3(L+2)(L+3) + 3(L+2)(L+1)((2L+1)/2+1) +
c                     3(L+2)(L+3) + 3(L+2)(L+1)((2L+1)/2+1)
c
c Each of the five routines will use a subset of the pointers. Blocks
c for each pointer are setup using blocks of hnd_giaoh11, as these are
c the largest blocks. hnd_giaotv10 needs one extra block (which will
c be block 12), which lies half way one of the other (L+2)(L+3) blocks.
c Hence, for hnd_giaotv10 use gh01 ordering 1,12,2,3,...,11
c
      gh01(1) = 1                  
      do i = 2, 11
         gh01(i) = gh01(i-1) + 3*((2*maxang+3)/2+1)*(maxang+2)**2
      enddo
      scrmem = gh01(11) + 3*((2*maxang+3)/2+1)*(maxang+2)**2
      mem_1e = max(mem_1e,scrmem)
      mem_1e = max(mem_1e,isz_1e)
c
      gh01(12) = gh01(1) + (maxang+2)*(maxang+3)
c
      return
      end
