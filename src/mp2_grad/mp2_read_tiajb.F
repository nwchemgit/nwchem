      subroutine mp2_read_tijab(nv_lo, nv_hi, irs, symia, 
     $     num_o, sym_hi_o, sym_lo_o, tunit, tunitptr, t)
      implicit none
#include "eaf.fh"
      integer nv_lo, nv_hi
      integer irs(*), num_o(0:*), sym_hi_o(0:*), sym_lo_o(0:*), symia
      integer tunit
      double precision tunitptr
      double precision t(*)
c
      integer count, b, symb, symj
c
      count = 0
      do b = nv_lo, nv_hi
         symb=irs(b)
         symj=ieor(symb,symia)
         if(num_o(symj).ne.0) count = count + 
     $              sym_hi_o(symj) - sym_lo_o(symj) + 1
      enddo
      if (eaf_read(tunit, tunitptr, t, 8*count).ne.0) 
     $     call errquit('mp2_pijab: problem reading t_a',0)
      tunitptr = tunitptr + 8*count
c
      end

