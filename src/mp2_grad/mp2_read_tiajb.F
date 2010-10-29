      subroutine mp2_read_tijab(nv_lo, nv_hi, irs, symia, 
     $     num_o, sym_hi_o, sym_lo_o, tunit, tunitptr, t)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "eaf.fh"
#include "util.fh"
      integer nv_lo, nv_hi
      integer irs(*), num_o(0:*), sym_hi_o(0:*), sym_lo_o(0:*), symia
      integer tunit
      double precision tunitptr
      double precision t(*)
c
c     Read t(j,b,i,a) all j, b for given i, a taking into
c     account symmetry blocking
c
      integer count, b, symb, symj
c
#include "bitops.fh"
c
      count = 0
      do b = nv_lo, nv_hi
         symb=irs(b)
         symj=ieor(symb,symia)
         if(num_o(symj).ne.0) count = count + 
     $              sym_hi_o(symj) - sym_lo_o(symj) + 1
      enddo
      if (count .gt. 0) then
         if (eaf_read(tunit, tunitptr, t, 8*count).ne.0) 
     $        call errquit('mp2_pijab: problem reading t_a',0, DISK_ERR)
         tunitptr = tunitptr + 8*count
         if (util_print('tdebug',print_debug)) then
            write(6,12) ddot(count, t, 1, t, 1)
 12         format(1x,12x,f20.8)
         endif
      endif
c
      end

