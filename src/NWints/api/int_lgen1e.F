      subroutine int_lgen1e(i_basis, j_basis, ish, jsh, zerotol,
     &       ilab, jlab, l1e, Gen1e, lscr, scr, numgen) 
c $Id: int_lgen1e.F,v 1.5 1994-04-26 19:46:05 d3g681 Exp $
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
c
c::external subroutines used
c dcopy
c errquit
c::function
      logical int_chk_init
      external int_chk_init
      logical int_chk_sh
      external int_chk_sh
c::passed
      integer i_basis, j_basis, ish, jsh, l1e, numgen, lscr
      integer ilab(l1e),jlab(l1e)
      double precision zerotol, scr(lscr)
      double precision Gen1e(l1e)
c::local
      integer ibas, jbas, icount, i, j
      logical ijbas
      logical shells_ok
      logical nonzero
c
c check initialization and shells
c
      if (.not.int_chk_init('int_lgen1e'))
     &       call errquit('int_lgen1e: int_init was not called' ,0)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_lgen1e: invalid contraction/shell',0)
c      
      ibas  = i_basis + BASIS_HANDLE_OFFSET
      jbas  = j_basis + BASIS_HANDLE_OFFSET
      ijbas = ibas .eq. jbas
c
      call dcopy(l1e,Gen1e,1,scr,1)
c
      numgen = 0
      icount = 0
      do 00100 i = ibs_cn2bfr(1,ish,ibas), ibs_cn2bfr(2,ish,ibas)
        do 00200 j = ibs_cn2bfr(1,jsh,jbas), ibs_cn2bfr(2,jsh,jbas)
          icount = icount + 1
          nonzero = abs(scr(icount)).ge.zerotol
          if (nonzero) then
            if ((.not.ijbas).or.i.ge.j) then
              numgen = numgen + 1
              Gen1e(numgen) = scr(icount)
              ilab(numgen) = i
              jlab(numgen) = j
            endif
          endif
c
00200   continue
00100 continue
c
      end
