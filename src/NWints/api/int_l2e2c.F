      subroutine int_l2e2c(brain, ish, ketin, jsh, 
     &       zerotol, leri, eri, nint, ilab, jlab, 
     &       lscr, scr)
c $Id: int_l2e2c.F,v 1.2 1994-04-26 19:40:24 d3g681 Exp $
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
c::external subroutines used
c dcopy
c errquit
c::functions
      logical int_chk_init
      external int_chk_init
      logical int_chk_sh
      external int_chk_sh
c::passed
      integer brain, ish, jsh, ketin, leri, lscr
      double precision zerotol, eri(leri), scr(lscr)
      integer ilab(leri), jlab(leri), nint
c::local
      logical shells_ok
      integer i,j,bra,ket,icount
c
c:statement function
c
c
c check initialization and shells 
c
      if (.not.int_chk_init('int_l2e2c'))
     &       call errquit('int_l2e2c: int_init was not called' ,0)
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(ketin,jsh)
      if (.not.shells_ok)
     &       call errquit('int_l2e2c: invalid contraction/shell',0)
c
c   check canonicalizations of input shells for canket
c      
      if (brain.eq.ketin) then
        shells_ok = ish.ge.jsh
      else
        shells_ok = .true.
      endif
      if (.not.shells_ok) then
        write(6,*)'int_2e2c: shells not canonical on input '
        write(6,*)'bra basis set handle:',brain
        write(6,*)'ket basis set handle:',ketin
        write(6,*)'                 ish:',ish
        write(6,*)'                 jsh:',jsh
        call errquit('int_l2e2c: shells not in canonical order',0)
      endif
c
c copy eri's to scr
      call scopy(leri,eri,1,scr,1)
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      nint = 0
      icount = 0
      do 00100 i = ibs_cn2bfr(1,ish,bra), ibs_cn2bfr(2,ish,bra)
        do 00200 j = ibs_cn2bfr(1,jsh,ket), ibs_cn2bfr(2,jsh,ket)
          icount = icount + 1
          if (abs(scr(icount)).ge.zerotol) then
            if (bra.eq.ket) then
              if (i.ge.j) then
                nint = nint + 1
                eri(nint) = scr(icount)
                ilab(nint) = i
                jlab(nint) = j
              endif
            else
              nint = nint + 1
              eri(nint) = scr(icount)
              ilab(nint) = i
              jlab(nint) = j
            endif
          endif
00200   continue
00100 continue
c
      end
