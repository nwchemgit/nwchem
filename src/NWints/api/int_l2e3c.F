      subroutine int_l2e3c(brain, ish, ketin, jsh, ksh,
     &       zerotol, canket, leri, eri, nint, ilab, jlab, klab,
     &       lscr, scr)
c $Id: int_l2e3c.F,v 1.2 1994-04-09 15:50:58 d3e129 Exp $
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
      integer brain, ish, jsh, ketin, ksh, leri, lscr
      double precision zerotol, eri(leri), scr(lscr)
      integer ilab(leri), jlab(leri), klab(leri), nint
      logical canket
c::local
      logical shells_ok
      integer i,j,k,bra,ket,icount
c
c check initialization and shells 
c
      if (.not.int_chk_init('int_l2e3c'))
     &       call errquit('int_l2e3c: int_init was not called' ,0)
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(ketin,jsh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,ksh)
      if (.not.shells_ok)
     &       call errquit('int_l2e3c: invalid contraction/shell',0)
c
c   check canonicalizations of input shells for canket
c      
      shells_ok = jsh.ge.ksh
      if (.not.shells_ok) then
        write(6,*)'int_2e3c: shells not canonical on input '
        write(6,*)'bra basis set handle:',brain
        write(6,*)'ket basis set handle:',ketin
        write(6,*)'                 ish:',ish
        write(6,*)'                 jsh:',jsh
        write(6,*)'                 ksh:',ksh
        call errquit('int_l2e3c: shells not in canonical order',0)
      endif
c
c copy eri's to scr
      call dcopy(leri,eri,1,scr,1)
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      nint = 0
      icount = 0
      do 00100 i = ibs_cn2bfr(1,ish,bra), ibs_cn2bfr(2,ish,bra)
        do 00200 j = ibs_cn2bfr(1,jsh,ket), ibs_cn2bfr(2,jsh,ket)
          do 00300 k = ibs_cn2bfr(1,ksh,ket), ibs_cn2bfr(2,ksh,ket)
            icount = icount + 1
            if (abs(scr(icount)).ge.zerotol) then
              if (canket)then
                if(j.ge.k) then
                  nint = nint + 1
                  eri(nint) = scr(icount)
                  ilab(nint) = i
                  jlab(nint) = j
                  klab(nint) = k
                endif
              else
                nint = nint + 1
                eri(nint) = scr(icount)
                ilab(nint) = i
                jlab(nint) = j
                klab(nint) = k
              endif
            endif
00300     continue
00200   continue
00100 continue
c
      end
