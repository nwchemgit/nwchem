      subroutine int_l2e2c(brain, ish, ketin, jsh, 
     &       zerotol, leri, eri, nint, ilab, jlab, 
     &       lscr, scr)
c $Id: int_l2e2c.F,v 1.5 1995-05-09 18:00:55 d3e129 Exp $
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
      logical int_chk_sh
      external int_chk_sh
c::passed
      integer brain, ish, jsh, ketin, leri, lscr
      double precision zerotol, eri(leri), scr(lscr)
      integer ilab(leri), jlab(leri), nint
c::local
      logical shells_ok
      integer i,j,bra,ket,icount
      integer numint, newlscr
c
c:statement function
c
c
c check shells 
c
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
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
c copy eri's to scr
*      call dcopy(leri,eri,1,scr,1)
*
c compute eris into scr instead of copying them
c
c 
      numint =         ibs_cn2bfr(2,ish,bra)-ibs_cn2bfr(1,ish,bra)+1
      numint = numint*(ibs_cn2bfr(2,jsh,ket)-ibs_cn2bfr(1,jsh,ket)+1)
      newlscr = lscr - numint
      call int_2e2c(brain,ish,ketin,jsh,
     &    newlscr,scr(numint+1),numint,scr)
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
