      subroutine int_l2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &       zerotol, canonicalize, leri, eri, nint, ilab, jlab, klab,
     &       llab, lscr, scr)
c $Id: int_l2e4c.F,v 1.8 1994-09-22 13:34:12 d3e129 Exp $
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
      integer brain, ish, jsh, ketin, ksh, lsh, leri, lscr
      double precision zerotol, eri(leri), scr(lscr)
      integer ilab(leri), jlab(leri), klab(leri), llab(leri), nint
      logical canonicalize
c::local
      logical braket, allshl, ignore_ij_kl, shells_ok
      integer ijsh, klsh
      integer isym2
      integer i,j,k,l,bra,ket,icount,ij,kl
c
c:statement function
c
      isym2(i,j)=(max(i,j)*(max(i,j)-1))/2 + min(i,j)
c
c check initialization and shells 
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(brain,jsh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,ksh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,lsh)
      if (.not.shells_ok)
     &       call errquit('int_l2e4c: invalid contraction/shell',0)
c
c   check canonicalizations of input shells
c      
      shells_ok = ish.ge.jsh
      shells_ok = shells_ok .and. ksh.ge.lsh
      ijsh = isym2(ish,jsh)
      klsh = isym2(ksh,lsh)
      if (brain.eq.ketin) shells_ok = shells_ok .and. ijsh.ge.klsh
      if (.not.shells_ok) then
        write(6,*)'int_2e4c: shells not canonical on input '
        write(6,*)'bra basis set handle:',brain
        write(6,*)'ket basis set handle:',ketin
        write(6,*)'                 ish:',ish
        write(6,*)'                 jsh:',jsh
        write(6,*)'                 ksh:',ksh
        write(6,*)'                 lsh:',lsh
        if (brain.eq.ketin) then
          write(6,*)'                ijsh:',ijsh
          write(6,*)'                klsh:',klsh
        endif
        call errquit('int_l2e4c: shells not in canonical order',0)
      endif
c
c copy eri's to scr
      call dcopy(leri,eri,1,scr,1)
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      braket = bra .eq. ket
      allshl = isym2(ish,jsh) .eq. isym2(ksh,lsh)
      ignore_ij_kl = (.not.braket) .or. (.not.allshl)

c
      nint = 0
      icount = 0
      do 00100 i = ibs_cn2bfr(1,ish,bra), ibs_cn2bfr(2,ish,bra)
        do 00200 j = ibs_cn2bfr(1,jsh,bra), ibs_cn2bfr(2,jsh,bra)
          ij = isym2(i,j)
          do 00300 k = ibs_cn2bfr(1,ksh,ket), ibs_cn2bfr(2,ksh,ket)
            do 00400 l = ibs_cn2bfr(1,lsh,ket), ibs_cn2bfr(2,lsh,ket)
              kl = isym2(k,l)
              icount = icount + 1
              if (abs(scr(icount)).ge.zerotol) then
                if (i.ge.j.and.k.ge.l.and.
     &                 (ignore_ij_kl.or.(ij.ge.kl)))then
                  nint = nint + 1
                  eri(nint) = scr(icount)
                  if (canonicalize.and.braket) then
                    call int_canon(i,j,k,l,
     &                     ilab(nint),jlab(nint),
     &                     klab(nint),llab(nint))
                  else
                    ilab(nint) = i
                    jlab(nint) = j
                    klab(nint) = k
                    llab(nint) = l
                  endif
                endif
              endif
00400       continue
00300     continue
00200   continue
00100 continue
c
      end
