      subroutine int_l2e4c(bra_basis, ish, jsh, ket_basis, ksh, lsh,
     &       zerotol, canonicalize, leri, eri, nint, ilab, jlab, klab,
     &       llab, scr, lscr)
c $Id: int_l2e4c.F,v 1.2 1994-04-04 20:27:38 d3e129 Exp $
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c::externals used
c dcopy
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer bra_basis, ish, jsh, ket_basis, ksh, lsh, leri, lscr
      double precision zerotol, eri(leri), scr(lscr)
      integer ilab(leri), jlab(leri), klab(leri), llab(leri), nint
      logical canonicalize
c::local
      logical braket, allshl, ignore_ij_kl
      integer isym2
      integer i,j,k,l,bra,ket,icount,ij,kl
c
c:statement function
c
      isym2(i,j)=(max(i,j)*(max(i,j)-1))/2 + min(i,j)
c
c check initialization
c
      if (.not.int_chk_init('int_l2e4c'))
     &       call errquit('int_l2e4c: int_init was not called' ,0)

c
c copy eri's to scr
      call dcopy(leri,eri,1,scr,1)
c
      bra = bra_basis + BASIS_HANDLE_OFFSET
      ket = ket_basis + BASIS_HANDLE_OFFSET
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
              if (abs(scr(icount).gt.zerotol)) then
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
