      subroutine int_l2e4c(bra_basis, ish, jsh, ket_basis, ksh, lsh,
     &       zerotol, canonicalize, leri, eri, nint, ilab, jlab, klab,
     &       llab, scr, lscr)
      implicit none
#include "apiP.fh"
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c::externals used
c dcopy
c::functions
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
                    call canonical(i,j,k,l,
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
      subroutine canonical(i, j, k, l, ii, jj, kk, ll)
      implicit none
      integer i,j,k,l,ii,jj,kk,ll,iii,jjj
c
c     return in ii,jj,kk,ll the cannonically ordered i,j,k,l
c
c     (i>=j) >= (k>=l)
c
      if (i .ge. j) then
         ii = i
         jj = j
      else
         ii = j
         jj = i
      endif
      if (k .ge. l) then
         kk = k
         ll = l
      else
         kk = l
         ll = k
      endif
      if (ii.lt.kk .or. (ii.eq.k .and. jj.lt.ll)) then
         iii = ii
         jjj = jj
         ii = kk
         jj = ll
         kk = iii
         ll = jjj
      endif
c
      end
