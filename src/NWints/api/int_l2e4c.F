      subroutine int_l2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &       zerotol, canonicalize, leri, eri, nint, ilab, jlab, klab,
     &       llab, lscr, scr)
c $Id: int_l2e4c.F,v 1.12 1995-12-22 01:39:07 d3e129 Exp $
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
c::external subroutines used
c errquit
c::functions
      logical  cando_sp, cando_nw, cando_txs
      logical  int_chk_sh
      external int_chk_sh
      external cando_sp, cando_nw, cando_txs
c::passed
      integer brain, ish, jsh, ketin, ksh, lsh, leri, lscr
      double precision zerotol, eri(leri), scr(lscr)
      integer ilab(leri), jlab(leri), klab(leri), llab(leri), nint
      logical canonicalize
c::local
      logical status_sp, status_txs
      logical braket, shells_ok
#if defined(INT_LIMIT_CHECK)
      logical allshl, ignore_ij_kl
#endif
*       integer ijsh, klsh
      integer isym2
      integer i,j,k,l,bra,ket,icount,ij,kl
      integer numint, newlscr
      logical dum_log_f, dum_log_t
      integer num_quart, intp 
      double precision roff(3), val, q4

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
*rak: no longer valid assumption :c
*rak: no longer valid assumption :c   check canonicalizations of input shells
*rak: no longer valid assumption :c      
*rak: no longer valid assumption :      shells_ok = ish.ge.jsh
*rak: no longer valid assumption :      shells_ok = shells_ok .and. ksh.ge.lsh
*rak: no longer valid assumption :      ijsh = isym2(ish,jsh)
*rak: no longer valid assumption :      klsh = isym2(ksh,lsh)
*rak: no longer valid assumption :      if (brain.eq.ketin) shells_ok = shells_ok .and. ijsh.ge.klsh
*rak: no longer valid assumption :      if (.not.shells_ok) then
*rak: no longer valid assumption :        write(6,*)'int_2e4c: shells not canonical on input '
*rak: no longer valid assumption :        write(6,*)'bra basis set handle:',brain
*rak: no longer valid assumption :        write(6,*)'ket basis set handle:',ketin
*rak: no longer valid assumption :        write(6,*)'                 ish:',ish
*rak: no longer valid assumption :        write(6,*)'                 jsh:',jsh
*rak: no longer valid assumption :        write(6,*)'                 ksh:',ksh
*rak: no longer valid assumption :        write(6,*)'                 lsh:',lsh
*rak: no longer valid assumption :        if (brain.eq.ketin) then
*rak: no longer valid assumption :          write(6,*)'                ijsh:',ijsh
*rak: no longer valid assumption :          write(6,*)'                klsh:',klsh
*rak: no longer valid assumption :        endif
*rak: no longer valid assumption :        call errquit('int_l2e4c: shells not in canonical order',0)
*rak: no longer valid assumption :      endif
c
      status_sp  = cando_sp(brain, ish, jsh)  .and.
     &    cando_sp(ketin, ksh, lsh)
      status_txs = cando_txs(brain, ish, jsh) .and.
     &    cando_txs(ketin, ksh, lsh)          .and.
     &    brain.eq.ketin                                 ! right now texas uses only one basis
c
      if (status_txs .and. (.not.status_sp)) then
        roff(1) = 0.0d00
        roff(2) = 0.0d00
        roff(3) = 0.0d00
        num_quart = 1
        dum_log_f = .false.
        dum_log_t = .true.
        q4 = 1.0d00
        call texas_hf2_m (brain, ish, jsh, ketin, ksh, lsh, num_quart,
     &      q4, dum_log_f,
     &      roff, roff, roff, roff, dum_log_f,
     &      eri, leri,
     &      ilab, jlab, klab, llab, nint, dum_log_t,
     &      dum_log_f, scr, lscr)
#if defined(INT_LIMIT_CHECK)
        braket = brain.eq.ketin
        allshl = isym2(ish,jsh) .eq. isym2(ksh,lsh)
        ignore_ij_kl = (.not.braket) .or. (.not.allshl)
        icount = 0
        do 10100 intp = 1,nint
          i = ilab(intp)
          j = jlab(intp)
          k = klab(intp)
          l = llab(intp)
          ij = isym2(i,j)
          kl = isym2(k,l)
          val = eri(intp)
          if (abs(val).ge.zerotol) then
            if (i.ge.j.and.k.ge.l.and.
     &          (ignore_ij_kl.or.(ij.ge.kl)))then
              icount = icount + 1
              if (icount.lt.intp) then
                eri(icount)  = val
                if (canonicalize.and.braket) then
                  call int_canon(i,j,k,l,
     &                ilab(icount), jlab(icount),
     &                klab(icount), llab(icount))
                else
                  ilab(icount) = i
                  jlab(icount) = j
                  klab(icount) = k
                  llab(icount) = l
                endif
              else if (canonicalize.and.braket) then
                call int_canon(i,j,k,l,
     &              ilab(intp), jlab(intp),
     &              klab(intp), llab(intp))
              endif
            endif
          endif
10100   continue
        nint = icount
#else
        braket = brain.eq.ketin
        icount = 0
        do 10100 intp = 1,nint
          i = ilab(intp)
          j = jlab(intp)
          k = klab(intp)
          l = llab(intp)
          ij = isym2(i,j)
          kl = isym2(k,l)
          val = eri(intp)
          if (abs(val).ge.zerotol) then
            icount = icount + 1
            eri(icount)  = val
            if (canonicalize.and.braket) then
              call int_canon(i,j,k,l,
     &            ilab(icount), jlab(icount),
     &            klab(icount), llab(icount))
            else if (canonicalize) then
              if (i.lt.j) then
                ilab(icount) = j
                ilab(icount) = i
              else
                ilab(icount) = i
                jlab(icount) = j
              endif
              if (k.lt.l) then
                klab(icount) = l
                llab(icount) = k
              else
                klab(icount) = k
                llab(icount) = l
              endif
            else
              ilab(icount) = i
              jlab(icount) = j
              klab(icount) = k
              llab(icount) = l
            endif
          endif
10100   continue
        nint = icount
#endif
      else
        
c
c compute eri's (instead of copy) in front of scr()
c
        bra = brain + BASIS_HANDLE_OFFSET
        ket = ketin + BASIS_HANDLE_OFFSET
        
        numint =         ibs_cn2bfr(2,ish,bra)-ibs_cn2bfr(1,ish,bra)+1
        numint = numint*(ibs_cn2bfr(2,jsh,bra)-ibs_cn2bfr(1,jsh,bra)+1)
        numint = numint*(ibs_cn2bfr(2,ksh,ket)-ibs_cn2bfr(1,ksh,ket)+1)
        numint = numint*(ibs_cn2bfr(2,lsh,ket)-ibs_cn2bfr(1,lsh,ket)+1)
        newlscr = lscr - numint
        call int_2e4c(brain,ish,jsh,ketin,ksh,lsh,
     &      newlscr,scr(numint+1),numint,scr)
c
#if defined(INT_LIMIT_CHECK)
c
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
     &                (ignore_ij_kl.or.(ij.ge.kl)))then
                    nint = nint + 1
                    eri(nint) = scr(icount)
                    if (canonicalize.and.braket) then
                      call int_canon(i,j,k,l,
     &                    ilab(nint),jlab(nint),
     &                    klab(nint),llab(nint))
                    else
                      ilab(nint) = i
                      jlab(nint) = j
                      klab(nint) = k
                      llab(nint) = l
                    endif
                  endif
                endif
00400         continue
00300       continue
00200     continue
00100   continue
#else
        braket = bra .eq. ket
        nint = 0
        icount = 0
        do 00100 i = ibs_cn2bfr(1,ish,bra), ibs_cn2bfr(2,ish,bra)
          do 00200 j = ibs_cn2bfr(1,jsh,bra), ibs_cn2bfr(2,jsh,bra)
            do 00300 k = ibs_cn2bfr(1,ksh,ket), ibs_cn2bfr(2,ksh,ket)
              do 00400 l = ibs_cn2bfr(1,lsh,ket), ibs_cn2bfr(2,lsh,ket)
                icount = icount + 1
                if (abs(scr(icount)).ge.zerotol) then
                  nint = nint + 1
                  eri(nint) = scr(icount)
                  if (canonicalize.and.braket) then
                      call int_canon(i,j,k,l,
     &                    ilab(nint),jlab(nint),
     &                    klab(nint),llab(nint))
                  elseif (canonicalize) then
                    if (i.lt.j) then
                      ilab(nint) = j
                      ilab(nint) = i
                    else
                      ilab(nint) = i
                      jlab(nint) = j
                    endif
                    if (k.lt.l) then
                      klab(nint) = l
                      llab(nint) = k
                    else
                      klab(nint) = k
                      llab(nint) = l
                    endif
                  else
                    ilab(nint) = i
                    jlab(nint) = j
                    klab(nint) = k
                    llab(nint) = l
                  endif
                endif
00400         continue
00300       continue
00200     continue
00100   continue
#endif
c
      endif
      end
