#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on 12/ 8/97 at 23:34:29
C****************************************************************
C
      subroutine drftb1 (n,c,ch,wa,ifac)
*      
* $Id$
*      
*      
* $Id$
*      
      double precision c(1), ch(1), wa(1)
      integer ifac(*)
c      
      integer i1
      doubleprecision d1, d2, d3, d4
      logical l3, l11, l12
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1 = 1, nf
         l11 = na .ne. 0
         ip = ifac(k1+2)
         l12 = ip .ne. 4
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip .ne. 4) go to 103
         ix2 = iw + ido
         ix3 = ix2 + ido
         if (na .ne. 0) go to 101
         call dradb4 (ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
         go to 102
  101    continue
         call dradb4 (ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
  102    continue
         na = 1 - na
         go to 115
c
  103    continue
         if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call dradb2 (ido, l1, c, ch, wa(iw))
         go to 105
  104    continue
         call dradb2 (ido, l1, ch, c, wa(iw))
  105    continue
         na = 1 - na
         go to 115
c
  106    continue
         if (ip .ne. 3) go to 109
         ix2 = iw + ido
         if (na .ne. 0) go to 107
         call dradb3 (ido, l1, c, ch, wa(iw), wa(ix2))
         go to 108
  107    continue
         call dradb3 (ido, l1, ch, c, wa(iw), wa(ix2))
  108    continue
         na = 1 - na
         go to 115
c
  109    continue
         if (ip .ne. 5) go to 112
         ix2 = iw + ido
         ix3 = ix2 + ido
         ix4 = ix3 + ido
         if (na .ne. 0) go to 110
         call dradb5 (ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
         go to 111
  110    continue
         call dradb5 (ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
  111    continue
         na = 1 - na
         go to 115
c
  112    continue
         if (na .ne. 0) go to 113
         call dradbg (ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
         go to 114
  113    continue
         call dradbg (ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
         l3 = ido .eq. 1
  114    continue
         if (ido .eq. 1) na = 1 - na
  115    continue
         l1 = l2
         iw = iw + (ip - 1)*ido
  116 continue
      if (na .eq. 0) return
      i1 = iand(max0(n,0),3)
      do i = 1, i1
         c(i) = ch(i)
      end do
      do 117 i = i1 + 1, n, 4
         d1 = ch(i)
         d2 = ch(i+1)
         d3 = ch(i+2)
         d4 = ch(i+3)
         c(i) = d1
         c(i+1) = d2
         c(i+2) = d3
         c(i+3) = d4
  117 continue
c      
      return
      end
