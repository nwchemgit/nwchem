#if defined(XLF14) || defined(XLFLINUX)
@PROCESS OPT(0) STRICT(ALL)
#endif
      subroutine dcfti1 (n,wa,ifac)
*
* $Id$
*
      double precision wa(1), arg, argh, argld, fi, tpi
      integer ifac(*), ntryh(4)
      data ntryh(1), ntryh(2), ntryh(3), ntryh(4) /3, 4, 2, 5/
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c
      nl = n
      nf = 0
      j = 0
c
  101 j = j+1
      if (j.le.4) ntry = ntryh(j)
      if (j.gt.4) ntry = ntry + 2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr.ne.0) go to 101
c
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
c
  107 if (nl .ne. 1) go to 104
c
      ifac(1) = n
      ifac(2) = nf
c
      argh = tpi/dble(n)
      i = 2
      l1 = 1
      do 110 k1=1,nf
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         idot = ido+ido+2
         ipm = ip-1
c
         do 109 j=1,ipm
            i1 = i
            wa(i-1) = 1.d0
            wa(i) = 0.d0
            ld = ld+l1
            fi = 0.d0
            argld = dble(ld)*argh
            do 108 ii=4,idot,2
               i = i+2
               fi = fi+1.d0
               arg = fi*argld
               wa(i-1) = dcos(arg)
               wa(i) = dsin(arg)
  108       continue
            if (ip .le. 5) go to 109
            wa(i1-1) = wa(i-1)
            wa(i1) = wa(i)
  109    continue
c
         l1 = l2
  110 continue
c
      return
      end
