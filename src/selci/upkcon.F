      subroutine selci_upkcon(norbs, iocc, mocc, nintpo, nbitpi)
#include "errquit.fh"
*
* $Id: upkcon.F,v 1.3 2003-10-17 22:58:17 carlfahl Exp $
*
      dimension iocc(norbs), mocc(nintpo)
#include "bitops.fh"
c
c     scalar version
c     unpack occupation iocc(*) bitwise from mocc, nbitpi bits
c     per integer word
c     iocc(1:norbs) = 0 (uocc), 1 (socc), 3 (docc)
c     mocc(1:nintpo)= packed representation of iocc, 2 bits per element
c     nintpo = no. of integers needed to pack the occupancy
c     nbitpi = no. of bits in an integer
c
c     no. of elements per integer
c
      nelpi = nbitpi / 2
      if (nintpo .ne. (norbs+nelpi-1)/nelpi) 
     &     call errquit('nintpo has strange value ',nintpo, INPUT_ERR)
c
      ilo = 1
      do 10 iword = 1, nintpo
         move = nbitpi - 2
         itemp = mocc(iword)
         ihi = min(ilo+nelpi-1,norbs)
         do 20 i = ilo,ihi
            iocc(i) = iand(ishft(itemp,-move),3)
            move = move - 2
 20      continue
         ilo = ihi + 1
 10   continue
c

      end
