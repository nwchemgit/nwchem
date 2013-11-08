C> \file upkcon.F
C> Unpack configuration
C>
C> \ingroup selci
C> @{
C>
C> \brief Unpack a configuration
C>
C> Unpacks a configuration and stores the orbital occupations in an
C> integer array `iocc`.
C>
      subroutine selci_upkcon(norbs, iocc, mocc, nintpo, nbitpi)
#include "errquit.fh"
*
* $Id$
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
C>
C> @}
