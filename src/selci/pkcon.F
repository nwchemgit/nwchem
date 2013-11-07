C> \ingroup selci
C> @{
C>
C> \brief Pack electron configurations
C>
C> Initially a configuration is stored in a long array of integers
C> where each integer corresponds to an orbital. The occupation of
C> the orbital is specified by the value of the corresponding integer:
C> - 0 (uocc) - is un-occupied
C> - 1 (socc) - is singly occupied
C> - 3 (docc) - is doubly occupied
C> As a result we need two bits per orbital to store a configuration.
C> This routine eliminates the superflous bits and compresses the
C> relevant data into a much more compact data structure `mocc`.
C>
      subroutine selci_pkcon(norbs, iocc, mocc, nintpo, nbitpi)
      implicit none
#include "errquit.fh"
*
* $Id$
*
      integer norbs        !< [Input] The number of orbitals
      integer nintpo       !< [Input] The number of integers per 
                           !< output configuration
      integer nbitpi       !< [Input] The number of bits per integer
      integer iocc(norbs)  !< [Input] The decompressed electron
                           !< configuration
      integer mocc(nintpo) !< [Output] The compressed electron 
                           !< configuration
#include "bitops.fh"
c
c     scalar version
c     pack occupation iocc(*) bitwise into mocc, nbitpi bits
c     per integer word
c     iocc(1:norbs) = 0 (uocc), 1 (socc), 3 (docc)
c     mocc(1:nintpo)= packed representation of iocc, 2 bits per element
c     nintpo = no. of integers needed to pack the occupancy
c     nbitpi = no. of bits in an integer
c
c     Local
c
      integer i, ilo, ihi, itemp, iword, move, nelpi
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
         itemp = 0
         ihi = min(ilo+nelpi-1,norbs)
         do 20 i = ilo,ihi
            itemp = ior(itemp,ishft(iocc(i),move))
            move = move - 2
 20      continue
         mocc(iword) = itemp
         ilo = ihi + 1
 10   continue
c
c      write(6,*) ' input occupation ',iocc
c      write(6,1) mocc
c 1    format(1x,5(z10,2x))
c
      end
C>
C> @}
