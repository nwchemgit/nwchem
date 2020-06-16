C> \ingroup selci
C> @{
      logical function selci_ointer(noconf,ioconf,joconf,nintpo)
*
* $Id$
*
      integer ioconf(nintpo*noconf),joconf(nintpo)
c     
c     if vector length is too long we have redundant computation
c     (32 bit integer, 16 orbs per integer, 16 integers
c     needed for 255 orbs)
c     
      parameter (ivl=16)
      integer iwrk1(ivl*16),iwrk2(ivl*16)
#include "bitops.fh"
c     
c     return true if orbital configuration joconf can interact
c     with any of the configurations in ioconf. false otherwise.
c     can interact if orbital occupation difference is <= 4
c     
c     use ieor and popcnt to compute the occupation difference
c     between the I and J occupations represented as bit patterns
c     
c     Search backwards since it is most likely to find something
c     that interacts near to the occupation that generated the
c     excitation being tested.
c     
      selci_ointer = .true.
      nleft = noconf
      ioff  = max(0,noconf-ivl)*nintpo
c     
 100  ndo = min(ivl,nleft)
      if (nintpo .eq. 1) then
         itest = joconf(1)
         do i = 1,ndo
            iwrk1(i) = ieor(itest,ioconf(ioff+i))
         enddo
         do i = 1,ndo
           iwrk2(i) = popcnt(iwrk1(i))
         enddo
         do i = 1, ndo
            if (iwrk2(i).le.4) return
         enddo
      else if (nintpo .eq. 2) then
         ndo2 = ndo + ndo
         itest1 = joconf(1)
         itest2 = joconf(2)
         do i = 1,ndo2,2
            iwrk1(i  ) = ieor(itest1,ioconf(ioff+i  ))
            iwrk1(i+1) = ieor(itest2,ioconf(ioff+i+1))
         enddo
         do i = 1,ndo2
           iwrk2(i) = popcnt(iwrk1(i))
         enddo
         do i = 1, ndo2, 2
            if ((iwrk2(i)+iwrk2(i+1)).le.4) return
         enddo
      else if (nintpo .eq. 3) then
         ndo3 = ndo + ndo + ndo
         itest1 = joconf(1)
         itest2 = joconf(2)
         itest3 = joconf(3)
         do i = 1,ndo3,3
            iwrk1(i  ) = ieor(itest1,ioconf(ioff+i  ))
            iwrk1(i+1) = ieor(itest2,ioconf(ioff+i+1))
            iwrk1(i+2) = ieor(itest3,ioconf(ioff+i+2))
         enddo
         do i = 1,ndo3
           iwrk2(i) = popcnt(iwrk1(i))
         enddo
         do i = 1, ndo3,3
            if ((iwrk2(i)+iwrk2(i+1)+iwrk2(i+2)).le.4) return
         enddo
      else
         k = 0
         do i = 1,ndo
            ioffk = ioff + k
            do iw = 1,nintpo
               iwrk1(k+iw) = ieor(joconf(iw),ioconf(ioffk+iw))
            enddo
            k = k + nintpo
         enddo
         do i = 1,ndo*nintpo
           iwrk2(i) = popcnt(iwrk1(i))
         enddo
         k = 0
         do i = 1, ndo
            isum = 0
            do iw = 1, nintpo
               isum = isum + iwrk2(k+iw)
            enddo
            if (isum .le. 4) return
            k = k + nintpo
         enddo
      endif
c     
      ioff = max(0,ioff - ndo*nintpo)
      nleft = nleft - ndo
      if (nleft.gt.0) goto 100
c     
      selci_ointer = .false.
c     
      end
C> @}
