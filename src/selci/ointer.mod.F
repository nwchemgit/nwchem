      logical function selci_ointer(noconf,ioconf,joconf,nintpo)
*
* $Id: ointer.mod.F,v 1.2 1997-10-31 23:42:17 d3e129 Exp $
*
      integer ioconf(nintpo*noconf),joconf(nintpo)
c     
c     ivl set so can use short vector lengths ... if vector
c     length is too long we have redundant computation
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
c     strip mine into blocks so can call vpop32 instead of
c     selci_popcnt32(popcnt32 -> 100,000/s, vpop32 -> 540,000/s at vl=32
c     480,000/s at vl=16, 740,000/s at vl=1000)
c     
c     Search backwards since it is most likely to find something
c     that interacts near to the occupation that generated the
c     excitation being tested.
c     
      ointer = .true.
      nleft = noconf
      ioff  = max(0,noconf-ivl)*nintpo
c     
 100  ndo = min(ivl,nleft)
      if (nintpo .eq. 1) then
         itest = joconf(1)
         do i = 1,ndo
            iwrk1(i) = ieor(itest,ioconf(ioff+i))
         enddo
         call selci_vpop32(ndo,iwrk1,iwrk2)
         do i = 1, ndo
            if (iwrk2(i).le.4) return
         enddo
c$$$      else if (nintpo .eq. 2) then
c$$$         ndo2 = ndo + ndo
c$$$         itest1 = joconf(1)
c$$$         itest2 = joconf(2)
c$$$         do i = 1,ndo2,2
c$$$            iwrk1(i  ) = ieor(itest1,ioconf(ioff+i  ))
c$$$            iwrk1(i+1) = ieor(itest2,ioconf(ioff+i+1))
c$$$         enddo
c$$$         call selci_vpop32(ndo2,iwrk1,iwrk2)
c$$$         do i = 1, ndo2, 2
c$$$            if ((iwrk2(i)+iwrk2(i+1)).le.4) return
c$$$         enddo
      else
         k = 0
         do i = 1,ndo
            ioffk = ioff + k
            do iw = 1,nintpo
               iwrk1(k+iw) = ieor(joconf(iw),ioconf(ioffk+iw))
            enddo
            k = k + nintpo
         enddo
         call selci_vpop32(ndo*nintpo,iwrk1,iwrk2)
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
      ioff = max(0,ioff - ndo)*nintpo
      nleft = nleft - ndo
      if (nleft.gt.0) goto 100
c     
      ointer = .false.
c     
      end
