      logical function selci_ointer(noconf,ioconf,joconf,nintpo,
     $     icon4, icon4new)
*
* $Id: ointer.parthack.F,v 1.2 1997-10-31 23:42:17 d3e129 Exp $
*
      integer ioconf(nintpo*noconf),joconf(nintpo)
      integer icon4(4,noconf), icon4new(4)
c     
c     icon4(1:4,*) are occupations of 4 orbital subspaces. 
c     First scan over these before doing a detailed scan of
c     what remains.
c
c     if vector length is too long we have redundant computation
c     (32 bit integer, 16 orbs per integer, 16 integers
c     needed for 255 orbs)
c     
      parameter (ivl=16)
      integer iwrk1(ivl*16),iwrk2(ivl*16), list(ivl)
      integer selci_popcnt32
      integer idiff(-10:10)
#include "bitops.fh"
      data idiff/10,9,8,7,6,5,4,3,2,1,0,1,2,3,4,5,6,7,8,9,10/
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
      selci_ointer = .true.
c
      ndo = 0
      do iconf = noconf, 1, -1
         itest = 
     $        idiff(icon4(1,iconf)-icon4new(1))  + 
     $        idiff(icon4(2,iconf)-icon4new(2))  + 
     $        idiff(icon4(3,iconf)-icon4new(3))  + 
     $        idiff(icon4(4,iconf)-icon4new(4))
         if (itest .le. 4) then
            ndo = ndo + 1
            list(ndo) = iconf
         endif
         if (ndo.eq.ivl .or. iconf.eq.1) then
            k = 0
            do ido = 1, ndo
               jconf = list(ido)
               joff  = (jconf-1)*nintpo
               do i = 1, nintpo
                  iwrk1(i+k) = ieor(ioconf(joff+i),joconf(i))
               enddo
               k = k + nintpo
            enddo
            call selci_vpop32(ndo*nintpo,iwrk1,iwrk2)
            k = 0
            do ido = 1, ndo
               isum = 0
               do i = 1, nintpo
                  isum = isum + iwrk2(k+i)
               enddo
               if (isum .le. 4) return
               k = k + nintpo
            enddo
            ndo = 0
         endif
      enddo
c     
      selci_ointer = .false.
c     
      end
