      integer function selci_iodiff(ioconf,joconf,nintpo)
      integer ioconf(*),joconf(*), selci_popcnt32
#include "bitops.fh"
c
c     use ieor and popcnt to compute the occupation difference
c     between the I and J occupations represented as bit patterns
c
      ipop = 0
      do 10 i = 1,nintpo
         ipop = ipop + selci_popcnt32(ieor(ioconf(i),joconf(i)))
 10   continue
c
      selci_iodiff = ipop
      end
      subroutine selci_viodiff(ioconf,joconf,nintpo,nconf,iexcit)
      integer ioconf(nintpo,nconf), iexcit(nconf)
      integer joconf(*)
#include "bitops.fh"
c
      parameter (maxtmp=256)
      integer tmp(maxtmp), out(maxtmp)
c
      ivl = maxtmp/nintpo
c
      do iconflo = 1, nconf, ivl
         iconfhi = min(iconflo+ivl-1,nconf)
         k = 0
         do iconf = iconflo, iconfhi
            do i = 1,nintpo
               k = k + 1
               tmp(k) = ieor(ioconf(i,iconf),joconf(i))
            enddo
         enddo
c
         call selci_vpop32(k, tmp, out)
c
         k = 0
         do iconf = iconflo, iconfhi
            ipop = 0
            do i = 1,nintpo
               k = k + 1
               ipop = ipop + out(k)
            enddo
            iexcit(iconf) = ipop
         enddo
      enddo
c
      end
