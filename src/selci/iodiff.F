      integer function selci_iodiff(ioconf,joconf,nintpo)
*
* $Id: iodiff.F,v 1.3 1997-10-31 23:42:06 d3e129 Exp $
*
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
      integer ioconf(nintpo*nconf), iexcit(nconf)
      integer joconf(*)
#include "bitops.fh"
c
      parameter (maxtmp=4095)
      integer tmp(maxtmp), out(maxtmp)
c
      ivl = maxtmp/nintpo
c
      do iconflo = 1, nconf, ivl
         iconfhi = min(iconflo+ivl-1,nconf)
         icbase  = (iconflo-1)*nintpo
         nk = (iconfhi-iconflo+1)*nintpo
c
         if (nintpo .eq. 1) then
            itest1 = joconf(1)
            do k = 1, nk
               tmp(k  ) = ieor(ioconf(icbase+k  ),itest1)
            enddo
         else if (nintpo .eq. 2) then
            itest1 = joconf(1)
            itest2 = joconf(2)
            do k = 1, nk, 2
               tmp(k  ) = ieor(ioconf(icbase+k  ),itest1)
               tmp(k+1) = ieor(ioconf(icbase+k+1),itest2)
            enddo
         else if (nintpo .eq. 3) then
            itest1 = joconf(1)
            itest2 = joconf(2)
            itest3 = joconf(3)
            do k = 1, nk, 3
               tmp(k  ) = ieor(ioconf(icbase+k  ),itest1)
               tmp(k+1) = ieor(ioconf(icbase+k+1),itest2)
               tmp(k+2) = ieor(ioconf(icbase+k+2),itest3)
            enddo
         else
            k = 0
            do iconf = iconflo, iconfhi
               do i = 1,nintpo
                  tmp(k+i) = ieor(ioconf(icbase+i),joconf(i))
               enddo
               k = k + nintpo
               icbase = icbase + nintpo
            enddo
         endif
c
         call selci_vpop32(nk, tmp, out)
c
         k = 1
         if (nintpo .eq. 1) then
            do iconf = iconflo, iconfhi
               iexcit(iconf) = out(iconf-iconflo+1)
            enddo
         else if (nintpo .eq. 2) then
            do iconf = iconflo, iconfhi
               iexcit(iconf) = out(k) + out(k+1)
               k = k + 2
            enddo
         else if (nintpo .eq. 3) then
            do iconf = iconflo, iconfhi
               iexcit(iconf) = out(k) + out(k+1) + out(k+2)
               k = k + 3
            enddo
         else
            do iconf = iconflo, iconfhi
               ipop = 0
               do i = 1,nintpo
                  ipop = ipop + out(k)
                  k = k + 1
               enddo
               iexcit(iconf) = ipop
            enddo
         endif
      enddo
c
      end
