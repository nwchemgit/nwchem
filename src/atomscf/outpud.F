      subroutine outpud(copn,cc,ntest,iwr)
C$Id$
      implicit none
c
      real *8 copn(*), cc(*)
      integer ntest, iwr
c
#include "cguess.fh"
c
      integer i, j, k
      integer nbc1, naddr, noddr, jfirst, jlast
c
      if (ntest.gt.0) write (iwr,6010) nitscf , energ , cin , vir
      if (ntest.gt.0) write (iwr,6030)
      nbc1 = nbc(1)
      do i = 2 , nsym
         nbc1 = max(nbc1,nbc(i))
      enddo
      naddr = 0
      noddr = 0
      do i = 1 , nsym
         do j = 1 , ncsh(i) + nosh(i)
            do k = 1 , nbc(i)
               naddr = naddr + 1
               noddr = noddr + 1
               copn(noddr) = cc(naddr)
            enddo
            if (nbc(i).lt.nbc1) then
               do k = nbc(i) + 1 , nbc1
                  noddr = noddr + 1
                  copn(noddr) = 0.0d0
               enddo
            end if
         enddo
         naddr = naddr + nbc(i)*(nbc(i)-ncsh(i)-nosh(i))
      enddo
      jfirst = 1
 70   jlast = min(jfirst+7,nsht)
      if (ntest.gt.0) then
         write (iwr,6020) (eps(j),j=jfirst,jlast)
         write (iwr,6030)
         do i = 1 , nbc1
            write (iwr,6040) (copn(i+(j-1)*nbc1),j=jfirst,jlast)
         enddo
      end if
      if (jlast.eq.nsht) return
      jfirst = jfirst + 8
      if (ntest.gt.0) write (iwr,6030)
      go to 70
 6010 format (/,8x,'final scf results at iteration',i4,/,8x,
     +        'total hf energy',4x,'kinetic energy',4x,'virial theorem',
     +        /4x,3(e19.10),//,8x,'orbital energies and eigenvectors')
 6020 format (/,2x,8(2x,f12.5),1x)
 6030 format (' ')
 6040 format (2x,8f14.6)
      end
