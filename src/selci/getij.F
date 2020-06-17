C> \ingroup selci
C> @{
      subroutine selci_getij(i,j,ioconf,joconf,nintpo,nbitpi,iocc)
      dimension ioconf(nintpo),joconf(nintpo),ind(2),ijocc(2),
     $     iocc(*)
#include "bitops.fh"
#include "errquit.fh"
c
c     (i->j)|I> = |J>. Determine i and j from the
c     bit packed occupations
c
      nelpi = nbitpi/2
      nbit2 = 2 - nbitpi
      ngot = 0
      ibase = 0
      do 10 iword = 1,nintpo
         itemp = ieor(ioconf(iword),joconf(iword))
         jbase = 0
 20      lz = leadz(itemp)
         lz = lz - mod(lz,2)
         if (lz.ne.nbitpi) then
            jbase = jbase  + lz/2 + 1
            ngot = ngot + 1
            itemp = ishft(itemp,lz)
            k = ibase + jbase
            ind(ngot) = k
            ijocc(ngot) = iocc(k) - ieor(ishft(itemp,nbit2),iocc(k))
            itemp =ishft(itemp,2)
            if (ngot.eq.2) goto 30
            goto 20
         endif
         ibase = ibase + nelpi
 10   continue
      call errquit('getij: fell out of do 10',ngot, UNKNOWN_ERR)
c
 30   continue
c
c     now have indices and delta occupations ... which is which ?
c
      if (ijocc(1).gt.ijocc(2)) then
         i = ind(1)
         j = ind(2)
      else
         i = ind(2)
         j = ind(1)
      endif
c
      end
C> @}
