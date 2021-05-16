C> \ingroup selci
C> @{
      subroutine selci_gtijkl(ind,ioconf,joconf,nintpo,nbitpi,
     $     iocc,icase)
*
* $Id$
*
      dimension ind(*),ioconf(nintpo),joconf(nintpo),iocc(*)
      dimension ijocc(4)
#include "bitops.fh"
#include "errquit.fh"
c
c     (ik->jl)|I> = |J>. Determine ijkl from bit packed occupations
c     and return in ind(*) in order ikjl. 
c     icase returns a value which labels the possible cases
c
c     ii -> jj  icase = 1, ind(1)=i, ind(2)= j
c     ik -> jj  icase = 2, ind(1)=i, ind(2)=k, ind(3)=j
c     ii -> jl  icase = 3, ind(1)=i, ind(2)=j, ind(3)=l
c     ik -> jl  icase = 4, ind(1)=i, ind(2)=k, ind(3)=j, ind(4)=l
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
            goto 20
         endif
         ibase = ibase + nelpi
 10   continue
c
c     Now have indices and delta occupations ... which is which ?
c     First order by iocc(i)-jocc(i)
c
      do 30 i = 1,ngot-1
         do 40 j = i+1,ngot
            if (ijocc(i).lt.ijocc(j)) then
               ii = ind(i)
               ind(i) = ind(j)
               ind(j) = ii
               ii = ijocc(i)
               ijocc(i) = ijocc(j)
               ijocc(j) = ii
            endif
 40      continue
 30   continue
c
      if (ngot.eq.4) then
         icase = 4
      else if (ngot.eq.3) then
         if (ijocc(2).gt.0) then
            icase = 2
         else
            icase = 3
         endif
      else if (ngot.eq.2) then
         icase = 1
c        if (ijocc(2).gt.ijocc(1)) then
c          ii = ind(1)
c          ind(1) = ind(2)
c          ind(2) = ii
c        endif
      else
         call errquit('gtijkl: ngot is wierd',ngot, INPUT_ERR)
      endif
c
      end
C> @}
