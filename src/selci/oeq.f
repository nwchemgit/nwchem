C> \ingroup selci
C> @{
      logical function selci_oeq(noconf,ioconf,joconf,nintpo)
*
* $Id$
*
      integer ioconf(nintpo,noconf),joconf(nintpo)
c     
c     ivl set so can use short vector lengths
c     
      parameter (ivl=32)
      logical owrk(ivl)
c     
c     return true if orbital configuration joconf is same as
c     any of the configurations in ioconf. false otherwise.
c     
c     strip mine into blocks to minimize redundant work but
c     keep vectorization
c
      nleft = noconf
      ioff = 0
c     
 100  ndo = min(ivl,nleft)
c     iwrk3 will have the difference occupations accumulated into it
cvd$  shortloop
      do 5 i = 1,ndo
         owrk(i) = .true.
 5    continue
      do 10 iw = 1,nintpo
         itest = joconf(iw)
cvd$  shortloop
         do 20 i = 1,ndo
            owrk(i) = (itest .eq. ioconf(iw,i+ioff)) .and. owrk(i)
 20      continue
 10   continue
cvd$  shortloop
      do 40 i = 1,ndo
         if (owrk(i)) then
            selci_oeq = .true.
            return
         endif
 40   continue
c     
      ioff = ioff + ndo
      nleft = nleft - ndo
      if (nleft.gt.0) goto 100
c     
      selci_oeq = .false.
c     
      end
C> @}
