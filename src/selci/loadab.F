C> \ingroup selci
C> @{
      subroutine selci_loadab(h, g, int12, int34)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      dimension h(nnorbs), g(numint), int12(nnorbs), int34(nnorbs)
      dimension gbuf(1024),ijkl(4,1024)
      integer buflen
#include "bitops.fh"
c
      call dfill(numint, 0.0d0, g, 1)
c
      rewind iflint
      read(iflint) ! titlE
      read(iflint) junk1, junk2, buflen
      read(iflint) repnuc
      read(iflint) ! map
      read(iflint) ! sym 
      read(iflint) ! eval
      read(iflint) ! occ
      read(iflint) h
c
      if (buflen .gt. 1024) call errquit
     $  ('selci_loadab: buflen > hard dimension', buflen, INPUT_ERR)
c
      write(6,19) repnuc
 19   format(' nuclear replusion energy is ',f20.12)
c
10    call moints_read_buf(iflint, n, gbuf, ijkl)
      do 20 iw = 1, n
         i = ijkl(1,iw)
         j = ijkl(2,iw)
         k = ijkl(3,iw)
         l = ijkl(4,iw)
         gg = gbuf(iw)
c
         if (abs(gg).lt.1.0d-10) goto 20
         if (i.lt.j) then
            ii = i
            i = j
            j = ii
         endif
         if (k.lt.l) then
            kk = k
            k = l
            l = kk
         endif
         ij = iky(i) + j
         kl = iky(k) + l
         if (ij.lt.kl) then
            ijij = ij
            ij = kl
            kl = ijij
         endif
         ijs = ieor(isym(i),isym(j))
         kls = ieor(isym(k),isym(l))
         if (ijs.ne.kls) then
            write(ifllog,*) ' gg,i,j,k,l,ijs,kls,',gg,i,j,k,l,ijs,kls
            call errquit('integral should be zero by symmetry',0,
     &       CALC_ERR)
         endif
         g(int12(ij)+int34(kl)) = gg
20    continue
      if (n .gt. 0) goto 10
c
      close(iflint,status='keep')
c
      end
C> @}
