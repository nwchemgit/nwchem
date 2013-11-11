C> \ingroup selci
C> @{
      subroutine selci_loadg(g,int12,int34)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      dimension g(numint), int12(nnorbs), int34(nnorbs)
c
c     gamess has 340 integrals per block, columbus 1364
      parameter(n2eabf=1364)
      dimension q(n2eabf),iq(n2eabf),jq(n2eabf),kq(n2eabf),lq(n2eabf)
      dimension sum(5,3)
c      real *8 x, x2, w1, w2
#include "bitops.fh"
c
c     if gamess
c       read in the two electron integrals from ed6 at block 1
c       GAMESS INTERFACE IS NO LONGER FUNCTIONAL
c     if columbus
c       read from moints
c
      call dfill(numint, 0.0d0, g, 1)
      call dfill(15, 0.0d0, sum, 1)
c
      nread = 0
*      if (inttyp.eq.1) call search(1,7)
 10   if (inttyp.eq.1) then
         call errquit('loadg: gamess interface not functional',0,
     &       CAPMIS_ERR)
*         call gtintb(7,q,iq,jq,kq,lq,nw)
      else if (inttyp.eq.2) then
         call selci_colint(iflint,q,iq,jq,kq,lq,nw)
      else
         call errquit('loadg: strange inttyp ',inttyp, INPUT_ERR)
      endif
      do 20 iw = 1,nw
         i = iq(iw)
         j = jq(iw)
         k = kq(iw)
         l = lq(iw)
         gg = q(iw)
         if (abs(gg).lt.1.0d-9) goto 20
c
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
c     allow for reduced space calculations
         if (max(i,k).gt.norbs) goto 20
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
            call errquit('loadg: integral should be zero by symmetry',0,
     &       CALC_ERR)
         endif
         g(int12(ij)+int34(kl)) = gg
         nread = nread + 1
 20   continue
      if (nw.gt.0) goto 10
c
      write(6,*) ' number of non-zero 2e- integrals ',nread
c      write(nlist,'(1x,1p5e15.7)')sum
      call util_flush(6)
      end
C> @}
