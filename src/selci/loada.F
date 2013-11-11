C> \ingroup selci
C> @{
      subroutine selci_loada(h, g, int12, int34)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      dimension h(nnorbs), g(numint), int12(nnorbs), int34(nnorbs)
#include "bitops.fh"
c
c     load the symmetry packed integrals in core
c     this version is reading from the ascii list of integrals
c
      call dfill(numint, 0.0d0, g, 1)
c
 10   read (iflint, *) gg, i, j, k, l
      if (k.gt.0) then
         if (abs(gg).lt.1.0d-10) goto 10
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
         if (max(i,k).gt.norbs) goto 10
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
         goto 10
      else if (i.gt.0) then
         if (abs(gg).lt.1.0d-10) goto 10
         if (i.lt.j) then
            ii = i
            i = j
            j = ii
         endif
         if (i.gt.norbs) goto 10
         ij = iky(i) + j
         ijs = ieor(isym(i),isym(j))
         if (ijs.ne.0) then
            write(ifllog,*) ' h,i,j,ijs',gg,i,j,ijs
            call errquit(' integral should be zero by symmetry',0,
     &       CALC_ERR)
         endif
         h(ij) = gg
         goto 10
      else if (i.eq.0) then
         repnuc = gg
         write(6,19) repnuc
 19      format(' nuclear replusion energy is ',f20.12)
         goto 10
      endif
c
      close(iflint,status='keep')
c
      end
C> @}
