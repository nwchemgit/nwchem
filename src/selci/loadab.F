      subroutine selci_loadab(h, g, int12, int34)
#include "implicit.fh"
#include "ciinfo.fh"
      dimension h(nnorbs), g(numint), int12(nnorbs), int34(nnorbs)
      common / crap /gbuf(1024),ijkl(1024)
#include "bitops.fh"
c
c     load the symmetry packed integrals in core
c     this version is reading from the binary version of ascii list
c     Is terminated by the repulsion energy.
c
      call dfill(numint, 0.0d0, g, 1)
c
10    call dfill(1536, 0.0d0, gbuf, 1)
      call selci_sread(iflint, gbuf, 1536)
      do 20 iw = 1,1024
         i = ishft(ijkl(iw),-24)
         j = iand(ishft(ijkl(iw),-16),255)
         k = iand(ishft(ijkl(iw),-8 ),255)
         l = iand(ijkl(iw),255)
         gg = gbuf(iw)
      if (k.gt.0) then
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
            call errquit('integral should be zero by symmetry',0)
         endif
c         write(6,*) i,j,k,l,gg
         g(int12(ij)+int34(kl)) = gg
         goto 20
      else if (i.gt.0) then
         if (abs(gg).lt.1.0d-10) goto 20
         if (i.lt.j) then
            ii = i
            i = j
            j = ii
         endif
         if (i.gt.norbs) goto 20
         ij = iky(i) + j
         ijs = ieor(isym(i),isym(j))
         if (ijs.ne.0) then
            write(ifllog,*) ' h,i,j,ijs',gg,i,j,ijs
            call errquit(' integral should be zero by symmetry',0)
         endif
c         write(6,*) i,j,gg
         h(ij) = gg
         goto 20
      else if (i.eq.0) then
         repnuc = gg
         write(6,19) repnuc
 19      format(' nuclear replusion energy is ',f20.12)
         goto 30
      endif
20    continue
      goto 10
30    continue
c
      close(iflint,status='keep')
c
      end
