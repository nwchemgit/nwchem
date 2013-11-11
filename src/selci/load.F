C> \ingroup selci
C> @{
      subroutine selci_load(h, g, int12, int34)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      dimension h(nnorbs), g(numint), int12(nnorbs), int34(nnorbs)
#include "bitops.fh"
c
c     this stuff for reading integrals
c
c
c     load the symmetry packed integrals in core
c     later will be able to include re-ordering vector here.
c
      call dfill(nnorbs, 0.0d0, h, 1)
      call dfill(numint, 0.0d0, g, 1)
      if (inttyp.eq.0) then
         call selci_loada(h, g, int12, int34)
      else if (inttyp.eq.1) then
         call errquit('load: gamess interface not functional',0,
     &       CAPMIS_ERR)
      else if (inttyp.eq.2) then
         call errquit('load: columbus interface not functional',0,
     &       CAPMIS_ERR)
      else if (inttyp.eq.3) then
         call selci_loadab(h, g, int12, int34)
      else
         call errquit('load: invalid inttyp',inttyp, INPUT_ERR)
      endif
c
      ij = 0
      do 10 i = 1,norbs
         do 20 j = 1,i
            ij = ij + 1
            if (abs(h(ij)).lt.1.0d-10) then
               h(ij) = 0.0d0
            else if (isym(i).ne.isym(j)) then
               write(6,*) i,j,ij,isym(i),isym(j),h(ij)
               call errquit('load: integral zero by symmetry',ij,
     &       CALC_ERR)
            endif
 20      continue
 10   continue
c
c      ngot = 0
c      do 30 i= 1,norbs
c         do 40 j = 1,i
c            ijsym = ieor(isym(i),isym(j))
c            ijoff = int12(iky(i)+j)
c            do 50 k = 1,i
c               ijksym = ieor(ijsym,isym(k))
c               lhi = k
c               if (k.eq.i) lhi = j
c               do 60 l = 1,lhi
c                  if (isym(l).ne.ijksym) goto 60
c                  ijkl = ijoff +int34(iky(k)+l)
c                  if (abs(g(ijkl)).gt.1.0d-10) ngot = ngot + 1
c 60            continue
c 50         continue
c 40      continue
c 30   continue
c
c      write(6,*) ' ngot in load ',ngot
      end
C> @}
