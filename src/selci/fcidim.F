C> \ingroup selci
C> @{
      double precision function selci_fcidim(nelec, multi, issss,
     $     norbs, isym)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
c
      dimension isym(norbs)
c
      parameter (maxa=40, maxb=40)
      dimension r1(0:7, 0:maxb, 0:maxa), r2(0:7, 0:maxb, 0:maxa)
#include "bitops.fh"
c
c     compute the spin adapted fullci dimension for nelec 
c     electrons in norbs orbitals, state symmetry issss (0:7) with
c     symmetry of each orbital as isym(*) (0:7).
c
c     construct the drt one level at a time using floating point
c     numbers to accumulate the weights.
c
c     define head of graph
c
*      write(6,*) ' fcidim: nelec, multi, issss, norbs',
*     & nelec, multi, issss, norbs
*      write(6,*) ' fcidim: isym', isym
      call dfill(8*(maxa+1)*(maxb+1), 0.0d0, r1, 1)
      ib = multi-1
      ia = (nelec-ib)/2
      if (ia.gt.maxa .or. ib.gt.maxb .or.
     $     issss.lt.0 .or. issss.gt.7) then
         write(6,*) ' nelec, multi, issss ', nelec, multi, issss
         write(6,*) ' ib, ia , maxb, maxa ',ib,ia,maxb, maxa
         call errquit(' fcidim: exceeded maxa or maxb or issss ',0,
     &       INPUT_ERR)
      endif
      r1(issss,ib,ia) = 1.0d0
c
c     now construct the drt one level at a time
c
      do 10 iorb = norbs,1,-1
*         do ia = 0,maxa
*            do ib = 0,maxb
*               if (r1(0,ib,ia).ne.0.0d0) write(6,*) iorb, ia, ib, 
*     $              r1(0,ib,ia)
*            enddo
*         enddo
         call dfill(8*(maxa+1)*(maxb+1),0.0d0,r2,1)
c
c     loop thru rows at this level
c
         do 20 ia = 0,maxa
            do 30 ib = 0,maxb
               do 40 is = 0,7
c
c     if this row has non-zero weight then propagate down
c
                  if (r1(is,ib,ia).gt.0.0d0) then
c
c     arc=0 : a->a, b->b
c
                     ian = ia
                     ibn = ib
                     isn = is
                     r2(isn,ibn,ian) = r2(isn,ibn,ian) + r1(is,ib,ia)
c
c     arc=1 : a->a, b->b-1, symmetry changes
c
                     ian = ia
                     ibn = ib-1
                     isn = ieor(isym(iorb),is)
                     if (ibn.ge.0) then
                        r2(isn,ibn,ian) = r2(isn,ibn,ian) + r1(is,ib,ia)
                     endif
c
c     arc=2 : a->a-1, b->b-1, symmetry changes
c
                     ian = ia-1
                     ibn = ib+1
                     isn = ieor(isym(iorb),is)
                     if (ian.ge.0 .and. ibn.le.maxb) then
                        r2(isn,ibn,ian) = r2(isn,ibn,ian) + r1(is,ib,ia)
                     else if (ibn.gt.maxb) then
                        call errquit('fcidim: ibn.gt.maxb',ibn,
     &       INPUT_ERR)
                     endif
c
c     arc=3 : a->a-1, b->b
c
                     ian = ia-1
                     ibn = ib
                     isn = is
                     if (ian.ge.0) then
                        r2(isn,ibn,ian) = r2(isn,ibn,ian) + r1(is,ib,ia)
                     endif
                  endif
 40            continue
 30         continue
 20      continue
c
         call dcopy(8*(maxa+1)*(maxb+1),r2,1,r1,1)
 10   continue
c
      selci_fcidim = r1(0,0,0)
c
      end
C> @}
