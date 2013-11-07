C> \ingroup selci
C> @{
C>
      subroutine selci_makef(f,h,g,int12,int34,iocc,listd,lists,ns,nd,
     $     odonly)
*
* $Id$
*
#include "implicit.fh"
#include "ciinfo.fh"
      dimension f(*),h(*),g(*),int12(*),int34(*),listd(*),lists(*),
     $     iocc(255)
      logical odonly
c
      dimension ind(255),rocc(255),indir(255),indrr(255)
c     
c     f(ij) = hij + sum(r) (nr*(ij|rr)-(ir|rj))
c     where nr is the occupancy (1 or 2) of orbital r
c     and ij is a triangular index
c     
c     listd(*) is the list of doubly occupied orbitals (nd of them)
c     lists(*) is the list of singly occupied orbitals (ns of them)
c     
c     12/13/89 ... only make for at least one index occupied
c              ... lots of other optimisations
c
c     2/3/96   ... only compute diag. if that is all that is needed
c     
#include "stfunc.fh"
c     
c     call dfill(nnorbs,0.0d0,f,1)
c
      do 5 i = 1,ns
         ind(i) = lists(i)
         rocc(i) = 1.0d0
 5    continue
      do 6 i = 1,nd
         ind(i+ns) = listd(i)
         rocc(i+ns) = 2.0d0
 6    continue
      nocc = ns+nd
      do 7 i = 1,nocc
         indrr(i) = iky(ind(i)) + ind(i)
 7    continue
c
      do 10 i = 1,norbs
         if (iocc(i).gt.0) then
            do 15 kr = 1,nocc
               indir(kr) = itrian(i,ind(kr))
 15         continue
c
            if (odonly) then
               ii = itrian(i,i)
               s = h(ii)
               do kr = 1,nocc
                  ir = ind(kr)
                  iirr = intadr(ii,indrr(kr))
                  irri = intadr(indir(kr),itrian(ir,i))
                  s = s + rocc(kr)*g(iirr) - g(irri)
               enddo
               f(ii) = s
            else
               do 20 j = 1,norbs
                  if (isym(i).eq.isym(j)) then 
                     ij = itrian(i,j)
                     s = h(ij)
                     do 30 kr = 1,nocc
                        ir = ind(kr)
                        ijrr = intadr(ij,indrr(kr))
                        irrj = intadr(indir(kr),itrian(ir,j))
                        s = s + rocc(kr)*g(ijrr) - g(irrj)
 30                  continue
                     f(ij) = s
                  endif
 20            continue
            endif
c
         endif
 10   continue
c     
      end
C>
C> @}
