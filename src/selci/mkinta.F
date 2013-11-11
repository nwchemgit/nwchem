C> \ingroup selci
C> @{
      subroutine selci_mkinta(int12, int34)
*
* $Id$
*
#include "implicit.fh"
#include "ciinfo.fh"
      dimension itemp(8), int12(*), int34(*)
#include "bitops.fh"
c
c     make the integral addressing arrays
c     int34(ij) = position of ij amoungst pairs of same symmetry
c     int12(ij) = sum of int12(kl) for kl<ij and same symmetry
c
      call ifill(nsym,0,itemp,1)
      call ifill(nnorbs,0,int12,1)
      call ifill(nnorbs,0,int34,1)
      ij = 0
      do 10 i = 1,norbs
         do 20 j = 1,i
            ijs = ieor(isym(i),isym(j)) + 1
            itemp(ijs) = itemp(ijs) + 1
            ij = ij + 1
            int34(ij) = itemp(ijs)
 20      continue
 10   continue
c
      ipt = 0
      kl = 0
      do 30 k = 1,norbs
         do 40 l = 1,k
            kl = kl + 1
            int12(kl) = ipt
            ipt = ipt + int34(kl)
 40      continue
 30   continue
c     
      numint = ipt + int34(kl)
c
      end
C> @}
