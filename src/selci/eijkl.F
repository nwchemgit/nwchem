C> \ingroup selci
C> @{
C>
      subroutine selci_eijkl(e,i,j,k,l,ns,indbar,iocc,w1,w2,
     $     work,numf,numf2)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      dimension e(numf,*),indbar(*),iocc(*),w1(nfmax,nfmax,nsmax),
     $     w2(nfmax2,nfmax,nsmax-1),work(numf,*)
c
      parameter (idimt=167*167) ! 167 = max nf for ns=12
      dimension temp(idimt)
c
      dimension nsdel(4)
c     save statement to force temp into static memory
      save temp
      data nsdel/0,-2,+2,0/
c
c     Compute e(u,v) = <Iu|eijkl|Jv> where u,v label the
c     spin functions and the orbital configuration I is specified
c     by indbar and iocc. 
c
c     Works for i.ne.j.ne.k.ne.l
c
c     Assumes that arrives with occupation of i .ge. that of k
c     !!! I don't think that this version does assume this ...
c     !!! the previous one did
c
      if (iocc(i).lt.iocc(k)) call errquit('eijkl: iocc ',iocc(i),
     &       UNKNOWN_ERR)
c
c     iocc(i) = 0 (uoc), 1 (soc), 3 (doc)
c
c     possible occupations
c
c     i  k  (i+k)/2           j  l  j*2+l+1
c     ----  -------           ----  -------
c     1  1     1              0  0     1
c     3  1     2              0  1     2
c     3  3     3              1  0     3
c                             1  1     4
c
c     this version of eijkl always uses eijkl = eij*ekl
c
      call selci_eij(work,i,j,ns,indbar,iocc,w1,w2)
c
      kbar = indbar(k)
      lbar = indbar(l)
      kb = kbar
      lb = lbar
c
      if (iocc(i).eq.1) then
         if (k.gt.i) kb = kb - 1
         if (l.gt.i) lb = lb - 1
      else
         if (k.gt.i) kb = kb + 1
         if (l.gt.i) lb = lb + 1
      endif
      if (iocc(j).eq.1) then
         if (k.gt.j) kb = kb - 1
         if (l.gt.j) lb = lb - 1
      else
         if (k.gt.j) kb = kb + 1
         if (l.gt.j) lb = lb + 1
      endif
c
      nsint = ns + nsdel(iocc(i)+iocc(j))
      ns2 = nsint+ nsdel(iocc(k)+iocc(l))
      numfin = nf(nsint)
      numf2 = nf(ns2)
      if (idimt.lt.numfin*numf2)
     $     call errquit('new eijkl: get rid of stupid local array',
     $     numfin*numf2, UNKNOWN_ERR)
c
      indbar(k) = kb
      indbar(l) = lb
      call selci_eij(temp,k,l,nsint,indbar,iocc,w1,w2)
      indbar(k) = kbar
      indbar(l) = lbar
c
      call selci_axb(work,numf,temp,numfin,e,numf,numf,numfin,numf2)
c
      end
C>
C> @}
