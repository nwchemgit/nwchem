C> \ingroup selci
C> @{
      subroutine selci_inicij(w1,w2)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
#include "ceij.fh"
c
      dimension w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1)
      data big /0.17976931348623132d+309/
c
c     make the coupling coefficient arrays <*|eij|*> for
c     up to 8 unpaired electrons (hard dimensions in ceij)
c
c     make sure of unambiguous disaster if unitialized
c     elements are accessed
c
      call dfill(nfc*nfc*nsc*(nsc+1),big,case1,1)
      call dfill(nfc*nfc2*nsc*nsc,big,case1,1)
      call dfill(nfc*nfc*nsc*(nsc+1),big,case4,1)
c
      do 5 ns = nsmax,0,-2
         if (ns.le.8) then
            nseij = ns
            ns4eij = min(ns,nsmax-2)
            goto 6
         endif
 5    continue
      write(6,*) ' nsmax, nseij, ns4eij ',nsmax, nseij, ns4eij
      call errquit('inicij: failed to initialize nseij etc. ',nseij,
     &       CALC_ERR)
 6    continue
      numf = nf(nseij)
      numf2 = nf(nseij-2)
      ns = nseij
      nsdiff = nsmax - ns
c
c     case 1: iocc(i)=1, iocc(j)=0
c
      do 10 iib = 1,ns
         do 20 jjb = 1,ns+1
            jb = jjb + nsdiff
            if (jjb.gt.iib) jb = jb - 1
            ib = iib + nsdiff
            call selci_axbt(w1(1,1,ib),nfmax, w1(1,1,jb),nfmax,
     $           case1(1,1,iib,jjb),nfc, numf, numf, numf)
 20      continue
 10   continue
c
c     case 2: iocc(i) = 1, iocc(j) = 1
c
c     case 3: iocc(i) = 3, iocc(j) = 0 is the same as swapping the
c             indices ib and jb and then transposing the matrix
c
      do 30 iib = 1,ns
         do 40 jjb = 1,ns
            if (iib.eq.jjb) goto 40
            jb = jjb + nsdiff
            if (jjb.gt.iib) jb = jb - 1
            ib = iib + nsdiff
            call selci_axbt(w1(1,1,ib),nfmax, w2(1,1,jb),nfmax2,
     $           case2(1,1,iib,jjb),nfc, numf, numf, numf2)
 40      continue
 30   continue
c
c     case 4: iocc(i) = 3, iocc(j) = 1
c
      ns = ns4eij
      numf = nf(ns4eij)
      numf2 = nf(ns+2)
      nsdiff = nsmax - (ns+2)
c
      do 50 iib = 1,ns+1
         do 60 jjb = 1,ns
            jb = jjb + nsdiff
            if (jjb.ge.iib) jb = jb + 1
            ib = iib + nsdiff
            call selci_axbt(w2(1,1,ib),nfmax2, w2(1,1,jb),nfmax2,
     $           case4(1,1,iib,jjb),nfc, numf, numf2, numf)
 60      continue
 50   continue
c
      end
C> @}
