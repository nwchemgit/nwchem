C> \file eij.F
C> Coupling coefficients eij
C>
C> \ingroup selci
C> @{
C>
C> \brief Coupling coefficients involving single excitations
C>    
C> Calculate coupling coefficients [1,2]
C> \f{eqnarray*}{
C>    e(u,v) &=& <Iu|E_{ij}|Jv>
C> \f}
C> where u, v label the spin functions
C> and the orbital occupation of I is specified by ns, iocc
C> and indbar. Only for i.ne.j.
C>    
C> \f{eqnarray*}{
C>    <Iu|Eij|Jv> &=& \sum_q <Iu|E_{ia}|Sq><Sq|E_{aj}|Jv>
C> \f}
C>    
C> possible cases are:
C>    
C> iocc(i) | iocc(j) | case
C> --------|---------|-----
C> 1       |  0      | 1
C> 1       |  1      | 2
C> 3       |  0      | 3
C> 3       |  1      | 4
C>
C> ### References ###
C>
C> [1] P.J. Knowles, H.-J. Werner, "An efficient method for the 
C>     evaluation of coupling coefficients in configuration interaction
C>     calculations", Chem. Phys. Lett. <b>145</b> (1988) 514-522, doi:
C>     <a href="https://doi.org/10.1016/0009-2614(88)87412-8">
C>     10.1016/0009-2614(88)87412-8</a>.
C>
C> [2] P.E.M. Siegbahn, "A new direct CI method for large CI expansions
C>     in a small orbital space", Chem. Phys. Lett. <b>109</b> (1984)
C>     417-423, doi: 
C>     <a href="https://doi.org/10.1016/0009-2614(84)80336-X">
C>     10.1016/0009-2614(84)80336-X</a>.
C>
      subroutine selci_eij(e,i,j,ns,indbar,iocc,w1,w2)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
#include "ceij.fh"
c     
      dimension e(*),indbar(*),iocc(*),
     $     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1)
c     
c     e(u,v) = <Iu|Eij|Jv>, where u,v label the spin functions
c     and the orbital occupation of I is specified by ns, iocc
c     and indbar. Only for i.ne.j.
c     
c     <Iu|Eij|Jv> = sum(q) <Iu|Eia|Sq><Sq|Eaj|Jv>
c     
c     possible cases are:
c     
c     iocc(i)  iocc(j)  case
c     1         0       1
c     1         1       2
c     3         0       3
c     3         1       4
c     
      if (i.eq.j) call errquit('eij, i=j',i, UNKNOWN_ERR)
c     
      numf = nf(ns)
      icase = iocc(i)+iocc(j)
      goto (10,20,30,40) icase
      call errquit('eij: fell thru goto',icase, UNKNOWN_ERR)
c     
c     ni=1, nj=0 : e(u,v) = sum(q) w1(u,q,ib)*w1(v,q,jb)
c     
 10   if (ns.le.nseij) then
         nsdiff = nseij - ns
         ib = indbar(i) + nsdiff
         jb = indbar(j) + nsdiff
         if (numf.eq.1) then
            e(1) = case1(1,1,ib,jb)
         else
            ie = 0
cvd$ shortloop
            do 11 is = 1,numf
cvd$ shortloop
               do 12 js = 1,numf
                  e(ie+js) = case1(js,is,ib,jb)
 12            continue
               ie = ie + numf
 11         continue
         endif
      else
         nsdiff = nsmax - ns
         ib = indbar(i) + nsdiff
         jb = indbar(j) + nsdiff
         if (j.gt.i) jb = jb - 1
c     call selci_mxma(w1(1,1,ib),1,nfmax,w1(1,1,jb),nfmax,1,e,1,numf,
c     $     numf,numf,numf)
         call selci_axbt(w1(1,1,ib),nfmax, w1(1,1,jb),nfmax, e,numf,
     $        numf, numf, numf)
      endif
      return
c     
c     ni=1, nj=1 : e(u,v) = sum(q) w1(u,q,ib)*w2(v,q,jb)
c     
 20   numf2 = nf(ns-2)
      if (ns.le.nseij) then
         nsdiff = nseij - ns
         ib = indbar(i) + nsdiff
         jb = indbar(j) + nsdiff
         ie = 0
cvd$  shortloop
         do 22 js = 1,numf2
cvd$  shortloop
            do 23 is = 1,numf
               e(ie+is) = case2(is,js,ib,jb)
 23         continue
            ie = ie + numf
 22      continue
      else
         nsdiff = nsmax - ns
         ib = indbar(i) + nsdiff
         jb = indbar(j) + nsdiff
         if (j.gt.i) jb = jb - 1
c     call selci_mxma(w1(1,1,ib),1,nfmax,w2(1,1,jb),nfmax2,1,e,1,numf,
c     $     numf,numf,numf2)
         call selci_axbt(w1(1,1,ib),nfmax, w2(1,1,jb),nfmax2, e,numf,
     $        numf, numf, numf2)
      endif
      return
c     
c     ni=2, nj=0 : e(u,v) = sum(q) w2(u,q,ib)*w1(v,q,jb)
c     
 30   numf2 = nf(ns+2)
      if (ns+2.le.nseij) then
         nsdiff = nseij - ns - 2
         ib = indbar(i) + nsdiff
         jb = indbar(j) + nsdiff
         if (i.gt.j) then
            ib = ib + 1
         else
            jb = jb + 1
         endif
         ie = 1
cvd$  shortloop
         do 32 is = 1,numf
            iee = ie
cvd$  shortloop
            do 33 js = 1,numf2
               e(iee) = case2(js,is,jb,ib)
               iee = iee + numf
 33         continue
            ie = ie + 1
 32      continue
      else
c     
         nsdiff = nsmax - ns - 2
         ib = indbar(i) + nsdiff
         jb = indbar(j) + nsdiff
         if (j.gt.i) jb = jb + 1
c     call selci_mxma(w2(1,1,ib),1,nfmax2,w1(1,1,jb),nfmax,1,e,1,numf,
c     $     numf,numf2,numf2)
         call selci_axbt(w2(1,1,ib),nfmax2, w1(1,1,jb),nfmax, e,numf,
     $        numf, numf2, numf2)
      endif
      return
c     
c     ni=2, nj=1 : e(u,v) = sum(q) w2(u,q,ib)*w2(v,q,jb)
c     
 40   continue
      if (ns.le.ns4eij) then
         nsdiff = ns4eij - ns 
         ib = indbar(i) + nsdiff
         jb = indbar(j) + nsdiff
         if (numf.eq.1) then
            e(1) = case4(1,1,ib,jb)
         else
            ie = 0
cvd$  shortloop
            do 42 is = 1,numf
cvd$  shortloop
               do 43 js =1,numf
                  e(ie+js) = case4(js,is,ib,jb)
 43            continue
               ie = ie + numf
 42         continue
         endif
      else
         numf2 = nf(ns+2)
         nsdiff = nsmax - ns - 2
         ib = indbar(i) + nsdiff
         jb = indbar(j) + nsdiff
         if (j.gt.i) jb = jb + 1
c     call selci_mxma(w2(1,1,ib),1,nfmax2,w2(1,1,jb),nfmax2,1,e,1,numf,
c     $     numf,numf2,numf)
         call selci_axbt(w2(1,1,ib),nfmax2, w2(1,1,jb),nfmax2, e,numf,
     $        numf, numf2, numf)
      endif
      return
c     
      end
C>
C> @}
