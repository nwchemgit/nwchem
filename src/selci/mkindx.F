      subroutine selci_mkindx(ioconf,indxci,ncinew,nnsmax)
*
* $Id: mkindx.F,v 1.2 1997-10-31 23:42:13 d3e129 Exp $
*
#include "implicit.fh"
#include "ciinfo.fh"
c
      dimension ioconf(nintpo,noconf), indxci(noconf+1)
c
      dimension iocc(255)
#include "bitops.fh"
c
      ncinew = 0
      nnsmax = 0
      do 10 iconf = 1,noconf
         call selci_upkcon(norbs, iocc, ioconf(1,iconf), nintpo, nbitpi)
         is = 0
         ns = 0
c$doit asis
         do 20 j = 1,norbs
            if (iocc(j).eq.1) then
               is = ieor(is,isym(j))
               ns = ns + 1
            endif
 20      continue
         if (is.ne.issss) then
            call selci_wrtcon(ifllog, iocc, norbs)
            call errquit(' mkindx: bad symmetry ',is)
         endif
         if (mod(ns,2).ne.mod(multi-1,2) .or. ns.lt.multi-1) then
            call selci_wrtcon(ifllog, iocc, norbs)
            call errquit(' mkindx: a reference has invalid ns',i)
         endif
         nnsmax = max(ns,nnsmax)
         if (nnsmax.gt.nsmax) then
            call selci_wrtcon(ifllog, iocc, norbs)
            call errquit('mkindx: wmatrx nsmax small, need:',nnsmax)
         endif
         indxci(iconf) = ncinew
         ncinew = ncinew + nf(ns)
 10   continue
      indxci(noconf+1) = ncinew
c
      end
