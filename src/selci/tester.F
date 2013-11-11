C> \ingroup selci
C> @{
      subroutine selci_tester(q,h, g, int12, int34, w1, w2, ioconf,
     $     indxci,
     $     roots, ci, nconmx, thresh, ept, enew, iconf, ncold,
     $     vc, iocc, itemp, irange, ptnorm,
     $     ptnorm_mp,ept_mp,roots_mp)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
#include "mptr.fh"
#include "global.fh"
#include "stdio.fh"
c
c     argument declarations
c
      dimension q(*)
      dimension h(nnorbs),g(numint),int12(nnorbs),int34(nnorbs),
     $     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1),
     $     ioconf(nintpo,nconmx),indxci(nconmx),roots(nroot),
     $     ci(nci,nroot),ept(nroot),enew(nroot),vc(nfmax,nroot),
     $     iocc(255), itemp(nintpo), irange(21), ptnorm(nroot)
      dimension ptnorm_mp(nroot),ept_mp(nroot),roots_mp(nroot)
c
      logical selci_ointer, selci_oeq
c
c     Do the actual PT selection of new configurations
c
c     new test configuration ... have we seen it before, or
c     is it one of the reference set already ? Does it
c     have the correct multiplicity?
c
      ns = 0
      do 5 i = 1,norbs
         if (iocc(i).eq.1) ns = ns + 1
 5    continue
      if (mod(ns,2).ne.mod(multi-1,2) .or. ns.lt.multi-1) return
c
      if(selci_ointer(iconf-1,ioconf,itemp,nintpo)) return
c
      if(selci_oeq(ncold-iconf+1,ioconf(1,iconf),itemp,nintpo)) return
c
c     got a brand new configuration ... make vc = <test|H|ci roots>
c     and <test|H|test> (returned in work1). numf is returned.
c     (note now only get the diag of <test|H|test>).
c
      iwork1 = selci_mptr(nfmax*nfmax)
      ide = selci_mptr(nroot)
      idnorm = selci_mptr(nroot)
      iexv = selci_mptr(mitod(ncold))
      call selci_mkvc(q, h,g,int12,int34,w1,w2,ioconf,indxci,q(iexv),
     $     iconf,ncold,ci,vc,iocc,itemp,numf,q(iwork1))
c
c     Update energy and ptnorm ... the last correction computed
c     is the one used for the selection process
c
      if(ga_nodeid().eq.0)call util_flush(luout)
      if (iwpt.eq.0) then       ! EN only
         call selci_mkdele(numf, nroot, q(ide), vc, roots, q(iwork1),
     $        nfmax, q(idnorm))
         call selci_vadd(nroot, q(ide), 1, ept, 1, ept, 1)
         call selci_vadd(nroot, q(idnorm), 1, ptnorm, 1, ptnorm, 1)
      else if (iwpt.eq.1) then  ! MP only
         call selci_mpdele(numf, nroot, q(ide), vc, roots_mp,nfmax,
     &      q(idnorm),iocc)
         call selci_vadd(nroot, q(ide), 1, ept_mp, 1, ept_mp, 1)
         call selci_vadd(nroot, q(idnorm), 1, ptnorm_mp, 1, ptnorm_mp,1)
      else if (iwpt.eq.2) then  ! Print MP, select on EN
         call selci_mpdele(numf, nroot, q(ide), vc, roots_mp,nfmax,
     &      q(idnorm),iocc)
         call selci_vadd(nroot, q(ide), 1, ept_mp, 1, ept_mp, 1)
         call selci_vadd(nroot, q(idnorm), 1, ptnorm_mp, 1, ptnorm_mp,1)
         call selci_mkdele(numf, nroot, q(ide), vc, roots, q(iwork1),
     $        nfmax, q(idnorm))
         call selci_vadd(nroot, q(ide), 1, ept, 1, ept, 1)
         call selci_vadd(nroot, q(idnorm), 1, ptnorm, 1, ptnorm, 1)
      else if (iwpt.eq.3) then   ! print EN, select on MP
         call selci_mkdele(numf, nroot, q(ide), vc, roots, q(iwork1),
     $        nfmax, q(idnorm))
         call selci_vadd(nroot, q(ide), 1, ept, 1, ept, 1)
         call selci_vadd(nroot, q(idnorm), 1, ptnorm, 1, ptnorm, 1)
         call selci_mpdele(numf, nroot, q(ide), vc, roots_mp,nfmax,
     &      q(idnorm),iocc)
         call selci_vadd(nroot, q(ide), 1, ept_mp, 1, ept_mp, 1)
         call selci_vadd(nroot, q(idnorm), 1, ptnorm_mp, 1, ptnorm_mp,1)
      else
         call errquit( ' tester: iwpt = ',iwpt, UNKNOWN_ERR)
      endif
c
      test = abs(q(ide-1+idamax(nroot, q(ide), 1)))
      if (test.gt.0.0d0) then
         ihash = 1 - nint(2.0d0*log10(test))
         ihash = min(ihash,21)
         ihash = max(ihash,1)
         irange(ihash) = irange(ihash) + numf
      endif
c
      if (test.ge.thresh) then
c
c     A reference configuration with too many open shells makes
c     a calculation with the current code too expensive due to 
c     memory and cpu consumption by coupling coefficients.
c     Limit the max number of open shells in a reference to 8.
c     
         if (ns .gt. 8) then
            write(luout,101) ns, test 
 101        format(' rejecting new reference because 8 < ns =',i3,
     $           ' : tester = ', 1pd9.2)
         else
c
c     configuration interacts more than threshold ... add to
c     list and update enew
c
            call selci_vadd(nroot, q(ide), 1, enew, 1, enew, 1)
            noconf = noconf + 1
            if (noconf.gt.nconmx)
     $           call errquit(' select: nconf > nconmx',nconmx,
     &       INPUT_ERR)
            call selci_icopy(nintpo,itemp,1,ioconf(1,noconf),1)
 21         format(4x,5f14.7/)
         endif
      endif
c
      junk = selci_mfree(iwork1)
c
      end
C> @}
