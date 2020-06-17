C> \file ci.F
C> Initialize CI module
C>
C> \ingroup selci
C> @{
C>
C> \brief CI initialization
C>
      subroutine selci_ci(q, lword)
*
* $Id$
*
#include "implicit.fh"
#include "ciinfo.fh"
#include "mptr.fh"
#include "cselcifiles.fh"
#include "global.fh"
      character *80 title
c
      integer lword !< [Input] The available amount of memory
      double precision q(lword) !< [Scratch] The main memory array
c
c     data initialization
c
      junk = selci_minit(lword)
      call selci_gettim(startc, startw)
c
      lencib = mitob(lencii) + mdtob(lencir)
      me = ga_nodeid()
      if (me.eq.0) then
        call selci_initci(title, 'ci')
        write (iflham) nci
      else
        iflham = 4
        open(iflham,form='unformatted',status='unknown',file=hamil)
      endif
      call ga_brdcst(1,repnuc,lencib,0)
c
c     make integral addressing arrays, and then load the integrals
c
      int12 = selci_mptr(mitod(nnorbs))
      int34 = selci_mptr(mitod(nnorbs))
      call selci_mkinta(q(int12), q(int34))
      ih = selci_mptr(nnorbs)
      ig = selci_mptr(numint)
      if (me.eq.0) then
        write(ifllog,7) numint
 7      format(' no. of integrals       ',i9)
        call selci_load(q(ih), q(ig), q(int12), q(int34))
        call util_flush(ifllog)
      endif
      call ga_brdcst(2,q(ih),mdtob(nnorbs),0)
      call ga_brdcst(3,q(ig),mdtob(numint),0)
c         some integral types don't make repnuc available until load
      call ga_brdcst(4,repnuc,mdtob(1),0)
c
c     get in the coupling coefficient arrays
c
      iw1 = selci_mptr(nfmax*nfmax*nsmax)
      iw2 = selci_mptr(nfmax2*nfmax*(nsmax-1))
      if (me.eq.0) then
         call selci_loadw(iflwmt,q(iw1),q(iw2),
     $        nfmax,nfmax2,nsmax)
         close(iflwmt,status='delete')
      endif
      call ga_brdcst(5,q(iw1),mdtob(nfmax*nfmax*nsmax),0)
      call ga_brdcst(6,q(iw2),mdtob(nfmax2*nfmax*(nsmax-1)),0)
      call selci_inicij(q(iw1),q(iw2))
c
c     allocate space for orbital configurations and
c     the index vector, then read them off the ciconf file
c
      ioconf = selci_mptr(mitod(nintpo*noconf))
      indxci = selci_mptr(mitod(noconf))
      if (me.eq.0) then
         call selci_rdconf(iflcon,q(ioconf),q(indxci),
     $        nintpo,noconf)
         close(iflcon,status='keep')
      endif
      call ga_brdcst(7,q(ioconf),mitob(nintpo*noconf),0)
      call ga_brdcst(8,q(indxci),mitob(noconf),0)
c      call selci_pconls(q(ioconf))
c
c     make the hamiltonian
c
      ihd = selci_mptr(nci)
      iwork1 = selci_mptr(max(nfmax*nfmax, nci))
      iwork2 = selci_mptr(nfmax*nfmax)
      iwork3 = selci_mptr(nfmax*nfmax)
      if = selci_mptr(nnorbs)
c
      call selci_makeh(q(ih), q(ig), q(int12), q(int34), q(iw1), q(iw2),
     $     q(ioconf), q(indxci), q(ihd), q(iwork1), q(iwork2),
     $     q(iwork3), q(if))
c
      if (ga_nodeid() .eq. 0) then
         junk = selci_mstats(0)
         call selci_gettim(endc, endw)
         write(6,1) 'CI', endc-startc, endw-startw
 1       format(/1x,a,': Total cpu',f10.1,' : Total wall',f10.1)
         call util_flush(6)
      endif
c
      call ga_sync()
c
      end
C>
C> @}
