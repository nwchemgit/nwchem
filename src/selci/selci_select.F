C> \ingroup selci
C> @{
C>
      subroutine selci_select(q, lword, nroot_arg, iwpt_arg, 
     $     thresh_arg, oupdate, energies)
*
* $Id$
*
#include "implicit.fh"
#include "ciinfo.fh"
#include "mptr.fh"
#include "cselcifiles.fh"
#include "global.fh"
#include "stdio.fh"
      character*80 title
      integer nconmx
c
      dimension q(lword)
      dimension irange(21), map(255)
      dimension energies(*)     ! Returns ci+pt energies for each root
c
c     selected ci+pt program ... parallel version 1.0 using message pass
c     ... each process needs to have own copy of all data including
c     all the integrals ... working on how to fix this.
c
c     data initialization
c
      nroot = nroot_arg         ! These used to be read in by process zero
      iwpt  = iwpt_arg
      thresh= thresh_arg
c
      call selci_gettim(startc, startw)
      lencib = mitob(lencii) + mdtob(lencir)
      if (ga_nodeid().eq.0) then
c
         call selci_initci(title, 'select')
c
         nconmx = 50000
         if (nconmx .lt. 4*noconf) nconmx = 4*noconf
c
         min2c = 5              ! not currently used
      endif
      call ga_brdcst(81,nconmx,mitob(1),0)
      call ga_brdcst(9,iwpt,mitob(1),0)
      call ga_brdcst(10,repnuc,lencib,0)
      call ga_brdcst(11,thresh,mdtob(1),0)
      call ga_brdcst(111,min2c,mitob(1),0)
      if (ga_nodeid().eq.0) then
c
c     iwpt=0 -> EN only (default)
c     iwpt=1 -> MP only
c     iwpt=2 -> EN select, print MP
c     iwpt=3 -> MP select, print EN
c
         if (iwpt.eq.0 .or. iwpt.eq.2) 
     $        write(ifllog,'(/1x,a,f15.10)')
     &        ' EN| Selection threshold ',thresh
         if (iwpt.eq.1 .or. iwpt.eq.3) 
     $        write(ifllog,'(/1x,a,f15.10)')
     &        ' MP| Selection threshold ',thresh
      endif
c
c     make integral addressing arrays, and then load the integrals
c     also load in the coupling coeffs and make the eia*eaj table
c     in the parallel version only process 0 does any file access
c     ... the other processes only communicate via message passing.
c
      int12 = selci_mptr(mitod(nnorbs))
      int34 = selci_mptr(mitod(nnorbs))
      call selci_mkinta(q(int12), q(int34))
      ig = selci_mptr(numint)
      ih = selci_mptr(nnorbs)
      iw1 = selci_mptr(nfmax*nfmax*nsmax)
      iw2 = selci_mptr(nfmax2*nfmax*(nsmax-1))
      ioconf = selci_mptr(mitod(nintpo*nconmx))
      indxci = selci_mptr(mitod(nconmx+1))
      iroots = selci_mptr(nroot)
      iroots_mp = selci_mptr(nroot)
      iept = selci_mptr(nroot)
      iept_mp = selci_mptr(nroot)
      ienew = selci_mptr(nroot)
      iptnrm = selci_mptr(nroot)
      iptnrm_mp = selci_mptr(nroot)
      ici = selci_mptr(nci*nroot)
      if (ga_nodeid() .eq. 0) then
         write(ifllog,7) numint
 7       format(' no. of integrals       ',i9)
         call selci_load(q(ih), q(ig), q(int12), q(int34))
         call selci_loadw(iflwmt, q(iw1), q(iw2), nfmax, nfmax2, nsmax)
         close(iflwmt,status='delete')
         call selci_rdconf(iflcon,q(ioconf),q(indxci),nintpo,noconf)
c
         call selci_orbanal(noconf, q(ioconf), nintpo, nbitpi, norbs,
     $        map)
c
         call selci_sread(iflciv,q(iroots),nroot)
         do 10 iroot = 0,nroot-1
            call selci_sread(iflciv,q(ici+iroot*nci),nci)
 10      continue
c
         if (iwpt.ge.1) then
c     
c......read in eigenvalues
c     
            call selci_geteig(eigins)
c     
c......calculate E0 for each root 
c     
            call selci_mke0mp(q(ioconf),q(indxci),q(ici),q(iroots_mp),
     $           noconf)
         endif
c     
      endif
      call ga_brdcst(12,q(ih),mdtob(nnorbs),0)
      call ga_brdcst(13,q(ig),mdtob(numint),0)
      call ga_brdcst(14,q(iw1),mdtob(nfmax*nfmax*nsmax),0)
      call ga_brdcst(15,q(iw2),mdtob(nfmax2*nfmax*(nsmax-1)),0)
      call ga_brdcst(16,q(ioconf),mitob(nintpo*nconmx),0)
      call ga_brdcst(17,q(indxci),mitob(nconmx),0)
      call ga_brdcst(18,q(iroots),mdtob(nroot),0)
      call ga_brdcst(19,q(ici),mdtob(nci*nroot),0)
      call ga_brdcst(191,map,mitob(norbs),0) ! sorted orbitals by occ
      if (iwpt.ge.1) then
         call selci_brdeig()
         call ga_brdcst(1901,q(iroots_mp),mdtob(nroot),0)
      endif
c         some integral types don't make repnuc available until load
      call ga_brdcst(191,repnuc,mdtob(1),0)
c
      call selci_inicij(q(iw1),q(iw2))
c
*      icon4 = selci_mptr(mitod(noconf*4))
*      call makeicon4(noconf,q(ioconf),nintpo,nbitpi,norbs,
*     $     map,q(icon4))
c
      call selci_gettim(tc,tw)
      if (ga_nodeid().eq.0) write(ifllog,11) tc-startc, tw-startw
11    format(/' time in initialisation: cpu=',f6.2,', wall=',f6.2/)
c
c     do the selection
c
      ncold = noconf
c     write(6,*) nodeid(), repnuc, multi, nelec, issss, norbs,
c    & nnorbs, nroot, nci, noconf, nintpo, nbitpi, numint, inttyp, 
c    & nsmax, nfmax, nfmax2, nsym, nbpsy, ih, iq, int12, int34,
c    & iw1, iw2, ioconf, indxci, iroots, ici, nconmx, ncold, thresh,
c    & iept, ienew, iptnrm
c     call util_flush(6)
c...... do selection
      call selci_slect(q,q(ih),q(ig),q(int12),q(int34),q(iw1),q(iw2),
     $     q(ioconf),q(indxci),q(iroots),q(ici),
     $     nconmx,ncold,thresh,q(iept),q(ienew),irange,q(iptnrm),
     &     q(iptnrm_mp),q(iept_mp),q(iroots_mp),min2c)
      call selci_gettim(fc,fw)
      if (ga_nodeid().eq.0) write(ifllog,21) fc-tc, fw-tw
21    format(/' time in select: cpu=',f9.2,', wall=',f9.2)
c
      call selci_ptidy(q(ioconf), nintpo, noconf, ncold, nconmx, nroot,
     &     q(iept), q(ienew), irange, q(iptnrm),
     &     q(iptnrm_mp),q(iept_mp))
c
      do i = 1, nroot
         energies(i) = q(iroots+i-1) + q(iept+i-1)
      enddo
c
      if (ga_nodeid().gt.0) goto 1010
c
c     make the new index for the ci
c
      call selci_mkindx(q(ioconf),q(indxci),ncinew,nnsmax)
c
c     print out analysis of interactions
c
      call selci_anal(irange, q(iroots), q(iept), q(iptnrm),
     $     q(ienew), ncinew, q(iptnrm_mp), q(iept_mp))
c
c     if we haven't found any more configuations signal
c     convergence by error termination ... this and the
c     .nonew file technique is not satisfactory ... need to
c     have a separate history file
c
      if (ncold.eq.noconf) then
         write(luout,81)
 81      format(//
     $        ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/
     $        ' ! select: no new configurations found !'/
     $        ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/)
         goto 1000
      endif
c
      if (oupdate) then
c     
c     write out new index file information
c     
         rewind iflcon
         call selci_mknsneed(q(ioconf),noconf,nintpo,nbitpi,norbs,
     $        nsneed)
         call selci_wthcon(iflcon, title, multi, nelec, issss, norbs,
     $        nnsmax, ncinew, noconf, nintpo, nbitpi, nbpsy, isym, nsym,
     $        inttyp,nsneed)
c     
         call selci_wtconf(iflcon,q(ioconf),nintpo*noconf,q(indxci),
     $        noconf)
c
c     remap the ci vectors to help diagonalisation ... note that
c     ci vectors were the last allocated so can realloc easily
c     
         junk = selci_mfree(ici)
         ici = selci_mptr(ncinew*nroot)
         call dfill(ncinew*nroot,0.0d0,q(ici),1)
c     
         rewind iflciv
         read (iflciv)
         call selci_sread(iflciv,q(iroots),nroot)
         do 20 iroot = 0,nroot-1
            call selci_sread(iflciv,q(ici+iroot*ncinew),nci)
 20      continue
c
         rewind iflciv
         write(iflciv) nroot, ncinew
         call selci_swrite(iflciv,q(iroots),nroot)
         do 30 iroot = 0,nroot-1
            call selci_swrite(iflciv,q(ici+iroot*ncinew),ncinew)
 30      continue
      endif
c
c     that's it folks
c
 1000 continue
c
c delete restart tape calc completed and restart not needed.
c
      call restk
c
 1010 if (ga_nodeid() .eq. 0) then
         close(iflciv,status='keep')
         close(iflcon,status='keep')
         junk = selci_mstats(0)
         call selci_gettim(endc, endw)
         write(luout,1) 'SELECT', endc-startc, endw-startw
 1       format(/1x,a,': Total cpu',f10.1,' : Total wall',f10.1)
         call util_flush(luout)
      endif
c     
      junk = selci_mfree(int12)
c
      call selci_check(ga_nodeid())
c
      call ga_sync()
c
      end
C>
C> @}
