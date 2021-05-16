C> \ingroup selci
C> @{
C>
C> \brief The configuration generator
C>
C> A crude version of a configuration generator.
C>
      subroutine selci_conf(rtdb, q, lword)
      implicit none
*
* $Id$
*
*include "implicit.fh"
#include "cconf.fh"
#include "mptr.fh"
#include "cselcifiles.fh"
#include "rtdb.fh"
#include "global.fh"
c
      integer rtdb  !< [Input] The RTDB handle
      integer lword !< [Input] The length of the main memory array
      double precision q(lword) !< [Scratch] Work space
c
      integer mdtoi, mitod
      external mdtoi, mitod
      external data_selci_conf  ! For T3D to link block data
c
c     Local
c
      integer junk
      integer indxci, irefo, jrefo
c
      double precision startc, startw
      double precision endc,   endw
c
c     Crude version of a configuration generation program
c     Needs much doing to it to be very useful
c
      character*80 title
      logical status
c
      if (ga_nodeid().ne.0) goto 99999
      status = rtdb_parallel(.false.)
      
      call selci_gettim(startc, startw)
      junk = selci_minit(lword)
c
c     Leave a quarter of the space for excit to function
c     (we could be more intelligent but I have a bad cold so
c     thinking hurts)
c
      nrefmx = mdtoi(lword)/4 + 1
      irefo  = selci_mptr(mitod(nrefmx))
      jrefo  = selci_mptr(mitod(nrefmx))
      indxci = selci_mptr(mitod(nrefmx+1))
c
      call selci_start(rtdb, title, q(irefo))
c
c     apply the refgens to the reference functions
c
      call selci_refgen(q, q(irefo),q(jrefo))
c
c     screen the references for symmetry and compute maximum
c     open shell dimension, the index vector and the ci dimension
c
      call selci_chkci(q(irefo),q(indxci),wmatrx)
c
c     write out the header for the conf file
c
      open(1,form='unformatted',status='unknown',file=ciconf)
      call selci_wthcon(1, title, multi, nelec, issss, norbs,
     &     nnsmax, nci, nref, nintpo, nbitpi, nbpsy, isym, nsym,
     &     inttyp, nsneed)
c
      call selci_wtconf(1,q(irefo),nintpo*nref,q(indxci),nref)
      close(1,status='keep')
c
      junk = selci_mfree(irefo)
c
      status = rtdb_parallel(.true.)
99999 call ga_sync()
      if (ga_nodeid() .eq. 0) then
         junk = selci_mstats(0)
         call selci_gettim(endc, endw)
         write(6,1) 'CONF', endc-startc, endw-startw
 1       format(/1x,a,': Total cpu',f10.1,' : Total wall',f10.1)
         call util_flush(6)
      endif
c
      end
C>
C> \brief Construct the configuration table
C>
      subroutine selci_start(rtdb,title, irefo)
#include "implicit.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "cconf.fh"
#include "cselcifiles.fh"
#include "bitops.fh"
      integer rtdb
      character*80 title
      integer irefo(nrefmx)
      parameter (ndd=8)
      dimension itemp(16)
      double precision occ(255)
      character*80 name
      parameter (max_inconf = 4095)
      integer inconf(max_inconf)
c
      character*10 intfmt
c
      logical cckohn_ci
c
      if (.not. rtdb_get(rtdb,'selci:cckohn_ci',mt_log,1,cckohn_ci))
     $   cckohn_ci=.false.
c
c     read in title, nelec, multi, symmetry, norbs, nbpsy, conf, refgen
c
      write(6,1)
 1    format(//'             Configuration Generator v1.0 8/31/89'/
     &         '             -----------------------'//)
c
c     Get most junk from the moints file
c
      title = ' '
c
      if (.not. rtdb_cget(rtdb,'fourindex:file format',1,intfmt))
     $   intfmt = 'nwchem'
c
      if (intfmt .eq. 'ascii') then
        inttyp = 0
      else if (intfmt .eq. 'nwchem') then
        inttyp = 3
      else
        call errquit('selci:conf:unknown integral file format',0,
     &       INPUT_ERR)
      endif
c
      iflint = 8
      if (intfmt .eq. 'ascii') then
        open(iflint, form='formatted',status='old',file=moints)
        read(iflint,'(a)') title
        read(iflint,*) norbs, nfzc
        read(iflint,*)            ! map
        read(iflint,*) (isym(i),i=1,norbs) ! syms
        read(iflint,*)            ! eigenvalues
        read(iflint,*) (occ(i),i=1,norbs) ! occupations
      else
        open(iflint, form='unformatted',status='old',file=moints)
        read(iflint) title
        read(iflint) norbs, nfzc
        read(iflint)
        read(iflint)            ! map
        read(iflint) (isym(i),i=1,norbs) ! syms
        read(iflint)            ! eigenvalues
        read(iflint) (occ(i),i=1,norbs) ! occupations
      endif
      close(iflint, status='keep')
c
c     nelectrons ... add up the occupation numbers ... 
c
      dnelec = 0.0d0
      do i = 1, norbs
         dnelec = dnelec + occ(i)
      enddo
      nelec = nint(dnelec)
c
c     multiplicity ... get this from the SCF/MCSCF last run
c
      if (.not. rtdb_get(rtdb,'scf:nopen', mt_int, 1, nopen))
     $     nopen = 0
      multi = nopen + 1
      if(cckohn_ci)then
         write(*,*)"one more electron are added to the scf case"
         multi=multi+1
         nelec=nelec+1
      endif
c
c     symmetry ... make nbpsy, isym, norbs, nelpi, nintpo
c
      nsym = 0
      do is = 0, 7
         nbpsy(is+1) = 0
         do i = 1, norbs
            if (isym(i).eq.is) nbpsy(is+1) = nbpsy(is+1)+1
         enddo
         if (nbpsy(is+1).gt.0) nsym = is+1
      enddo
      if (nsym .le. 2) then
         continue
      else if (nsym .le. 4) then
         nsym = 4
      else if (nsym .le. 8) then
         nsym = 8
      else
         call errquit('selci_conf: symmetry nonsense',nsym, INPUT_ERR)
      endif
      nelpi = nbitpi / 2
      nintpo = (norbs + nelpi - 1) / nelpi
c
c     conf <nosort> ... construct a default reference from the occupation
c     pattern heading the moints.  Eventually we'll also want to use
c     default from MCSCF.
c
      if (.not. rtdb_get(rtdb, 'selci:osort', mt_log, 1, osort))
     $     osort = .true.
      nref = 0
      call ifill(max_inconf, -1, inconf, 1)
      if (rtdb_get(rtdb, 'selci:conf', mt_int, max_inconf, inconf)) then
c
c     list on database is terminated with a 0
c
         ind = 1
 10      if (inconf(ind) .ge. 0) then
            ns = inconf(ind)
            ind = ind + 1
            if (mod(ns,2).ne.mod(multi-1,2) .or. ns.lt.multi-1)
     $           call errquit('conf: ns and multi do not match',ns,
     &       INPUT_ERR)
            call ifill(norbs, 0, iocc, 1)
            do j = 1,ns
               if (inconf(ind).le.0) 
     $              call errquit('selci: bad socc in conf',inconf(ind),
     &       INPUT_ERR)
               iocc(inconf(ind)) = 1
               ind = ind + 1
            enddo
            do j = ns+1,nelec,2
               if (inconf(ind).le.0) 
     $              call errquit('selci: bad docc in conf',inconf(ind),
     &       INPUT_ERR)
               iocc(inconf(ind)) = 3
               ind = ind + 1
            enddo
            if (nref.eq.0) then
               nref = 1
               call selci_pkcon(norbs, iocc, irefo, nintpo, nbitpi)
            else
               if (osort) then
                  call selci_pkcon(norbs, iocc, itemp, nintpo, nbitpi)
                  call selci_bserch(nref, itemp, irefo, nintpo, ipos)
                  if (ipos.le.0)
     $                 call errquit('conf: replication in conf',ipos,
     &       INPUT_ERR)
                  call selci_insert(nref, itemp, irefo, nintpo, ipos)
               else
                  ipt = nref*nintpo+1
                  if (ipt+nintpo-1.gt.nrefmx)
     $                 call errquit('conf: too many references',nref+1,
     &       INPUT_ERR)
                  call selci_pkcon(norbs,iocc,irefo(ipt),nintpo,nbitpi)
               endif
            endif
            goto 10
         endif
      else
         ns = 0
         nd = 0
         do i = 1, norbs
            iocc(i) = nint(occ(i))
            if (iocc(i) .eq. 0) then
               continue
            else if (iocc(i) .eq. 1) then
               ns = ns + 1
            else if (iocc(i) .eq. 2) then
               nd = nd + 1
               iocc(i) = 3
            else
               call errquit('conf: bad default occ?',0, INPUT_ERR)
            endif
         enddo
         if (mod(ns,2).ne.mod(multi-1,2) .or. ns.lt.multi-1)
     $        call errquit('conf: ns and multi do not match',ns,
     &       INPUT_ERR)
         if ((ns+2*nd) .ne. nelec)
     $        call errquit('conf: default nelec don not match',2*nd+ns,
     &       INPUT_ERR)
         call selci_pkcon(norbs, iocc, irefo, nintpo, nbitpi)
         nref = 1
      endif
c
c     Determine state symmetry from the first conf
c
      call selci_upkcon(norbs, iocc, irefo, nintpo, nbitpi)
      issss = 0
      do i = 1, norbs
         if (iocc(i) .eq. 1) issss = ieor(issss,isym(i))
      enddo
c
c     refgen [excit <iexcit>]
c
      if (.not. rtdb_get(rtdb, 'selci:exci',mt_int, 1, iexcit))
     $     iexcit = 0
      if (.not. rtdb_get(rtdb, 'selci:ngen',mt_int, 1, ngen))
     $     ngen = 0
      if (ngen .gt. ngenmx)
     $     call errquit('conf: too many refgen sets',ngen, INPUT_ERR)
      call ifill(2*255*ngen, -1, igen, 1)
      do jj = 1, ngen
         name = ' '
         write(name,'(''selci:refgen '',i2)') jj
         write(6,*) ' FMT = ', name
         if (.not.rtdb_get(rtdb,name,mt_int,2*255,igen(1,1,jj))) then
            write(6,*) ' name = ', name
            call errquit('conf: could not find refgen set',jj, RTDB_ERR)
         endif
         do i = 1, 255
            if (igen(1,i,jj) .le. 0) then
               nigen(jj) = i-1
               goto 160
            endif
         enddo
 160     continue
      enddo
c
c     enter
c
      if (nelec.le.0)
     $     call errquit('conf: invalid no. electrons',nelec, INPUT_ERR)
      if (multi.le.0 .or. (mod(multi-1,2).ne.mod(nelec,2)))
     $     call errquit('conf: multi and nelec do not agree',multi,
     &       INPUT_ERR)
      if (norbs.le.0)
     $     call errquit('conf: invalid norbs',norbs, INPUT_ERR)
      if (nsym.le.0 .or. nsym.gt.8)
     $     call errquit('conf: ivalid nsym',nsym, INPUT_ERR)
      if (nref.le.0)
     $     call errquit('conf: no reference functions',nref, INPUT_ERR)
      if (issss+1.le.0 .or. issss.gt.nsym)
     $     call errquit('conf: invalid state symmetry',issss, INPUT_ERR)
c
c     work out dimension of the full ci
c
      fci = selci_fcidim(nelec, multi, issss, norbs, isym)
c
      write(6,2) title,multi,issss+1,nelec,norbs,fci,nref,ngen,nsym,
     $     (nbpsy(i),i=1,nsym)
 2    format(1x,a//
     $       ' state multiplicity      ',i7/
     $       ' state symmetry          ',i7/
     $       ' no. of electrons        ',i7/
     $       ' no. of orbitals         ',i7/
     $       ' fullci dimension        ',d22.14/
     $       ' no. of input references ',i7/
     $       ' no. of refgen sets      ',i7/
     $       ' no. of symmetry types   ',i7/
     $       ' no. of orbitals per sym   ',8i5//
     $       ' symmetry of the orbitals:'/)
      write(6,3) (i,isym(i)+1,i=1,norbs)
 3    format(9(i5,'(',i1,')':))
c
      write(6,4)
 4    format(//' Input reference functions '//)
      ipt = 1
      do 185 iref = 1,nref
         call selci_prtcon(6, norbs, irefo(ipt), nintpo, nbitpi)
         ipt = ipt + nintpo
 185  continue
c
      do 190 i = 1,ngen
         write(6,5) i
 5       format(//' Creation/Annihilation operators set no.',i3/)
         if (norbs.gt.100) then
            write(6,6) (igen(1,k,i),igen(2,k,i),k=1,nigen(i))
 6          format(6('  (',i3,',',i3,')':))
         else
            write(6,7) (igen(1,k,i),igen(2,k,i),k=1,nigen(i))
 7          format(8('  (',i2,',',i2,')':))
         endif
 190  continue
c
      end
      subroutine selci_refgen(q, irefo, jrefo)
#include "implicit.fh"
#include "errquit.fh"
#include "cconf.fh"
#include "mptr.fh"
      dimension q(*)
      dimension irefo(nrefmx), jrefo(nrefmx)
      dimension itemp(16)
c
c     apply the refgen sets to the list of functions
c
      if (ngen.gt.0 .or. iexcit.gt.0) write(6,1)
 1    format(//' Apply reference generators '/)
      do 10 ig = 1,ngen
         call selci_icopy(nintpo*nref,irefo,1,jrefo,1)
c
c     loop through functions from old list
c
         nrefo = nref
         ipto = 1
         do 20 iro = 1,nrefo
            call selci_upkcon(norbs, iocc, jrefo(ipto), nintpo, nbitpi)
            ns = 0
            do 21 k = 1,norbs
               if (iocc(k).eq.1) ns = ns + 1
 21         continue
c
c     loop thru operators in this refgen set
c
            do 30 k = 1,nigen(ig)
               call selci_upkcon(norbs,iocc,jrefo(ipto),nintpo,nbitpi)
               nns = ns
               iadd = igen(1,k,ig)
               ikil = igen(2,k,ig)
c
c     apply annihilation operator
c
               if (iocc(ikil).eq.0) then
                  goto 30
               else if (iocc(ikil).eq.1) then
                  iocc(ikil) = 0
                  nns = nns - 1
               else if (iocc(ikil).eq.3) then
                  iocc(ikil) = 1
                  nns = nns + 1
               else
                  call errquit('refgen: odd iocc(ikil)',iocc(ikil),
     &       UNKNOWN_ERR)
               endif
c
c     apply creation operator
c
               if (iocc(iadd).eq.3) then
                  goto 30
               else if (iocc(iadd).eq.1) then
                  iocc(iadd) = 3
                  nns = nns - 1
               else if (iocc(iadd).eq.0) then
                  iocc(iadd) = 1
                  nns = nns + 1
               else
                  call errquit('refgen: odd iocc(iadd)',iocc(iadd),
     &       UNKNOWN_ERR)
               endif
c
c     are there enough unpaired electrons for the multiplicity
c
               if (nns.lt.multi-1) goto 30
c
c     find position in sorted list and check if there already
c
               call selci_pkcon(norbs, iocc, itemp, nintpo, nbitpi)
               call selci_bserch(nref, itemp, irefo, nintpo, ipos)
               if (ipos.le.0) goto 30
               if ((nref+1)*nintpo.gt.nrefmx)
     $              call errquit('refgen: too many confs',nref+1,
     &       INPUT_ERR)
               call selci_insert(nref, itemp, irefo, nintpo, ipos)
 30         continue
            ipto = ipto + nintpo
 20      continue
         write(6,2) ig,nref
 2       format(' After refgen ',i2,' no. of references',i6)
 10   continue
c
c     now we have applied the specific refgen sets apply any
c     excitation levels ... these apply to all orbitals
c
      if (iexcit.gt.0) then
         ioccrf = selci_mptr(mitod(norbs*nref))
         ixstak = selci_mptr(mitod((norbs+2)*nref))
         call selci_excit(irefo, q(ioccrf), q(ixstak))
         junk = selci_mfree(ioccrf)
      endif
c
c      write(6,4)
c 4    format(//' List of reference functions after refgen'//)
c      ipt = 1
c      do 50 iref = 1,nref
c         call selci_prtcon(6, norbs, irefo(ipt), nintpo, nbitpi)
c         ipt = ipt + nintpo
c 50   continue
c
      end
      subroutine selci_chkci(irefo, indxci,wmatrx)
#include "implicit.fh"
#include "errquit.fh"
#include "cconf.fh"
      character*(*) wmatrx
      integer irefo(nrefmx), indxci(nrefmx+1)
#include "bitops.fh"
c
c     check the symmetry of the reference functions and
c     also compute the maximum open shell dimension, the
c     index vector and the dimension of the ci
c
c     Generate the coupling coefficient file which contains
c     the number of spin functions for each number of open shells
c     
      call selci_mknsneed(irefo,nref,nintpo,nbitpi,norbs,nsneed)
      call util_flush(6)
      call selci_couple(multi,nsneed(1),0,wmatrx)
c
c     read the info off the wmatrix file
c
      call ifill(33,0,nf,1)
      open(2,form='formatted',status='old',file=wmatrx)
      call selci_rdhwmt(2,mmulti, nsmax, nf, nfmax, nfmax2)
      close(2,status='delete')
      if (multi.ne.mmulti)
     $     call errquit('conf: multi of wmatrx incorrect',mmulti,
     &       INPUT_ERR)
      write(6,7) (nf(i),i=mod(nsmax,2),nsmax,2)
 7    format(/' no. of spin functions    ',10i5/)
      nnsmax = 0
      nrefo = nref
      nref = 0
      ipt = 1
      ipto = 1
      nci = 0
      do 10 i = 1,nrefo
         call selci_upkcon(norbs, iocc, irefo(ipto), nintpo, nbitpi)
         is = 0
         ns = 0
c$doit asis
         do 20 j = 1,norbs
            if (iocc(j).eq.1) then
               is = ieor(is,isym(j))
               ns = ns + 1
            endif
 20      continue
         if (is.eq.issss) then
            do 30 j = 0,nintpo-1
               irefo(ipt+j) = irefo(ipto+j)
 30         continue
            if (mod(ns,2).ne.mod(multi-1,2) .or. ns.lt.multi-1) then
               call selci_wrtcon(6, iocc, norbs)
               call errquit('conf: a reference has invalid ns',i,
     &       INPUT_ERR)
            endif
            ipt = ipt + nintpo
            nref = nref + 1
            nnsmax = max(ns,nnsmax)
            if (nnsmax.gt.nsmax)
     $           call errquit('conf: wmatrx nsmax small, need:',nnsmax,
     &       INPUT_ERR)
            indxci(nref) = nci
            nci = nci + nf(ns)
         endif
         ipto = ipto + nintpo
 10   continue
      indxci(nref+1) = nci
c
      write(6,1) nref,nci
 1    format(/' no. of references after symmetry screen',i6/
     $        ' dimension of the CI is ',i7)
      write(6,4)
 4    format(//' Final list of reference functions'//)
      ipt = 1
      do 50 iref = 1,nref
         call selci_prtcon(6, norbs, irefo(ipt), nintpo, nbitpi)
         ipt = ipt + nintpo
 50   continue
c
      end
      subroutine selci_excit(irefo, ioccrf, ixstak)
      implicit integer (a-h,o-z)
#include "errquit.fh"
#include "cconf.fh"
c
c     generate orbital configurations of up to a certain level
c     of excitation from the current reference set.
c     does this by selectively walking through the full-CI
c     orbital graph
c
      integer irefo(*)
      dimension ioccrf(norbs,nref),ixstak(0:norbs+1,nref)
      integer nstak(0:255),nsstak(0:255)
#include "bitops.fh"
c
      nrefo = nref
      ipto = 1
      do 10 ir = 1,nrefo
         call selci_upkcon(norbs,ioccrf(1,ir),irefo(ipto),nintpo,nbitpi)
         ipto = ipto + nintpo
         do 20 i = 1,norbs
            if (ioccrf(i,ir).eq.3) ioccrf(i,ir) = 2
 20      continue
 10   continue
      call ifill(nref*(norbs+1),0,ixstak,1)
      call ifill(norbs+1,0,nstak,1)
      call ifill(norbs+1,0,nsstak,1)
c
      ixcit2 = iexcit*2
      iorb = 1
      iocc(1) = -1
c
c     new level ... if no such level then have finished
c
 30   if (iorb.lt.1) goto 100
      call ifill(norbs-iorb,-1,iocc(iorb+1),1)
c
c     new occupation at same level
c
 40   iocc(iorb) = iocc(iorb) + 1
      if (iocc(iorb).eq.1) then
         nsstak(iorb) = nsstak(iorb-1) + 1
      else
         nsstak(iorb) = nsstak(iorb-1)
      endif
      nstak(iorb) = nstak(iorb-1) + iocc(iorb)
      if (iocc(iorb).gt.2 .or. nstak(iorb).gt.nelec) then
         iorb = iorb - 1
         goto 30
      endif
      if (nelec-nstak(iorb).gt.2*(norbs-iorb)) goto 40
      minex = 2*nelec
      do 50 iref = 1,nrefo
         ixstak(iorb,iref) = ixstak(iorb-1,iref) +
     $        abs(iocc(iorb)-ioccrf(iorb,iref))
         minex = min(minex,ixstak(iorb,iref))
 50   continue
      if (minex.gt.ixcit2) goto 40
c
      if (iorb.eq.norbs) then
         if (nstak(iorb).ne.nelec) goto 40
         if (nsstak(iorb).lt.multi-1) then
            iorb = iorb - 1
            goto 30
         endif
c     if zero excitation then have it already
         if (minex.eq.0) then
            iorb = iorb - 1
            goto 30
         endif
         ksym = 0
         do 70 i = 1,norbs
            if (iocc(i).eq.2) iocc(i) = 3
            if (iocc(i).eq.1) ksym = ieor(ksym,isym(i))
 70      continue
c        check the symmetry
         if (ksym .ne. issss) goto 81
         nref = nref + 1
         ipt = (nref-1)*nintpo + 1
         if (ipt+nintpo.gt.nrefmx) then
            call errquit('excit: ipt+nintpo>nrefmx',ipt+nintpo,
     &       INPUT_ERR)
         else
            call selci_pkcon(norbs, iocc, irefo(ipt), nintpo, nbitpi)
         endif
 81      do 80 i = 1,norbs
            if (iocc(i).eq.3) iocc(i) = 2
 80      continue
         iorb = iorb - 1
         goto 30
      else
         iorb = iorb + 1
         iocc(iorb) = -1
         goto 30
      endif
c
 100  write(6,2) iexcit,nref
 2    format(' After excit ',i3,' no. of references',i6)
c
      end
      block data data_selci_conf
#include "implicit.fh"
#include "cconf.fh"
c
      data nelec/-1/,multi/-1/,norbs/-1/,nsym/-1/,nbpsy/8*0/,nref/-1/,
#if defined(EXT_INT)
     $     issss/-1/,ngen/0/,nbitpi/64/,nintpo/-1/,iexcit/0/,
#else
     $     issss/-1/,ngen/0/,nbitpi/32/,nintpo/-1/,iexcit/0/,
#endif
     $     inttyp/-1/
c
      end
C>
C> @}
