      program conf
c
#include "implicit.fh"
#include "mptr.fh"
c
      parameter (lword = 1250000,iqoff=1)
      dimension q(lword)
c
      call pbeginf()
      call selci_gettim(startc, startw)
      junk = selci_minit(lword)
c
      call selci_confmain(q(iqoff+1), lword)
c
      if (nodeid() .eq. 0) then
         junk = selci_mstats(0)
         call selci_gettim(endc, endw)
         write(6,1) 'CONF', endc-startc, endw-startw
 1       format(/1x,a,': Total cpu',f10.1,' : Total wall',f10.1)
         call flush(6)
      endif
c
      call pend()
c
      end
      subroutine selci_confmain(q, lword)
#include "implicit.fh"
#include "cconf.fh"
#include "mptr.fh"
      dimension q(lword)
c
c     Crude version of a configuration generation program
c     Needs much doing to it to be very useful
c
      character*80 title,ciconf, confin
c
      if (nodeid().ne.0) goto 99999
c
      confin = 'confin'
      ciconf = 'ciconf'
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
      open(5,file=confin,status='old',form='formatted')
      call inibit
      call selci_start(title, q(irefo))
c
c     apply the refgens to the reference functions
c
      call selci_refgen(q, q(irefo),q(jrefo))
c
c     screen the references for symmetry and compute maximum
c     open shell dimension, the index vector and the ci dimension
c
      call selci_chkci(q(irefo),q(indxci))
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
99999 call synch(911)
c
      end
      subroutine selci_start(title, irefo)
#include "implicit.fh"
#include "cconf.fh"
      character*80 title
      integer irefo(nrefmx)
      parameter (ndd=8)
      dimension ydd(ndd),itemp(16)
      data ydd/'titl','nele','mult','symm','conf','refg','inte',
     $     'ente'/
c
c     read in title, nelec, multi, symmetry, norbs, nbpsy, conf, refgen
c
      call setinp(5,6)
      write(6,1)
 1    format(//'             Configuration Generator v1.0 8/31/89'/
     &         '             -----------------------'//)
c
 10   call input
c
 20   call inpa4(ytest)
 21   i=locatc(ydd,ndd,ytest)
      goto (110,120,130,140,150,160,170,180) i
      call errquit('conf: unrecognised directive',i)
c
c     title
c
 110  call input
      call getlin(title)
      goto 10
c
c     nelectrons
c
 120  call inpi(nelec)
      goto 10
c
c     multiplicity
c
 130  call inpi(multi)
      goto 10
c
c     symmetry ... make nbpsy, isym, norbs, nelpi, nintpo
c
 140  if (inttyp.eq.-1)
     $     call errquit('conf: integrals must precede symmetry',inttyp)
      call inpi(nsym)
      call inpi(issss)
      issss = issss - 1
      norbs = 0
      if (inttyp.eq.2) then
c
c     for columbus style integrals read in the no. of orbitals of each s
c
         call input
         do 131 i = 1,nsym
            call inpi(nbpsy(i))
            ioff = norbs+1
            norbs = norbs + nbpsy(i)
            if (norbs.gt.255)
     $           call errquit('conf: too many orbitals',norbs)
            call ifill(nbpsy(i),i-1,isym(ioff),1)
 131     continue
         nelpi = nbitpi / 2
         nintpo = (norbs + nelpi - 1) / nelpi
      else
c
c     for the ascii, binasc and gamess integrals we need to read in
c     the symmetry of each orbital ... the input is terminated
c     by the next directive
c
         call ifill(8,0,nbpsy,1)
 132     call input
c     is this the end of the input ?
         call inpa4(ytest)
         i=locatc(ydd,ndd,ytest)
         if (i.ne.0) goto 134
         call setfld(ifield()-1)
         do 133 i = 1,nfield()
            call inpi(itest)
            if (itest.lt.1 .or. itest.gt.nsym)
     $           call errquit('conf: invalid orbital symmetry ',itest)
            norbs = norbs + 1
            if (norbs.gt.255)
     $           call errquit('conf: too many orbitals ',norbs)
            isym(norbs) = itest - 1
            nbpsy(itest) = nbpsy(itest) + 1
 133     continue
         goto 132
 134     nelpi = nbitpi / 2
         nintpo = (norbs + nelpi - 1) / nelpi
         goto 21
      endif
c
      goto 10
c
c     conf <nosort>
c
 150  if (norbs.eq.-1) call errquit('conf: symm must precede conf',-1)
      if (nelec.eq.-1) call errquit('conf: nele must precede conf',-1)
      if (multi.eq.-1) call errquit('conf: mult must precede conf',-1)
      if (nintpo.le.0)call errquit('conf: nintpo invalid value',nintpo)
      osort = .true.
      call inpa4(ytest)
      if (ytest.eq.'noso') osort = .false.
      nref = 0
 151  call input
c     are we at the end of this directive
      call inpa4(ytest)
      if (ytest.eq.' ') goto 151
      i=locatc(ydd,ndd,ytest)
      if (i.ne.0) goto 21
      call setfld(ifield()-1)
c
      call ifill(norbs,0,iocc,1)
      call inpi(ns)
      if (mod(ns,2).ne.mod(multi-1,2) .or. ns.lt.multi-1)
     $     call errquit('conf: ns and multi do not match',ns)
      do 152 j = 1,ns
         call inpi(iorb)
         iocc(iorb) = 1
 152  continue
      do 153 j = ifield()+1,nfield()
         ns = ns + 2
         call inpi(iorb)
         iocc(iorb) = 3
 153  continue
      if (ns.ne.nelec) call errquit('conf: wrong no. of electrons',ns)
      if (nref.eq.0) then
         nref = 1
         call selci_pkcon(norbs, iocc, irefo, nintpo, nbitpi)
      else
         if (osort) then
            call selci_pkcon(norbs, iocc, itemp, nintpo, nbitpi)
            call selci_bserch(nref, itemp, irefo, nintpo, ipos)
            if (ipos.le.0)
     $           call errquit('conf: replication in conf',ipos)
            call selci_insert(nref, itemp, irefo, nintpo, ipos)
         else
            ipt = nref*nintpo+1
            if (ipt+nintpo-1.gt.nrefmx)
     $           call errquit('conf: too many ref',nref+1)
            call selci_pkcon(norbs, iocc, irefo(ipt), nintpo, nbitpi)
            nref = nref + 1
         endif
      endif
      ipt = ipt + nintpo
      goto 151
c
c     refgen
c
 160  call inpa(ytest)
      if (ytest.eq.'exci') then
         call inpi(iexcit)
         goto 10
      endif
      ngen = ngen + 1
      if (ngen.gt.ngenmx)
     $     call errquit('conf: too many refgen sets',ngen)
      nigen(ngen) = 0
c
 161  call input
c     are we at the end of this directive
      call inpa4(ytest)
      if (ytest.eq.' ') goto 161
      i=locatc(ydd,ndd,ytest)
      if (i.ne.0) goto 21
      call setfld(ifield()-1)
      if (mod(nfield(),2).ne.0)
     $     call errquit('conf: odd no. of values in refgen line',
     $     nfield())
      do 162 i=1,nfield()/2
         nigen(ngen) = nigen(ngen) + 1
         call inpi(igen(1,nigen(ngen),ngen))
         call inpi(igen(2,nigen(ngen),ngen))
 162  continue
      goto 161
c
c     integrals gamess (assumes 1e on ed3 1, section 171. 2e on ed6 1)
c     integrals ascii
c     integrals binasc
c     integrals columbus
c
 170  call inpa4(ytest)
      if (ytest.eq.'asci') then
         inttyp = 0
      else if (ytest.eq.'game') then
         call errquit('selci_conf: gamess integrals no longer',0)
         inttyp = 1
      else if (ytest.eq.'colu') then
         inttyp = 2
      else if (ytest.eq.'bina') then
         inttyp = 3
      else
         call errquit('conf: unknown integral source specified',inttyp)
      endif
      goto 10
c
c     enter
c
 180  continue
      if (nelec.le.0)
     $     call errquit('conf: invalid no. electrons',nelec)
      if (multi.le.0 .or. (mod(multi-1,2).ne.mod(nelec,2)))
     $     call errquit('conf: multi and nelec do not agree',multi)
      if (norbs.le.0)
     $     call errquit('conf: invalid norbs',norbs)
      if (nsym.le.0 .or. nsym.gt.8)
     $     call errquit('conf: ivalid nsym',nsym)
      if (nref.le.0)
     $     call errquit('conf: no reference functions',nref)
      if (issss+1.le.0 .or. issss.gt.nsym)
     $     call errquit('conf: invalid state symmetry',issss)
c
c     work out dimension of the full ci
c
      fci = selci_fcidim(nelec, multi, issss, norbs, isym)
c     fci = 0.0d0
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
#include "cconf.fh"
#include "mptr.fh"
      dimension q(*)
      dimension irefo(nrefmx), jrefo(nrefmx)
      dimension itemp(16)
c
c     apply the refgen sets to the list of functions
c
      write(6,1)
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
                  call errquit('refgen: odd iocc(ikil)',iocc(ikil))
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
                  call errquit('refgen: odd iocc(iadd)',iocc(iadd))
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
     $              call errquit('refgen: too many confs',nref+1)
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
      subroutine selci_chkci(irefo, indxci)
#include "implicit.fh"
#include "cconf.fh"
      character*80 wmatrx
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
      call flush(6)
      wmatrx = 'wmatrix'
      call selci_couple(multi,nsneed(1),0,wmatrx)
c
c     read the info off the wmatrix file
c
      call ifill(33,0,nf,1)
      open(2,form='formatted',status='old',file=wmatrx)
      call selci_rdhwmt(2,mmulti, nsmax, nf, nfmax, nfmax2)
      close(2,status='delete')
      if (multi.ne.mmulti)
     $     call errquit('conf: multi of wmatrx incorrect',mmulti)
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
               call selci_wrtcon(ifllog, iocc, norbs)
               call errquit('conf: a reference has invalid ns',i)
            endif
            ipt = ipt + nintpo
            nref = nref + 1
            nnsmax = max(ns,nnsmax)
            if (nnsmax.gt.nsmax)
     $           call errquit('conf: wmatrx nsmax small, need:',nnsmax)
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
            call errquit('excit: ipt+nintpo>nrefmx',ipt+nintpo)
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
      block data
#include "implicit.fh"
#include "cconf.fh"
c
      data nelec/-1/,multi/-1/,norbs/-1/,nsym/-1/,nbpsy/8*0/,nref/-1/,
     $     issss/-1/,ngen/0/,nbitpi/32/,nintpo/-1/,iexcit/0/,
     $     inttyp/-1/
c
      end
