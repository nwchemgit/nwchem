C> \file davids.F
C> The Davidson diagonalizer
C>
C> \ingroup selci
C> @{
C>
      subroutine selci_davids(q, lword, nroot, thresh, roots)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "cbuff.fh"
#include "mptr.fh"
#include "cselcifiles.fh"
#include "global.fh"
#include "stdio.fh"
      parameter (iflham = 1, iflhdg=2, iflcon=3, iflciv=4, maxvec = 10)
c
      dimension q(lword)
      dimension roots(*)        ! Return eigenvalues to calling code
c
      dimension nbpsy(8), isym(255), nsneed(3)
**      character*10 carg
      character*80 title
c
c     diagonalize the hamiltonian. Off diagonal elements are
c     stored on file hamil, the diagonals on hamdg.
c
c     open files and read header info
c
      call selci_gettim(startc, startw)
      junk = selci_minit(lword)
      me = ga_nodeid()
      nproc = ga_nnodes()
      if (me.eq.0) then
        open(iflham,status='old',form='unformatted',file=hamil)
        open(iflhdg,status='old',form='unformatted',file=hamdg)
        open(iflcon,status='old',form='unformatted',file=ciconf)
        open(iflciv,file=civecs,form='unformatted',status='unknown')
        call selci_rdhcon(iflcon, title, multi, nelec, issss, norbs,
     &       nnsmax, nnci, noconf, nintpo, nbitpi, nbpsy, isym, nsym,
     $       inttyp, nsneed)
        read (iflham) nci
        if (nci.ne.nnci) call errquit('nci diff hamil,ciconf',nci,
     &       UNKNOWN_ERR)
      else
        open(iflham,status='old',form='unformatted',file=hamil)
      endif
      call ga_brdcst(1, nnci, mitob(1), 0)
      call ga_brdcst(2, nci,  mitob(1), 0)
      call ga_brdcst(3, noconf,  mitob(1), 0)
      call ga_brdcst(4, nintpo,  mitob(1), 0)
      call ga_brdcst(5, nbitpi,  mitob(1), 0)
c
c     how many roots are wanted and what threshold? ... now from 
c     arguments ... historically process zero read them in from disk
c
      call ga_brdcst(10, nroot, mitob(1), 0)
      call ga_brdcst(11, thresh, mdtob(1), 0)
c
      if (me.eq.0) then
        write(luout,2) title,multi,nelec,issss+1, norbs, nci, noconf,
     $       nroot, thresh
 2      format(//' Davidson CI diagonalisation '/
     &           ' ---------------------------'//
     &         1x,a//
     &         ' state multiplicity      ',i7/
     &         ' no. of electrons        ',i7/
     &         ' state symmetry          ',i7/
     &         ' no. of orbitals         ',i7/
     &         ' no. of configurations   ',i7/
     &         ' no. of orbital conf.    ',i7/
     $         ' no. of roots sought     ',i7/
     $         ' convergence           ',d9.2/)
      endif
      if (nroot.lt.0 .or. nroot.gt.nci) call errquit(' nroot ? ',nroot,
     &       UNKNOWN_ERR)
      if (thresh.lt.0.0d0) call errquit(' thresh -ve? ',0, UNKNOWN_ERR)
c
c     allocate core
c
      nrm = nroot*maxvec
      nnrm = nrm*(nrm+1)/2
      iwork = selci_mptr( max(nnrm+2*nroot*maxvec, nroot*nci) )
      iroots = selci_mptr(nroot)
      ihd = selci_mptr(nci)
      ialpha = selci_mptr(nrm*nrm)
      iredh = selci_mptr(nnrm)
      iconv = selci_mptr(nroot)
      ide = selci_mptr(nroot)
      ici = selci_mptr(nroot*nci*maxvec)
      isi = selci_mptr(nroot*nci*maxvec)
      indxci = selci_mptr(mitod(noconf+1))
      ioconf = selci_mptr(mitod(noconf*nintpo))
c
c     call actual work routine
c
      call selci_ddavid(q(iwork), q(iroots), q(ihd),q(ialpha),q(iredh),
     $     q(iconv), q(ide), q(ici), q(isi), q(indxci), q(ioconf),
     $     nroot, nci, noconf, nintpo, norbs, nbitpi, thresh)
c
      call dcopy(nroot, q(iroots), 1, roots, 1) ! Return eigenvalues
c
      close(iflham,status='keep')
      if (ga_nodeid() .eq. 0) then
         close(iflciv,status='keep')
         junk = selci_mstats(0)
         call selci_gettim(endc, endw)
         write(luout,1) 'DAVIDS', endc-startc, endw-startw
 1       format(/1x,a,': Total cpu',f10.1,' : Total wall',f10.1)
         call util_flush(luout)
      endif
      call ga_sync()
c
      end
      subroutine selci_ddavid(work, roots, hd, alpha, redh,
     $     conv, de, ci, si, indxci, ioconf, nroot, nci, noconf, nintpo,
     $     norbs, nbitpi, thresh)
#include "implicit.fh"
#include "cbuff.fh"
#include "mptr.fh"
#include "global.fh"
#include "stdio.fh"
      parameter (iflham = 1, iflhdg=2, iflcon=3, iflciv=4, maxvec = 10)
c
      dimension work(*), roots(nroot), hd(nci), alpha(*), redh(*),
     $     conv(nroot), de(nroot), ci(nci, nroot, maxvec),
     $     si(nci, nroot, maxvec), indxci(noconf+1),
     $     ioconf(nintpo,noconf)
 
c
c     read in the diagonals and the configuration info
c
      me = ga_nodeid()
      if (me.eq.0) then
        call selci_sread(iflhdg,hd,nci)
        close(iflhdg)
*        write(6,*)' hd(1)   = ',hd(1),' ',ga_nodeid()
*        write(6,*)' hd(782) = ',hd(782),' ',ga_nodeid()
        call selci_rdconf(iflcon,ioconf,indxci,nintpo,noconf)
        close(iflcon)
        indxci(noconf+1) = nci
c
c     make the starting guesses
c
*        write(6,*) ' Diagonal elements '
*        call selci_analci(hd,indxci,ioconf,nci,noconf,
*     $       norbs,nintpo,nbitpi,1)
        call selci_makges(nci,nroot,ci,hd,iflciv,roots)
        do 4 ir = 1,nroot
           write(luout,*)
           write(luout,*) ' Guess vector ',ir
           write(luout,*)
           call selci_analci(ci(1,ir,1),indxci,ioconf,nci,noconf,
     $          norbs,nintpo,nbitpi,0)
 4      continue
      endif
c
c     iterate the ci
c
      iter = 1
      itermx = 150
      nvec = 1
 20   continue
      if(ga_nodeid().eq.0)call util_flush(luout)
      call ga_brdcst(101, ci(1,1,nvec), mdtob(nroot*nci), 0)
      call selci_sigma(iflham,nci,nroot,hd,ci(1,1,nvec),si(1,1,nvec))
c      do 41 ir = 1,nroot
c         write(6,*)
c         write(6,*) ' sigma vectors ',ir,nvec
c         write(6,*)
c         call selci_analci(si(1,ir,nvec),indxci,ioconf,nci,noconf,
c     $        norbs,nintpo,nbitpi,0)
c 41   continue
c      call output(si(1,1,nvec),1,nci,1,nroot,nci,nroot,1)
      if (me.eq.0) then
        call selci_davidd(iter,itermx,nvec,maxvec,nci,nroot,ci,si,hd,
     &   alpha,redh,work,roots,conv,de,.true.)
        ndone = 0
        do 30 i = 1,nroot
           if (conv(i).lt.thresh) ndone = ndone + 1
 30     continue
      endif
      call ga_brdcst(102, ndone, mitob(1), 0)
      if (ndone.lt.nroot) goto 20
c
c     converged ... make solution vectors in work
c
      if (me.eq.0) then
        nvec = nvec - 1
        call selci_mxma(ci,1,nci,alpha,1,nroot*nvec,work,1,nci,
     $          nci,nroot*nvec,nroot)
c
c     analyse the eigenvectors
c
        do 50 iroot = 1,nroot
           write(luout,7) iroot,roots(iroot)
 7         format(//' Root ',i2,' final energy = ',f19.10/)
           call selci_analci(work((iroot-1)*nci+1),indxci,ioconf,
     $          nci,noconf, norbs,nintpo,nbitpi,0)
 50     continue
c
c     write the solution vectors to disk
c
        call selci_wrtciv(iflciv,nroot,nci,roots,work)
      endif
c
      end
      subroutine selci_wrtciv(iflciv,nroot,nci,roots,ci)
#include "implicit.fh"
      dimension roots(nroot),ci(nci,nroot)
c
      rewind iflciv
      write(iflciv) nroot,nci
      call selci_swrite(iflciv,roots,nroot)
      do 10 iroot = 1,nroot
         call selci_swrite(iflciv,ci(1,iroot),nci)
 10   continue
c
      end
      subroutine selci_makges(nci, nroot, ci, hd, iflciv, roots)
      implicit real*8 (a-h,o-z)
      integer selci_idmin
#include "errquit.fh"
#include "stdio.fh"
      dimension ci(nci,nroot), hd(nci), temp(50), itemp(50),
     $     roots(nroot)
c
      if (nroot.gt.50) call errquit(' nroot rather large!',nroot,
     &       INPUT_ERR)
c
c     are there any trial vectors sitting on disk from
c     a previous (related?) calculation?
c
      rewind iflciv
      read(iflciv,end=5) nnroot, nnci
      if (nnroot.ge.nroot .and. nnci.eq.nci) then
         write(luout,1)
 1       format(/' Using trial vectors from disk '/)
         call selci_sread(iflciv,roots,nroot)
         write(luout,2) (roots(i),i=1,nroot)
 2       format(' Eigenvalues from disk '//10(1x,5f19.8/))
         do 3 i = 1,nroot
            call selci_sread(iflciv, ci(1,i), nci)
 3       continue
         rewind iflciv
         return
      else
         rewind iflciv
      endif
c
 5    continue
      do 10 i = 1,nroot
         itemp(i) = selci_idmin(nci,hd,1)
         temp(i) = hd(itemp(i))
         hd(itemp(i)) = 1.0d20 ! Use a big number so this element will not
                               ! be found again!
 10   continue
c
      call dfill(nroot*nci,0.0d0,ci,1)
c     Running the next loop backwards ensures that the original values
c     of the diagonal will be restored no matter what! (Including asking
c     for more roots than there are bound states).
      do 20 i = nroot,1,-1
         hd(itemp(i)) = temp(i)
         ci(itemp(i),i) = 1.0d0
 20   continue
c
      end
      subroutine selci_davidd(iter,itermx,nvec,maxvec,n,nroot,ci,si,hd,
     &     alpha,redh,work,roots,conv,de,oprint)
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
#include "stdio.fh"
c
c     multi-root in-core davidson routine
c
c     iter = iteration no. ... set to 1 or less on first entry
c     itermx = max iters
c     nvec = no. of active vectors ... update selci_sigma(*,*,nvec) on exit
c     maxvec = max no. of vectors to hold
c     n = length of ci vectors
c     nroot = no. of roots required
c     ci = ci vectors
c     si = sigma vectors
c     hd = diagonals of hamiltonian ... or anything u fancy
c     alpha = solution coefficients ... size (nroot*maxvec)**2
c     redh = reduced hamiltonian ... size triangle of (nroot*maxvec)
c     work = workspace...size=max(triangle of nroot*maxvec + 2*maxvec,
c     nroot*n)
c     roots = returns roots
c     conv = returns norm of update vectors
c     de = returns incremental change in roots
c
      dimension ci(n,nroot,maxvec),si(n,nroot,maxvec),hd(n),
     &     alpha(maxvec*nroot*maxvec*nroot),redh(*),work(*),
     &     roots(nroot),conv(nroot),de(nroot)
      logical oprint
c
      if (iter.le.1) then
         iter = 1
         nvec = 1
         if (maxvec.lt.2)
     $        call errquit(' stupid maxvec in davids ',maxvec,
     &       INPUT_ERR)
         if (itermx.lt.2)
     $        call errquit(' stupid itermx in davids ',itermx,
     &       INPUT_ERR)
         call dfill(nroot,0.0d0,roots,1)
      endif
c
c     special case n = 1
c
      if (n.eq.1) then
         if (nroot.gt.1)
     $        call errquit(' n=1, nroot > 1 in davids ',nroot,
     &       INPUT_ERR)
         conv(1) = 0.0d0
         roots(1) = si(1,1,1)
         alpha(1) = 1.0d0
         de(1) = si(1,1,1)
         nvec = 2
         nnvec = 1
         goto 100
      endif
c
c     update the reduced hamiltonian and diagonalise it
c
      nnvec = nvec
 3    nrnv1 = nroot * (nvec-1)
      ih = nrnv1*(nrnv1+1)/2
      do 5 jroot = 1,nroot
         do 10 ivec = 1,nvec
            itop = nroot
            if (ivec.eq.nvec) itop = jroot
            do 20 iroot = 1,itop
               ih = ih + 1
               redh(ih) = ddot(n,ci(1,iroot,ivec),1,
     &              si(1,jroot,nvec),1)
 20         continue
 10      continue
 5    continue
c
*      write(6,*) ' reduced h matrix '
*      call selci_prtri(redh,nvec*nroot)
      nrnv = ih
      iredh = 1
      ibig = iredh + nrnv
      iibig = ibig + nvec*nroot
c
      call dcopy(nrnv,redh,1,work(iredh),1)
      call selci_yacobi(work(iredh),alpha,nroot*nvec,
     $     work(ibig),work(iibig))
*      write(6,*) ' reduced h e-vectors'
*      call output(alpha,1,nroot*nvec,1,nroot*nvec,nroot*nvec,
*     $     nroot*nvec,1)
*      write(6,*) ' reduced h evals '
*      call output(work(ibig),1,nroot*nvec,1,1,nroot*nvec,1,1)
      do 25 iroot = 1,nroot
         rnew = work(ibig-1 + iroot)
         if (nvec.eq.nnvec) then
            de(iroot) = rnew - roots(iroot)
         endif
         roots(iroot) = rnew
 25   continue
c
c     if nvec = maxvec then take lowest nroot eigen vectors
c     as the new expansion set, put nvec = 1 and rediagonalise
c
      if (nvec.eq.maxvec) then
         call selci_mxma(ci,1,n,alpha,1,nroot*nvec,work,1,n,
     $        n,nroot*nvec,nroot)
         call dcopy(n*nroot,work,1,ci(1,1,1),1)
         call selci_mxma(si,1,n,alpha,1,nroot*nvec,work,1,n,
     $        n,nroot*nvec,nroot)
         call dcopy(n*nroot,work,1,si(1,1,1),1)
         nvec = 1
         goto 3
      endif
c
c     now form the updates
c
      call selci_mxma(ci,1,n,alpha,1,nroot*nvec,work,1,n,n,
     $     nroot*nvec,nroot)
      do 30 iroot = 1,nroot
         ioff = (iroot-1)*n + 1
         call dscal(n,-roots(iroot),work(ioff),1)
 30   continue
c      call dgemm('n','n',n, nroot, nroot*nvec, 1.0d0, si, n,
c     *     alpha, nroot*nvec, 1.0d0, work, n)
      call selci_mxmb(si,1,n,alpha,1,nroot*nvec,work,1,n,n,
     $     nroot*nvec,nroot)
      do 40 iroot = 1,nroot
         ioff = (iroot-1)*n
         s = 0.0d0
         do 50 i = 1,n
            ediff = roots(iroot) - hd(i)
            if (abs(ediff).lt.0.1d0) ediff = sign(0.1d0,ediff)
            work(ioff + i) = work(ioff + i)/ediff
            s = s + work(ioff + i)*work(ioff + i)
 50      continue
         conv(iroot) = sqrt(s)
 40   continue
c
c     work now contains the non-orthogonalised update vectors
c
      nvec = nvec + 1
      do 60 iroot = 1,nroot
         ioff = (iroot-1)*n + 1
         call selci_renorm(n,work(ioff),1,dw)
         if (dw.eq.0.0d0) goto 60
         call dcopy(n,work(ioff),1,ci(1,iroot,nvec),1)
         npass = 0
 65      continue
         do 70 jvec = 1,nvec
            jtop = nroot
            if (jvec.eq.nvec) jtop = iroot-1
            do 80 jroot = 1,jtop
               s = ddot(n,ci(1,iroot,nvec),1,ci(1,jroot,jvec),1)
               call daxpy(n,-s,ci(1,jroot,jvec),1,ci(1,iroot,nvec),1)
               call selci_renorm(n,ci(1,iroot,nvec),1,dci)
               if (dci.lt.1.0d-2) then
                  npass = npass + 1
                  if (npass.le.2 .and. dci.ne.0.0d0) goto 65
c     can't othogonalise to other vectors ... linear dependance
c     may be due to exhaustion of expansion space or error
c     print out message, then resum and proceed from there
                  write(luout,81)
 81               format(/' orthogonalisation failure ... if space',
     &                 ' exhausted then OK '/)
                  call selci_mxma(ci,1,n,alpha,1,nroot*nnvec,work,1,n,n,
     $                 nroot*nnvec,nroot)
                  call dcopy(n*nroot,work,1,ci(1,1,1),1)
                  call selci_mxma(si,1,n,alpha,1,nroot*nnvec,work,1,n,n,
     $                 nroot*nnvec,nroot)
                  call dcopy(n*nroot,work,1,si(1,1,1),1)
                  nvec = 1
                  goto 3
c                  do 85 irr = 1,nroot
c                     conv(irr) = 0.0d0
c 85               continue
c                  goto 90
               endif
 80         continue
 70      continue
 60   continue
c
 90   continue
c
 100  continue
      if (oprint) then
         write(luout,1) iter,nnvec,(roots(i),de(i),conv(i),i=1,nroot)
 1       format(' iter=',i2,'  nvec=',i2,
     &        10('  root=',f15.8,'  de=',f15.8,'  conv=',f13.8:/17x))
         call util_flush(luout)
      endif
c
      iter = iter + 1
      if (iter.ge.itermx+1) 
     $     call errquit(' too many iters, davids ',iter,
     &       INPUT_ERR)
c
      end
      subroutine selci_analci(ci,indxci,ioconf,nci,noconf,
     $     norbs,nintpo,nbitpi,iorder)
      implicit real*8 (a-h,o-z)
      integer selci_idmin
#include "errquit.fh"
      parameter (maxtmp = 20)
      dimension ci(nci),indxci(noconf+1),ioconf(nintpo,noconf),
     $     temp(maxtmp),itemp(maxtmp)
      dimension irange(21)
      integer selci_isum
c
c     print out maxtmp largest non-zero components of the ci vector
c     along with orbital occupation information etc
c
c     iorder = 0 -> select by abs max components and also print
c                   statistics about size of components
c     iorder = 1 -> select by min components
c
      myiorder = iorder
      if (myiorder.ne.0 .and. myiorder.ne.1) then
         write(6,*) ' myiorder odd in analci ... default to zero '
         myiorder = 0
      endif
      write(6,1)
 1    format(1x,'    coeff.      conf  spin     conf norm'/
     &       1x,'------------  ------  ----  ------------')
      ndo = min(maxtmp,nci)
      ngot = 0
c
c     find the ndo largest coefficients
c
      do 10 i = 1,ndo
         if (myiorder.eq.0) then
            it = idamax(nci,ci,1)
         else
            it = selci_idmin(nci,ci,1)
         endif
         temp(i) = ci(it)
         itemp(i) = it
c
         if (temp(i).eq.0.0d0) goto 50
         ngot = ngot + 1
         ci(it) = 0.0d0
 10   continue
c
 50   do 40 i = 1,ngot
         ci(itemp(i)) = temp(i)
 40   continue
c
      do 60 i = 1,ngot
         it = itemp(i)
         do 20 io = 1,noconf+1
            if (it.le.indxci(io)) goto 30
 20      continue
         call errquit('analci: fell out of do 20',io, INPUT_ERR)
 30      iconf = io-1
         ispin = it - indxci(iconf)
         nf = indxci(iconf+1)-indxci(iconf)
         cnorm = dnrm2(nf,ci(indxci(iconf)+1),1)
         call selci_prtanl(iconf,ispin,temp(i),cnorm,norbs,
     $        ioconf(1,iconf), nintpo,nbitpi)
 60   continue
c
      if (myiorder.eq.0) then
         call ifill(21,0,irange,1)
         do 70 i = 1,nci
            test = abs(ci(i))
            if (test.gt.0.0d0) then
               ihash = 1 - nint(2.0d0*log10(test))
               ihash = min(ihash,21)
               ihash = max(ihash,1)
               irange(ihash) = irange(ihash) + 1
            endif
 70      continue
         write(6,101) (10.0d0**(dble(1-i)*0.5d0), irange(i),
     $        i=1,21)
 101     format(/
     $     ' Analysis of CI Coefficients by Size '/
     $     ' ----------------------------------- '//
     $     '         Value            Count '/
     $     '     -------------        ----- '/
     $     21 (:5x,f13.10, 5x,i8/))
         write(6,1011) selci_isum(21,irange,1)
 1011    format(5x,13x,4x,'---------'/5x,13x,4x,i9)
      endif
c
      end
      subroutine selci_prtanl(iconf,ispin,c,d,norbs,ioconf,
     $     nintpo,nbitpi)
      implicit real*8 (a-h,o-z)
      dimension ioconf(nintpo),iocc(255),list(255)
c
      call selci_upkcon(norbs, iocc, ioconf, nintpo, nbitpi)
      igot = 0
      do 10 i = 1,norbs
         if (iocc(i).gt.0) then
            igot = igot + 1
            list(igot) = i
            if (iocc(i) .eq. 1) then
               iocc(igot) = 1
            else
               iocc(igot) = 2
            endif
         endif
 10   continue
c
      if (norbs.ge.100) then
         write(6,1) c,iconf,ispin,d,(list(i),iocc(i),i=1,igot)
 1       format(1x,f12.6,2x,i6,2x,i4,2x,f12.6,4x,9(i4,'(',i1,')':))
      else
         write(6,2) c,iconf,ispin,d,(list(i),iocc(i),i=1,igot)
 2       format(1x,f12.6,2x,i6,2x,i4,2x,f12.6,4x,11(i3,'(',i1,')':))
      endif
c
      end
C>
C> @}
