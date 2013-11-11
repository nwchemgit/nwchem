C> \ingroup selci
C> @{
      subroutine selci_mkvc(q, h,g,int12,int34,w1,w2,ioconf,indxci,
     $     iexv,iconf,ncold,ci,vc,iocc,ipack,numf,work1)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
#include "mptr.fh"
      dimension q(*)
      dimension h(nnorbs),g(numint),int12(nnorbs),int34(nnorbs),
     $     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1),
     $     ioconf(nintpo,ncold),indxci(ncold),iexv(ncold),
     $     ci(nci,nroot),vc(nfmax,nroot),ipack(nintpo),
     $     work1(nfmax,nfmax),iocc(255)
c
c     for orbital configuration specified by packed occupancy
c     ipack and unpacked occupancy iocc, form the interaction
c     with the nroot ci vectors, starting from orbital
c     configuration iconf. In work1 also return the diagonal
c     block of the hamiltonian involving the test configuration
c
      dimension indbar(255), lists(255), listd(255),
     $     ind(4)
      logical omadef
c
      iwork2 = selci_mptr(nfmax*nfmax)
      iwork3 = selci_mptr(nfmax*nfmax)
      if = selci_mptr(nnorbs)
c
      call selci_mkindb(norbs, iocc, indbar, listd, lists, ns, nd)
c
c     Now defer making f until we know if it is for single excitations
c     or diagonals (in which case we need much less).
c
      omadef = .false.
c
      numf = nf(ns)
c
      call dfill(nfmax*nroot,0.0d0,vc,1)
c
      call selci_viodiff(ioconf(1,iconf),ipack,nintpo,
     $     ncold-iconf+1,iexv)
c
      do 10 jconf = iconf,ncold
c
         iexcit = iexv(jconf-iconf+1)
         if (iexcit .gt. 4) goto 10
c
         if (iexcit.eq.0) then
            call selci_prtcon(ifllog, norbs, ioconf(1,jconf), nintpo,
     $           nbitpi)
            call selci_prtcon(ifllog, norbs, ipack, nintpo, nbitpi)
            call errquit('mkvc: iexcit=0 ... screw up',jconf,
     &       INPUT_ERR)
         else  if (iexcit.eq.2) then
            if (.not. omadef) then
               omadef = .true.
               call selci_makef(q(if), h, g, int12, int34, iocc,
     $              listd, lists, ns, nd, .false.)
            endif
            call selci_getij(i,j,ipack,ioconf(1,jconf),
     $           nintpo,nbitpi,iocc)
            call selci_makehs(work1,i,j,ns,indbar,iocc,lists,
     $           w1,w2,q(iwork2),q(iwork3),q(if),g,int12,int34,
     $           numf,numf2)
         else if(iexcit.eq.4) then
            call selci_gtijkl(ind,ipack,ioconf(1,jconf),
     $           nintpo,nbitpi,iocc,icase)
            call selci_makhdb(work1,ind,icase,ns,indbar,iocc,
     $           w1,w2,q(iwork2),q(iwork3),g,int12,int34,numf,numf2)
         else
            call errquit('strange excitation value ',iexcit, INPUT_ERR)
         endif
         jbase = indxci(jconf) + 1
c         call dgemm('n', 'n', numf, nroot, numf2, 1.0d0, work1, numf,
c     *        ci(jbase,1), nci, 1.0d0, vc, nfmax )
         call selci_mxmb(work1, 1, numf, ci(jbase,1), 1, nci,
     $        vc, 1, nfmax, numf, numf2, nroot)
 10   continue
c
      if (.not. omadef) call selci_makef(q(if), h, g, int12, int34,
     $     iocc, listd, lists, ns, nd, .true.)
      call selci_makehd(work1,.true.,
     $     ns,nd,indbar,iocc,lists,listd,w1,w2,
     $     q(iwork2),q(iwork3),q(if),h,g,int12,int34,numf)
c
      junk = selci_mfree(iwork2)
c
      end
C> @}
