C> \ingroup selci
C> @{
C>
      subroutine selci_brdeig
*
* $Id$
*
#include "implicit.fh"
#include "ciinfo.fh"
#include "cmp2.fh"
      call ga_brdcst(20,eps,mdtob(norbs),0)
      return
      end
      subroutine selci_geteig(eigins)
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
#include "cmp2.fh"
      character*80 eigins
c$$$      logical isitthere
c
c subroutine to get orbital eigenvalues to be used in MP pt.
c currently reads ascii eigenvalues from a file called 'eigins'
c
      call errquit('select: MP denominators not functional',0,
     &       INPUT_ERR)
c$$$      if (norbs.gt.norbsmx) then
c$$$         call errquit('norbs > norbsmx = ',norbsmx)
c$$$      endif
c$$$      isitthere = .false.
c$$$      inquire(file=eigins,exist=isitthere)
c$$$      if (.not.isitthere) then
c$$$         call errquit(' file eigins missing ',911)
c$$$      endif
c$$$      open(unit=ifleig,form='formatted',access='sequential',
c$$$     &     status='old',file=eigins)
c$$$      read(ifleig,*)(eps(i),i=1,norbs)
c$$$      call dscal(norbs,0.5d00,eps,1)
c$$$c
c$$$c:debug 2
c$$$      write(6,*) ' Orbital eigenvalues '
c$$$      write(ifllog,10000) (eps(i),i=1,norbs)
c$$$10000 format(1x,6f12.6)
c$$$c
c$$$      return
      end
      double precision function selci_e0mp2(iocc)
#include "implicit.fh"
#include "ciinfo.fh"
#include "cmp2.fh"
c
      integer iocc(norbs)
c
c     return the zeroth order energy of this orbital occupation
c     using the mp2 partioning
c
      e = 0.0d0
      selci_e0mp2 = e
      do 10 i = 1, norbs
         if (iocc(i).ne.0) then
            e = e + eps(i)
            if (iocc(i).eq.3) e = e + eps(i)
         endif
 10   continue
c
      selci_e0mp2 = e
c:debug      write(6,*)' e,e0mp2,iocc = ',e,e0mp2,iocc
      end
      subroutine selci_mke0mp(ioconf, indxci, ci, erefmp, nconf)
#include "implicit.fh"
#include "ciinfo.fh"
c
      dimension ioconf(nintpo, nconf), indxci(nconf+1),
     $     ci(nci, nroot), erefmp(nroot)
      dimension iocc(255)
c
c     return in erefmp(*) the zeroth order energy of each reference
c     function with the mp2 partitioning
c
      call dfill(nroot, 0.0d0, erefmp, 1)
c
c     loop thru orbital configurations in reference space
c
c      write(6,*)' nconf = ',nconf
      do 10 iconf = 1, nconf
         call selci_upkcon(norbs, iocc, ioconf(1,iconf), nintpo, nbitpi)
         e = selci_e0mp2(iocc)
c
c     loop thru roots
c
         do 20 iroot = 1, nroot
c
c     compute norm of this orbital config for this root and
c     add in energy contribution
c
            csq = 0.0d0
c            write(6,*)' icsf(start) = ',indxci(iconf)+1
c            write(6,*)' icsf(end  ) = ',indxci(iconf+1)
            do 30 icsf = indxci(iconf)+1, indxci(iconf+1)
               csq = csq + ci(icsf, iroot)**2
 30         continue
c            write(6,*)' e = ',e,' c**2 = ',csq
            erefmp(iroot) = erefmp(iroot) + csq*e
 20      continue
 10   continue
c
      write(6,1)
 1    format(/' moller-plesset reference energies '/)
      call output(erefmp, 1, nroot, 1, 1, nroot, 1, 1)
c
      end
C>
C> @}
