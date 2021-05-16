C> \ingroup selci
C> @{
      subroutine selci_initci(title, type)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      character*(*) type, title
c
c     initialize these common variables to avoid block data
c
      nseij = -1
      ns4eij = -1
c
c     read headers from files and generally get all the info
c     needed to get going ... basically load up ciinifo
c
      iflcon = 1
      iflint = 2
      iflwmt = 3
      ifllog = 6
      iflham = 4
      iflhdg = 7
      iflciv = 8
      ifleig = 79
c
c     open and read info from ciconf, moints, wmatrx, civecs
c
      call selci_rdhead(title,type)
c
c     make symmetry info in ciinfo ... now already on conf file
c
c      call mksym
c
c     make array for addressing lower triangles
c
      nnorbs = (norbs*(norbs+1))/2
      do 10 i = 1,norbs
         iky(i) = (i*(i-1))/2
 10   continue
c
      if (type .eq. 'ci') then
         write(ifllog,1)
 1       format(/
     & '                Configuration Driven CI Program v1.0 18/8/89'/
     & '                --------------------------------------------'//)
      else if (type.eq.'select') then
         write(ifllog,11)
 11       format(/
     & '                Perturbation Selection Program v1.0 9/28/89'/
     & '                -------------------------------------------'//)
      else
         call errquit('initci: unknown run type',0, INPUT_ERR)
      endif
      write(ifllog,2) multi, nelec, issss+1, norbs, nci, noconf,
     &     nintpo, nbitpi, nsmax, nfmax, nfmax2
      write(ifllog,3) (nf(i),i=multi-1, nsmax, 2)
      write(ifllog,4) nsym, (nbpsy(i),i=1,nsym)
      if (type.eq.'select') then
         write(ifllog,22) nroot
 22      format(' no. of roots sought     ',i7)
      endif
      write(ifllog,5) (isym(i)+1,i=1,norbs)
 2    format(' state multiplicity      ',i7/
     &       ' no. of electrons        ',i7/
     &       ' state symmetry          ',i7/
     &       ' no. of orbitals         ',i7/
     &       ' no. of configurations   ',i7/
     &       ' no. of orbital conf.    ',i7/
     &       ' no. ints per occ.       ',i7/
     &       ' no. bits per int.       ',i7/
     &       ' max no. open shells     ',i7/
     &       ' nfmax                   ',i7/
     &       ' nfmax2                  ',i7)
 3    format(' no. of spin functions     ',10i5)
 4    format(' no. of symmetry types   ',i7/
     &       ' no. of orbitals per sym   ',8i5/)
 5    format(' symmetry of the orbitals'//10(26i3/:))
      end
C> @}
