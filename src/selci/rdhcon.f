      subroutine selci_rdhcon(iflcon, title, multi, nelec, issss, norbs,
     &     nnsmax, nci, noconf, nintpo, nbitpi, nbpsy, isym, nsym,
     &     inttyp, nsneed)
*
* $Id: rdhcon.f,v 1.2 1997-10-31 23:42:23 d3e129 Exp $
*
      character*80 title
      dimension nbpsy(8), isym(255), nsneed(3)
c
c     read header of the ciconf file
c
      read(iflcon) title, multi, nelec, issss, norbs, nnsmax, 
     &     nci, noconf, nintpo, nbitpi, nbpsy, isym, nsym, inttyp,
     &     nsneed
c
      end
