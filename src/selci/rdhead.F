      subroutine selci_rdhead(title,type)
#include "implicit.fh"
#include "ciinfo.fh"
      character*(*) type
      character*80 title, dtitle(10)
      character*80 civecs, hamil, hamdg, ciconf, wmatrx
c
      call selci_getnam('civecs',civecs)
      call selci_getnam('hamil',hamil)
      call selci_getnam('hamdg',hamdg)
      call selci_getnam('ciconf',ciconf)
      wmatrx = 'wmatrix'
c
c     open these here for tidiness
c
      open(iflham,form='unformatted',status='unknown',file=hamil)
      open(iflhdg,form='unformatted',status='unknown',file=hamdg)
c
c     open and read info from the CICONF file header
c
      open(iflcon, form='unformatted',status='old',file=ciconf)
      call selci_rdhcon(iflcon, title, multi, nelec, issss, norbs,
     $     nnsmax, 
     &     nci, noconf, nintpo, nbitpi, nbpsy, isym, nsym, inttyp,
     &     nsneed)
c
c     open and read info from the MOINTS file header
c
      call selci_rdhint(iflint,ntitle,dtitle,repnuc,nnorbs,inttyp)
c
c     make, open and read info from the WMATRX file header
c 
      nsmax = nsneed(2)
      if (type.eq.'select') nsmax=nsneed(3)
      call couple(multi,nsmax)
      open(iflwmt, form='formatted',status='old',file=wmatrx)
      call selci_rdhwmt(iflwmt,mmulti, nsmax, nf, nfmax, nfmax2)
c
c     check the files for consistency
c
      if (norbs.gt.nnorbs) then
         write(ifllog,*) ' norbs = ',norbs,' nnorbs=',nnorbs
         call errquit(' info from files does not tally ',0)
      endif
c
c     if actually do PT selection need to open civecs file
c
      if (type.eq.'select') then
         open(iflciv,file=civecs,form='unformatted',status='old')
         read (iflciv) nroot,nnci
         if (nnci.ne.nci)
     $        call errquit('nci ciconf, civecs mismatch',nnci)
      endif
c
      end
