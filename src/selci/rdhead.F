C> \ingroup selci
C> @{
      subroutine selci_rdhead(title,type)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
#include "cselcifiles.fh"
      character*(*) type
      character*80 title, dtitle(10)
c
c     open these here for tidiness
c
      if (type .eq. 'ci') then
c     closed in makeh
         open(iflham,form='unformatted',status='unknown',file=hamil)
         open(iflhdg,form='unformatted',status='unknown',file=hamdg)
      endif
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
      call selci_rdhint(iflint,ntitle,dtitle,repnuc,nnorbs,inttyp,
     $     moints)
c
c     make, open and read info from the WMATRX file header
c 
      nsmax = nsneed(2)
      if (type.eq.'select') nsmax=nsneed(3)
      call selci_couple(multi,nsmax,0,wmatrx)
      open(iflwmt, form='formatted',status='old',file=wmatrx)
      call selci_rdhwmt(iflwmt,mmulti, nsmax, nf, nfmax, nfmax2)
c
c     check the files for consistency
c
      if (norbs.gt.nnorbs) then
         write(ifllog,*) ' norbs = ',norbs,' nnorbs=',nnorbs
         call errquit(' info from files does not tally ',0, INPUT_ERR)
      endif
c
c     if actually do PT selection need to open civecs file
c
      if (type.eq.'select') then
         open(iflciv,file=civecs,form='unformatted',status='old')
         read (iflciv) nroot,nnci
         if (nnci.ne.nci)
     $        call errquit('nci ciconf, civecs mismatch',nnci,
     &       INPUT_ERR)
      endif
c
      end
C> @}
