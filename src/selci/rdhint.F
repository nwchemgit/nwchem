      subroutine selci_rdhint(iflint,ntitle,dtitle,repnuc,nnorbs,inttyp)
#include "implicit.fh"
c
c     common needed for gamess stuff
c
**#include "sizes.fh"
*      common /gmcrap/
*     *     pop(maxorb),potn,core,
*     *     ncolo,nbas,newb,ncore,
*     *     mapcie(maxorb),ilifc(maxorb),nval,mapaie(maxorb),
*     $     ilifa(maxorb),iqsec,jdsec
c
      character*80 dtitle(*), moints
      dimension nbpsy(8)
      real*8 repnuc
c
c     read header from integral/dump file
c
      if (inttyp.eq.1) then
c
c     gamess integrals ... dumpfile ed3 1, mainfile ed6 1
c     assume fock matrix in section 171
c
         call errquit('rdhint: gamess interface not functional',0)
*         call initb
*         call secini(1,4)
*         call secsum
*         call secget(171,1004,iblock)
*         nav = 2
*         nread = maxorb + 2 + (4*maxorb+5+nav)/nav
*         call read(pop,nread,iblock,4)
*         repnuc = core
*         ntitle = 0
*         nnorbs = nbas
      else if (inttyp.eq.0) then
c
c     ascii list of integrals
c
         call selci_getnam('moints',moints)
         open(iflint, form='formatted',status='old',file=moints)
         repnuc = 0.0d0
         ntitle = 0
         nnorbs = 255
      else if (inttyp.eq.3) then
c
c     binary version of ascii list of integrals
c
         call selci_getnam('moints',moints)
         open(iflint, form='unformatted',status='old',file=moints)
         repnuc = 0.0d0
         ntitle = 0
         nnorbs = 255
      else if (inttyp.eq.2) then
c
c     columbus integrals
c
         call selci_getnam('moints',moints)
         open(iflint, form='unformatted',status='old',file=moints)
         read(iflint) ntitle,(dtitle(i),i=1,ntitle),repnuc,nnorbs,nsym,
     &        (nbpsy(i),i=1,nsym)
      else
         call errquit('rdhead: unknown integral source ',inttyp)
      endif
c
      end
