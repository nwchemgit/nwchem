*fordeck wr_motra_info $Revision: 7.7 $
      Subroutine WR_MOTRA_Info(Lu,iOpt,iDisk,
     &                         TCONEMO,nTCONEMO,ECOR,NSYM,
     &                         NBAS,NORB,NFRO,NDEL,MxSym,BSLBL,nBSLBL)
      Implicit Real*8 (a-h,o-z)
*#include <SysDef.fh>
      Integer*8 TCONEMO(nTCONEMO),nSym, nBas(MxSym), nOrb(MxSym),
     & nFro(MxSym), nDel(MxSym)
      Character BSLBL(nBSLBL)*1
*
      NTEST = 0
      Call iDafile(Lu,iOpt,TCONEMO,nTCONEMO,iDisk)
      if(NTEST.ge.10) then
        write(6,*) ' ***************************'
        write(6,*) ' nTCONEMO:',nTCONEMO
        write(6,*) ' TCONEMO :',TCONEMO
      End IF
      Call dDafile(Lu,iOpt,ECor,   1,       iDisk)
      if(NTEST.ge.10) then
        write(6,*) ' ECor    : ',ECor
      End IF
      Call iDafile(Lu,iOpt,nSym,   1,       iDisk)
      if(NTEST.ge.10) then
        write(6,*) ' nSym    : ',nSym
      End IF
      Call iDafile(Lu,iOpt,nBas,   MxSym,   iDisk)
      if(NTEST.ge.10) then
        write(6,*) ' nBas    : ',nBas
      End IF
      Call iDafile(Lu,iOpt,nOrb,   MxSym,   iDisk)
      if(NTEST.ge.10) then
        write(6,*) ' nOrb    : ',nOrb
      End IF
      Call iDafile(Lu,iOpt,nFro,   MxSym,   iDisk)
      if(NTEST.ge.10) then
        write(6,*) ' nFro    : ',nFro
      End IF
      Call iDafile(Lu,iOpt,nDel,   MxSym,   iDisk)
      if(NTEST.ge.10) then
        write(6,*) ' nDel    : ',nDel
      End IF
      Call cDafile(Lu,iOpt,BSLBL,  nBSLBL,  iDisk)
      if(NTEST.ge.10) then
        write(6,*) ' nBSLBL  : ', nBSLBL
        write(6,*) ' BSLBL   : '
        write(6,*) BSLBL
        write(6,*) ' ***************************'
      End IF
*
      Return
      End
c $Id$
