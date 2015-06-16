*fordeck cdafile $Revision: 7.7 $
      Subroutine cDaFile(Lu,iOpt,Buf,lBuf_,iDisk)
************************************************************************
*                                                                      *
*     (c) Copyright. All rights reserved                               *
*                                                                      *
*     This code or parts thereof may not be copdied or redistributed,  *
*     without the written permission of the author. The use is re -    *
*     stricted to research purposes only, and the matrial may not be   *
*     included in any commercial product.                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     purpose:                                                         *
*     Control direct access I/O operations                             *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
*     iOpt    : integer, input                                         *
*               option code                                            *
*               iOpt= 0 dummy write                                    *
*               iOpt= 1 synchonous write                               *
*               iOpt= 2 synchonous read                                *
*               iOpt= 5 synchonous rewind                              *
*               iOpt= 6 asynchronous write                             *
*               iOpt= 7 asynchronous read                              *
*               iOpt=10 asynchronous rewind                            *
*               Note: At present the asynchronous modes are not        *
*                     supported and work identically the synchonous    *
*                     modes                                            *
*     Buf     : array of characters, input/output                      *
*               Buffer carrying/receiving the data to write/read       *
*     lBuf    : integer, input                                         *
*               length of the buffer Buf                               *
*     iDisk   : integer, input/output                                  *
*               disk address                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
*     L Serrano-Andres,University of Lund, Sweden, 1996                *
*     R. Lindh, Univeristy of Lund, Sweden, 2001                       *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

*#include <SysDef.fh>
      Include 'blksize.fh'
      Include 'fio.inc'

      Character*1 Buf(lBuf_)
*      Include 'fnc.inc'

      If ( Query ) Call LUCIAqEnter('cDaFile')
      min_Block_Length=MBL(Lu)
      lBuf=lBuf_/ItoB
      If ( Trace ) then
        Write (6,*) ' >>> Enter cDaFile <<<'
        Write (6,*) ' unit      :',Lu
        Write (6,*) ' option    :',iOpt
        Write (6,*) ' length    :',lBuf
        Write (6,*) ' disk adr. :',iDisk
      End If

      Call iDaFile_(Lu,iOpt,Buf,lBuf,iDisk)

      If ( Trace ) Write (6,*) ' >>> Exit cDaFile <<<'
      If ( Query ) Call LUCIAqExit('cDaFile')
      Return
      End
c $Id$
