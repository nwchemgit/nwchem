*fordeck idafile $Revision: 7.7 $
      Subroutine iDaFile(Lu,iOpt,Buf,lBuf,iDisk)
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
*               If Lu={40,50,60,70,80,90} we are possibly concerned    *
*               with a multi file unit (c.f. allocdisk)                *
*     iOpt    : integer, input                                         *
*               option code                                            *
*               iOpt= 0 dummy write                                    *
*               iOpt= 99 dummy read (return buf(1)=1 in success        *
*               iOpt= 1 synchonous write                               *
*               iOpt= 2 synchonous read                                *
*               iOpt= 5 synchonous rewind                              *
*               iOpt= 6 asynchronous write                             *
*               iOpt= 7 asynchronous read                              *
*               iOpt=10 asynchronous rewind                            *
*               Note: At present the asynchronous modes are not        *
*                     supported and work identically the synchonous    *
*                     modes                                            *
*     Buf     : array of integers, input/output                        *
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
      Include 'fio.inc'
      Include 'blksize.fh'

*      Dimension Buf(lBuf)
      INTEGER*8 Buf(lBuf)

*      Include 'fnc.inc'
      NTEST = 0 
      If ( Query ) Call LUCIAqEnter('iDaFile')
      min_Block_Length=MBL(Lu)
      ItoB = ItoBsize

      If (NTEST.ge.10 ) then
        Write (6,*) ' >>> Enter iDaFile <<<'
        Write (6,*) ' unit      :',Lu
        Write (6,*) ' option    :',iOpt
        Write (6,*) ' length    :',lBuf
        Write (6,*) ' disk adr. :',iDisk
        write (6,*) ' ItoB      :', ItoB
      End If

*      If ( Multi_File(Lu) .and. MaxFileSize.ne.0 ) then
*         If (iOpt.eq.0) Then
*           Dummy write. No I/O is made. Disk address is updated.
*            Addr(Lu)   = iDisk+ExtAdr(BufLen(lBuf))
*            iDisk      = Addr(Lu)
*         Else If ( (iOpt.eq.5) .or. (iOpt.eq.10) ) then
*            Addr(Lu)   = 0
*            iDisk      = Addr(Lu)
*         Else
*            jDisk       = iDisk
*            kDisk       = iDisk
*            Call MpDaFile(Lu,MaxFileSize,iOpt,Buf,lBuf,kDisk)
*            Addr(Lu)    = jDisk+ExtAdr(BufLen(lBuf))
*            iDisk       = Addr(Lu)
*         End If
*      Else
         Call iDaFile_(Lu,iOpt,Buf,lBuf,iDisk)
*      End If

      If (NTEST.ge.10 ) Write (6,*) ' >>> Exit iDaFile <<<'
      If ( Query ) Call LUCIAqExit('iDaFile')
      Return
      End
c $Id$
