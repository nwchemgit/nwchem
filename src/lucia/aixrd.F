*fordeck aixrd $Revision: 7.7 Patch(7.7): 005_unfixme $
************************************************************************
*                                                                      *
*                             A I X - I / O                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* rc=AixRd(Handle,Buf,nBuf,iDisk,iErrSkip)                             *
*                                                                      *
* A buffer is read from a file associated with the file handle. The    *
* operation is asynchronous, and must be followed by a call to AixWt   *
* to ensure that data is in memory.                                    *
*                                                                      *
* Input:  Handle   - This is the unique file identifier associated     *
*                    with the file. It is created by AixOpn, and must  *
*                    be used on subsequent references to the file.     *
*         Buf      - The buffer that is to be written to disk.         *
*         nBuf     - Length of the buffer in words.                    *
*         iDisk    - External disk address.                            *
*         iErrSkip - if 0 : stop in the error case                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          S&TC, ACIS, IBM Sweden                                      *
* Written: November 1990                                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*                                                                      *
************************************************************************
      Integer Function AixRd(handle,Buf,nBuf,iDisk,iErrSkip)
      Include 'ioparm.inc'
*     Include 'SysDef.inc'
*#include <SysDef.fh>
*      Include 'blksize.inc'
      Include 'blksize.fh'
      Include 'switch.inc'
      Include 'ctl.inc'
*      Dimension Buf(*)
      INTEGER*8 Buf(nBuf)
      Character*80 ErrTxt
      Character*16 TheName
*      Include 'fnc.inc'
*#include "warnings.fh"
      Data TheName/'AixRd'/
      NTEST = 0
*----------------------------------------------------------------------*
* Entry to AixRd                                                       *
*----------------------------------------------------------------------*
      AixRd=0
*----------------------------------------------------------------------*
* Check if file is opened.                                             *
*----------------------------------------------------------------------*
      IF (NTEST.ge.10) then
      write(6,*) ' *********************** '
      write(6,*) ' Entering AixRd function '
      write(6,*) ' INPUT Parameters:       '
      write(6,*) ' 1) handle = ', handle
      write(6,*) ' 2) nBuf   = ', nBuf
      write(6,*) ' 2) ItoB   = ', ItoB
      END IF
      n=1
100   If(CtlBlk(pHndle,n).ne.handle) Then
         n=n+1
         If(n.gt.MxFile) Then
            AixRd=eNtOpn
            Return
         End If
         Go To 100
      End If
      nFile=n
      IF (NTEST.ge.10) then
       write(6,*) ' 3) nFile  = ', nFile
      END IF
      desc=CtlBlk(pDesc,nFile)
      IF (NTEST.ge.10) then
        write(6,*) ' 4) desc   = ', desc
      END IF
*----------------------------------------------------------------------*
* Position file pointer                                                *
*----------------------------------------------------------------------*
      pDisk=IntAdr(iDisk)
      IF (NTEST.ge.10) then
        write(6,*) ' 5) pDisk  = ', pDisk
      END IF
*
      If(CtlBlk(pWhere,nFile).ne.pDisk) Then
         rc=c_lseek(desc,pDisk)
         If(rc.lt.0) Then
            if(iErrSkip.eq.1) then
             AixRd=99
             return
            endif
*            idum=iPrintLevel(3)
*            Call FASTIO('STATUS')
*            AixRd=AixErr(ErrTxt)
*            Call SysWarnFileMsg(TheName, FCtlBlk(nFile),
*     *                                'MSG: seek', ErrTxt)
            STOP 'MSG: seek'
*           Call SysCondMsg('rc < 0', rc, '<', 0)
         Else If(rc.ne.pDisk) Then
            if(iErrSkip.eq.1) then
             AixRd=99
             return
            endif
*            idum=iPrintLevel(3)
*            Call FASTIO('STATUS')
*            AixRd=eInErr
*            Call SysWarnFileMsg(TheName, FCtlBlk(nFile),
*     *                                'MSG: seek', ' ')
             STOP 'MSG: seek'
*           Call SysCondMsg('rc != pDisk', rc, '!=', pDisk)
         End If
      End If
      CtlBlk(pWhere,nFile)=pDisk
*----------------------------------------------------------------------*
* Read from file                                                       *
*----------------------------------------------------------------------*
      LenBuf=BufLen(nBuf)
      IF (NTEST.ge.10) then
        write(6,*) ' 7) LenBuf = ', LenBuf
        write(6,*) ' *********************** '
      END IF
      rc=c_read(desc,Buf,LenBuf)

      If(rc.lt.0) Then
            if(iErrSkip.eq.1) then
             AixRd=99
             return
            endif
*         idum=iPrintLevel(3)
*         Call FASTIO('STATUS')
*         AixRd=AixErr(ErrTxt)
*            Call SysQuitFileMsg(_RC_IO_ERROR_READ_,
*     *                                TheName,FCtlBlk(nFile),
*     *      'Premature abort while reading buffer from disk', ErrTxt)
         STOP 'Premature abort while reading buffer from disk!'
      Else If(rc.ne.LenBuf) Then
            if(iErrSkip.eq.1) then
             AixRd=99
             return
            endif
*         idum=iPrintLevel(3)
*         Call FASTIO('STATUS')
*         AixRd=eEof
*            Call SysQuitFileMsg(_RC_IO_ERROR_READ_,
*     *            TheName,FCtlBlk(nFile),
*     *            'Premature abort while reading buffer from disk:',
*     *      '\n End of file reached ')
         STOP 'Premature abort while reading buffer from disk!!'
      End If
      CtlBlk(pWhere,nFile)=CtlBlk(pWhere,nFile)+LenBuf
      iDisk=iDisk+ExtAdr(LenBuf)
*----------------------------------------------------------------------*
* Finished so return to caller                                         *
*----------------------------------------------------------------------*
      Return
      End

      Function ExtAdr(t)
      INTEGER t, ExtAdr
      Include 'blksize.fh'
*
      ExtAdr=1+(t-1)/min_Block_Length
      Return
      End 

      Function IntAdr(t)
      INTEGER t
      Include 'blksize.fh'
*
      IntAdr=min_Block_Length*t
      Return
      End 

      Function BufLen(t)
      INTEGER t, BufLen
      Include 'blksize.fh'
*
      BufLen=ItoB*t
      Return
      End 



c $Id$
