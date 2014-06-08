*----------------------------------------------------------------------*
* rc=AixOpn(Handle,FileName)                                           *
*                                                                      *
* A file is opened for read/write operations. If the file does not     *
* exist it is automatically created.                                   *
*                                                                      *
* Input:  FileName - A character string specifying a complete path     *
*                    name or relative path name for the file. The      *
*                    name must be shorter than 128 characters.         *
*                                                                      *
* Output: Handle   - When a file is sucessfully opened, a unique file  *
*                    id is supplied by the routine. This is used for   *
*                    specifying the file to other routines.            *
*----------------------------------------------------------------------*
* History:                                                             *
* 911021 - If return code is 13 when opening file, try to open it as a *
*          read only file. Per-Olof Widmark.                           *
* 931217 - Flags obtained from system calls.                           *
* 010710 - Change using iRand to incremental number                    *
*                                                                      *
************************************************************************
      Integer Function AixOpn(handle,name,translate)
      Include 'ioparm.inc'
      Include 'switch.inc'
      Include 'ctl.inc'
      Character*(*) name
      Character*256 tmp
      Logical Translate
c      Integer Length,Prog_Length
      Integer StrnLn
      External StrnLn
      Character*256 tmp1
      Character*80 ErrTxt
      External Get_Progname
      Character*100 Get_Progname
      save NVV
      data NVV /666/
      NTEST = 0
*----------------------------------------------------------------------*
* Entry to AixOpn                                                      *
*----------------------------------------------------------------------*
      AixOpn=0
*----------------------------------------------------------------------*
* Check if slot in table is available                                  *
*----------------------------------------------------------------------*
      n=1
100   If(CtlBlk(pStat,n).ne.vNaF) Then
         n=n+1
         If(n.gt.MxFile) Then
            AixOpn=eTmF
*            Call SysWarnMsg('Aixopn','Too many opened files\n',
*     *      'try to increase MxFile')
            STOP 'Aixopn: Too many opened files. Try to increase MxFile'
            Return
         End If
         Go To 100
      End If
      nFile=n
*----------------------------------------------------------------------*
* Strip file name and append string terminator                         *
*----------------------------------------------------------------------*
      n=Len(name)
200   If(name(n:n).eq.' ') Then
         n=n-1
         If(n.le.0) Then
            AixOpn=eBlNme
            Return
         End If
         Go To 200
      End If
      n=n+1
      If(n.ge.Len(tmp)) Then
         AixOpn=eTlFn
         Return
      End If
      tmp=name
      tmp(n:n)=Char(0)
*----------------------------------------------------------------------*
* Attempt to open file.                                                *
*----------------------------------------------------------------------*
      rc = 0
      tmp1=tmp
      ltmp=StrnLn(tmp1)
*      if(translate) then
cGLM        call PrgmTranslate(tmp1,tmp,ltmp)
*      endif
       tmp=tmp(1:ltmp)
       tmp(ltmp+1:ltmp+1)=Char(0)
       rc=c_open(tmp)
      If(rc.lt.0) Then
*         rc=AixErr(ErrTxt)
*        Call SysWarnFileMsg('AixOpn',name,
*     *            'MSG: open',ErrTxt)
        STOP 'MSG: open'
*       call SysPutsEnd()
*       Call Abend()
      End If
      desc=rc
*----------------------------------------------------------------------*
* Attempt sucessful, update control blocks.                            *
*----------------------------------------------------------------------*
c      handle=iRand()
       NVV=NVV+100
       handle=NVV
      CtlBlk(pHndle,nFile)=handle
      CtlBlk(pDesc ,nFile)=desc
      CtlBlk(pStat ,nFile)=vOpn
      CtlBlk(pWhere,nFile)=vNull
      FCtlBlk(nFile)=name
      IF (NTEST.ge.10) then
        write(6,*) ' *************************************************'
        write(6,*) ' In AixOpn:                   '
        write(6,*) ' nFile: ', nFile
        write(6,*) ' pHndle, CtlBlk(pHndle,nFile) : ', pHndle, handle
        write(6,*) ' pDesc , CtlBlk(pDesc ,nFile) : ', pDesc, desc
        write(6,*) ' pStat , CtlBlk(pStat ,nFile) : ', pStat, vOpn
        write(6,*) ' pWhere, CtlBlk(pWhere,nFile) : ', pWhere, vNull
        write(6,*) ' **************************************************'
      End if
*----------------------------------------------------------------------*
* Finished so return to caller                                         *
*----------------------------------------------------------------------*
      Return
      End
c $Id$
