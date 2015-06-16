      Subroutine iDaFile_(Lu,iOpt,Buf,lBuf,iDisk)
*----------------------------------------------------------------------*
*     purpose:                                                         *
*     Control direct access I/O operations                             *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
*     iOpt    : integer, input                                         *
*               option code                                            *
*               iOpt= 1 synchonous write                               *
*               iOpt= 2 synchonous read                                *
*               iOpt= 5 synchonous rewind                              *
*     Buf     : array of integers, input/output                        *
*               Buffer carrying/receiving the data to write/read       *
*     lBuf    : integer, input                                         *
*               length of the buffer Buf                               *
*     iDisk   : integer, input/output                                  *
*               disk address                                           *
************************************************************************
      Implicit Integer (A-Z)

      Include 'blksize.fh'
      Include 'fio.inc'
      INTEGER*8 Buf(lBuf)
      Character*80 Text
      Character*16 TheName
      Data TheName/'iDaFile_'/
      Integer iU,iUinUse
      common /traceiU/iU,iUinUse

      iUinUse=0
      iU=6

      NTEST = 0
      If ( Query ) Call LUCIAqEnter(TheName)

      If (NTEST.ge.10) then
        write(6,*) ' ********************** ' 
        Write(6,*) ' >>> Enter iDaFile_ <<<'
        Write(6,*) ' unit      :',Lu
        Write(6,*) ' option    :',iOpt
        Write(6,*) ' length    :',lBuf
        Write(6,*) ' disk adr. :',iDisk
        write(6,*) ' ItoB      ;',ItoB
        write(6,*) ' Buf', Buf
        write(6,*) ' ********************** ' 
      End If

*     Check arguments
      If ( (Lu.le.0) .or. (Lu.gt.MxFile) )
     & STOP 'MSG: unit problem!'

      If ( (iOpt.lt.0) .or. (iOpt.gt.10.and.iOpt.ne.99) ) then
        write (6,*) 'Invalid action code ',iOpt
        goto 1000
      endif

      If (iOpt.eq.0.or.iOpt.eq.1.or.iOpt.eq.2.or.iOpt.eq.6.or.
     *   iOpt.eq.7.or.iOpt.eq.99) then
        If ( lBuf.lt.0 ) then
          write (6,*) 'Invalid buffer size ',lBuf
          goto 1000
        endif
      End If

      If ( iDisk.lt.0 ) then
      write (6,*) 'Invalid disk address ',iDisk
      goto 1000
      endif

      Addr(Lu) = iDisk
      lDisk    = iDisk
*     Dummy write. No I/O is made. Disk address is updated.
      If ( iOpt.eq.0 ) then
        Addr(Lu)   = Addr(Lu)+ExtAdr(BufLen(lBuf))
        MxAddr(Lu) = Max(MxAddr(Lu),Addr(Lu))
        iDisk      = Addr(Lu)
*     Write to disk.
      Else If ( (iOpt.eq.1) .or. (iOpt.eq.6) ) then
cGLM        iRc = AixWr(FSCB(Lu),Buf,lBuf,lDisk)
        If ( iRc.ne.0 ) then
*          iRc = AixErr(Text)

         write (6,*) 'Premature abort while writing buffer from disk'
         write (6,*) Text
         write (6,*) ' Unit      :',Lu
         Write (6,*) ' Option    :',iOpt
         Write (6,*) ' Buffer    :',lBuf
         write (6,*) ' Address   :',iDisk
         STOP 'Premature abort while writing buffer from disk'
*         call quit(_RC_IO_ERROR_WRITE_)
        End If
        Addr(Lu)    = Addr(Lu)+ExtAdr(BufLen(lBuf))
        MxAddr(Lu)  = Max(MxAddr(Lu),Addr(Lu))
        Count(3,Lu) = Count(3,Lu)+1
        Count(4,Lu) = Count(4,Lu)+dble(lBuf*ItoB)/1024.0d0
        iDisk       = Addr(Lu)
*     Read from disk
      Else If ( (iOpt.eq.2) .or. (iOpt.eq.7) .or. (iOpt.eq.99)) then
        if(iOpt.ne.99) then
          iRc = AixRd(FSCB(Lu),Buf,lBuf,lDisk,0)
          IF(NTEST.ge.10) then
            write(6,*) 'lDisk = ', lDisk
            write(6,*) 'in idafile_ after AixRd! Buf:'
            write(6,*) Buf
          END IF 
        endif
        If ( iRc.ne.0 ) then
*         iRc = AixErr(Text)
          write (6,*) 'Premature abort while reading buffer from disk'
          write (6,*) Text
          write (6,*) ' Unit      :',Lu
          Write (6,*) ' Option    :',iOpt
          Write (6,*) ' Buffer    :',lBuf
          write (6,*) ' Address   :',iDisk
          STOP 'Premature abort while reading buffer from disk'
*         call Quit(_RC_IO_ERROR_READ_)
        End If
        Addr(Lu)    = Addr(Lu)+ExtAdr(BufLen(lBuf))
        MxAddr(Lu)  = Max(MxAddr(Lu),Addr(Lu))
        Count(1,Lu) = Count(1,Lu)+1
        Count(2,Lu) = Count(2,Lu)+dble(lBuf*ItoB)/1024.0d0
cGLM    iDisk       = Addr(Lu)
        iDisk       = lDisk
        If (NTEST.ge.10) then
          write(6,*) ' ************************** ' 
          write(6,*) ' >>>>> in iDaFile_ <<<< '
          write(6,*) ' Lu         : ', Lu
          write(6,*) ' Addr(Lu)   : ', Addr(Lu)
          write(6,*) ' MxAddr(Lu) : ', MxAddr(Lu)
          write(6,*) ' Count(1,Lu): ', Count(1,Lu)
          write(6,*) ' Count(2,Lu): ', Count(2,Lu)
          write(6,*) ' iDisk      : ', iDisk
          write(6,*) 'in idafile_ after AixRd! Buf:'
          write(6,*) (Buf(i),i=1,lBuf)
          write(6,*) ' ************************** ' 
        END IF 
*     Gather data specified in the gather/scatter list generated
*     by a call to GSlist and write to disk
      End If
      If (NTEST.ge.10) then
        Write(6,*) ' unit      :',Lu
        Write(6,*) ' option    :',iOpt
        Write(6,*) ' length    :',lBuf
        Write(6,*) ' disk adr. :',iDisk
      End If

      If ( Trace ) Write (6,*) ' >>> Exit DaFile_ <<<'
      If ( Query ) Call LUCIAqExit(TheName)

      Return
1000  write (6,*) 'I/O error in ',TheName
      write (6,*) 'Unit = ',Lu
*      Call QTrace()
*      Call Abend()
      STOP 'I/O error in '
      End
c $Id$
