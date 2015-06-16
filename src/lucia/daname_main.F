      Subroutine DaName_Main(Lu,String,mf,wa)

      Implicit Integer (A-Z)

      Include 'fio.inc'
      Include 'blksize.fh'
      Character*(*) String
      Logical mf, wa
*      Character*8 StdNam
      Character*80 Text

      Character*16 TheName
      Data TheName/'DaName_Main'/

      If ( Query ) Call LUCIAqEnter(TheName)

      If ( Trace ) then
        Write (6,*) ' >>> Enter DaName_Main <<<'
        Write (6,*) ' unit :',Lu
        Write (6,*) ' name :',String, mf, wa
      End If

*     tmp=Lu
*     Lu=isfreeunit(tmp)
*     Check calling arguments
      If ( (Lu.le.0) .or. (Lu.gt.MxFile) )
     & STOP 'MSG: unit problem!'
*     Check for consistency
*      If ( isOpen(Lu).ne.0 )
*     & STOP 'MSG: unit used!'

*     Reformat file name to standard notation
*     (capital letters, no leading blank)
*      Call StdFmt(String,StdNam)

*     If no file name has been given link it to the default name
*      If ( StdNam.eq.'        ' )
*     &   Write(StdNam,'(A,I2.2,A)')'FT',i,'F001'

*     Open file
      temp = 0
      iRc = AixOpn(temp,String,.true.)
      If ( iRc.ne.0 ) then
*       iRc = AixErr(Text)
        STOP 'MSG: unit problem!' 
      End If
      isOpen(Lu)  = 1
      FSCB(Lu)    = temp
      LuName(Lu)    = String
      Call SetLuMark(Lu)
      Addr(Lu)    = 0
      MxAddr(Lu)  = 0
      MPUnit(0,Lu)=Lu
      Count(1,Lu) = 0
      Count(2,Lu) = 0
      Count(3,Lu) = 0
      Count(4,Lu) = 0
        Multi_File(Lu)=.False.

      MBL(Lu)=MBl_nwa

      If ( Trace ) then
        Write (6,*) ' >>> Exit DaName_Main <<<'
      End If

      If ( Query ) Call LUCIAqExit(TheName)

      Return
      End
      Subroutine SetLuMark(Lu)
      Include 'fio.inc'
      LuMark(Lu)=0
      return
      end
c $Id$
