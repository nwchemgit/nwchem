      Logical function Open_3cI_File(LU, Name, Spin, NA, NI, G_Out,
     , filnam)
C$Id$
      Implicit NONE
#include "errquit.fh"
      Integer LU, Spin, NA, NI, G_Out
      Character*(*) Name
      Character*(*) FilNam! actualname
C
#include "global.fh"
#include "mafdecls.fh"
#include "ftnio.fh"
#include "stdio.fh"
#include "util.fh"
#include "eaf.fh"
C
C
C
#ifdef RIMP2_USE_DA
#include "dra.fh"
C
      Integer LNBlnk
      External LNBlnk
C
      Integer Type, AILen, FitLen
      Double Precision FileL, FileLMB
      Character*1 Sp
      Character*255  tmpname
C
C     Use the global array to figure out how large the fitting basis is
C
      Call GA_Inquire(G_Out, Type, AILen, FitLen)
C
      If ( Util_Print('integral i/o', Print_High)
     $   .AND. GA_NodeID() .eq. 0) then
         Write( LuOut, 9000) Name(:LNBlnk(Name))
C
         FileL = NI * NA * FitLen
         FileLMB = Dble( MA_SizeOf(Type, Int(FileL), MT_Byte) )
     $      / (1024 * 1024)
C
         Write( LuOut, 9010) FileL, FileLMB
      EndIf
C
 9000 Format( 1X, 'Transformed three-center integral file ', A,
     $   ' is a disk-resident array (DRA)')
 9010 Format( 1X, 'holding ', 1PE10.3, ' integrals (', F10.3, ' MB)')
C
C     Setup filename
C     This depends on ASCII character set, but is that a real limitation
C     these days?
C
      Sp = Char( Spin + IChar('a') - 1 )
C
      tmpname = name
      tmpname(lnblnk(tmpname)+1:) = sp
      call util_file_name(tmpname, .true., .true., filnam)
C
C     Try to open existing file first, then create it.  CAUTION:  This
C     may not work properly if there happens to be a file with the same
C     name, but the "wrong" contents!
C
      If ( DRA_Open(FilNam(:LNBlnk(FilNam)), DRA_RW, LU) .eq. 0) then
         Open_3ci_File = .TRUE.
      ElseIf ( DRA_Create(Type, NI*NA, FitLen,
     $   'MO 3-center integrals', FilNam(:LNBlnk(FilNam)),
     $   DRA_RW, NA, FitLen, LU) .eq. 0) then
         Open_3ci_File = .TRUE.
      Else
         Call Errquit('Open_3cI_File (DRA) failed', Spin, CALC_ERR)
      EndIf
C
#else
C
      Integer LNBlnk
      External LNBlnk
      integer ierr
      character*80 errmsg
C
      Integer RecL, Code, IALo, IAHi, SLo, SHi, NS, RecLB
      Character*1 Sp
      Character*255  String, tmpname
C
C     Use the global array to figure out how much of the s (fitting basis)
C     space we have on this node.
C
      Call GA_Distribution( G_Out, GA_NodeID(), IAlo, IAHi, Slo, SHi)
      NS = SHi - SLo + 1
C
C     Compute the record size in appropriate units for OPEN and for
C     user output, in bytes.
C
#if defined(EXT_INT)
C
C     When integers are forced to 64 bits with "-i8", record lengths
C     are unaffected, consequently the MT_DAFRecL is twice too large.
C     There seems no general solution for this.
C
      RecL = 2* MA_SizeOf( MT_Dbl, NS*NA, MT_DAFRecL)
      RecLB = MA_SizeOF(MT_DAFRecL, RecL/2, MT_Byte)
#else
      RecL = MA_SizeOf( MT_Dbl, NS*NA, MT_DAFRecL)
      RecLB = MA_SizeOF(MT_DAFRecL, RecL, MT_Byte)
#endif
C
      If ( Util_Print('integral i/o', Print_High) ) then
         Call Begin_Seq_Output
         If ( GA_NodeID() .eq. 0) then
C
C           Internal writes (to a string) must be done a line at a
C           time according to the standard.
C
            Write( String, 9000) Name(:LNBlnk(Name))
            Call Write_Seq( LuOut, String(:LNBlnk(String)) )
            Write( String, 9001) 
            Call Write_Seq( LuOut, String(:LNBlnk(String)) )
            Write( String, 9002) 
            Call Write_Seq( LuOut, String(:LNBlnk(String)) )
         EndIf
C
         If ( RecLB .gt. 0) then
            Write (String, 9010) GA_NodeID(), RecLB,
     $         Dble( NI * ReclB ) / ( 1024 * 1024 )
            Call Write_Seq(LuOut, String(:LNBlnk(String)) )
         EndIf
C
         Call End_Seq_Output
      EndIf
C
 9000 Format( 1X, 'Transformed three-center integral file ', A)
 9001 Format( 1X, ' Node', 3X, ' Rec. Len.', 3X, ' File Len.')
 9002 Format( 1X, '     ', 3X, '   (Bytes)', 3X, '  (MBytes)')
 9010 Format( 1X, I5, 3X,            I10, 3X,       F10.3)
C
C     Setup filename
C
      Sp = Char( Spin + IChar('a') - 1 )
C
* rjh
      tmpname = name
      tmpname(lnblnk(tmpname)+1:) = sp
      call util_file_name(tmpname, .true., .true., filnam)
C
#ifdef DEBUG
        write(luout,*) ga_nodeid(),' recl ',recl
        call util_flush(luout)
#endif
      If ( RecL .gt. 0) then
        ierr=eaf_open(FilNam(:LNBlnk(FilNam)), eaf_rw, LU)
#ifdef DEBUG
        write(luout,*) ga_nodeid(),' opened FD =',LU,
     ,     FilNam(:LNBlnk(FilNam))  
        call util_flush(luout)
#endif
        if (ierr .ne. 0) then
           call eaf_errmsg(ierr,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           call errquit('open_3ci: eaf_open failed', ierr, DISK_ERR)
        endif
      else
c
c     LU set to -1 to avoide close problems
c
         LU=-1
      endif
C
      Open_3ci_File = .TRUE.
      Return
C
 8000 Continue
      Write (LuErr, *) 'Open_3cI_File: Error opening file: ',
     $     FilNam(:LNBlnk(FilNam)), Code
      Call Errquit('Open_3cI_File', Code, DISK_ERR)
C
      Open_3ci_File = .FALSE.
C
#endif
C
      Return
      End
