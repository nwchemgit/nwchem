      Logical function Open_3cI_File(LU, Name, Spin, NA, NI, G_Out)
C$Id: open_3ci.F,v 1.1 1995-12-13 01:36:52 d3g681 Exp $
      Implicit NONE
      Integer LU, Spin, NA, NI, G_Out
      Character*(*) Name
C
#include "global.fh"
#include "mafdecls.fh"
#include "ftnio.fh"
#include "stdio.fh"
#include "util.fh"
C
C
C
#ifdef RIMP2_USE_DA
#include "da.fh"
C
      Integer LNBlnk
      External LNBlnk
C
      Integer Type, AILen, FitLen
      Double Precision FileL, FileLMB
      Character*1 Sp
      Character*256 FilNam
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
     $   ' is a disk array (DA)')
 9010 Format( 1X, 'holding ', 1PE10.3, ' integrals (', F10.3, ' MB)')
C
C     Setup filename
C     This depends on ASCII character set, but is that a real limitation
C     these days?
C
      Sp = Char( Spin + IChar('a') - 1 )
C
      FilNam = Name(:LNBlnk(Name))
      FilNam(LNBlnk(FilNam)+1:) = Sp
C
C     By opening with DA_RW we avoid having to worry about whether the
C     file already exists or not.
C
      If ( DA_Create(Type, NI*NA, FitLen,
     $   'MO 3-center integrals', FilNam(:LNBlnk(FilNam)),
     $   DA_RW, NA, FitLen, LU) .ne. 0)
     $   Call Errquit('Open_3cI_File (DA) failed', Spin)
C
      Open_3ci_File = .TRUE.
C
#else
C
      Integer LNBlnk
      External LNBlnk
C
      Integer RecL, Code, IALo, IAHi, SLo, SHi, NS, Proc, RecLB
      Character*4 PID
      Character*1 Sp
      Character*256 FilNam, String
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
#if defined(SGITFP) || defined(DECOSF)
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
C     Unfortunately, we cant force leading zeros on output with the
C     BZ format spec, so we have to account for each case.
C
      Proc = GA_NodeID()
      If ( Proc .le. 9) then
         Write(PID, '(''000'',I1)') Proc
      ElseIf ( Proc .le. 99) then
         Write(PID, '(''00'',I2)') Proc
      ElseIf ( Proc .le. 999) then
         Write(PID, '(''0'', I3)') Proc
      ElseIf ( Proc .le. 9999) then
         Write(PID, '(I4)') Proc
      Else
         Call ErrQuit('Open_3cI_File: Can''t create file name for node',
     $      Proc)
      EndIf
C
C     This depends on ASCII character set, but is that a real limitation
C     these days?
C
      Sp = Char( Spin + IChar('a') - 1 )
C
      FilNam = Name(:LNBlnk(Name))
      FilNam(LNBlnk(FilNam)+1:) = Sp//PID
C
      If ( RecL .gt. 0)
     $   Open(UNIT=LU, FILE=FilNam(:LNBlnk(FilNam)),
     $   ACCESS='DIRECT', RECL=RecL, ERR=8000, IOSTAT=Code)
C
      Open_3ci_File = .TRUE.
      Return
C
 8000 Continue
      Write (LuErr, *) 'Open_3cI_File: Error opening file: ',
     $     FilNam(:LNBlnk(FilNam)), Code
      Call Errquit('Open_3cI_File', Code)
C
      Open_3ci_File = .FALSE.
C
#endif
C
      Return
      End
