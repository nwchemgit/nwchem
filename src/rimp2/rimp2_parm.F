      Subroutine RIMP2_Get_Param( RTDB, Thresh, MxTBatch, Mem_Opt,
     $   MxEBatch, Approx, MinEVal)
C$Id$
      Implicit NONE
#include "errquit.fh"
      Integer RTDB, MxTBatch, MxEBatch(2)
      Double precision Thresh, MinEVal
      Character*(*) Approx
      Character*(*) Mem_Opt
C
      Integer Type, Junk
      Character*(26) Date
      Logical Fail
C
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
C
C     **********************
C     * Integral Threshold *
C     **********************
C
      If ( .NOT. RTDB_Get( RTDB, 'mp2:aotol2e', MT_Dbl, 1, Thresh) )
     $     Thresh = 1.0d-9      ! consistent with semi-direct MP2
C
C     ***************
C     * Batch sizes *
C     ***************
C
      If ( .NOT. RTDB_Get( RTDB, 'mp2:transformation batch size',
     $   MT_Int, 1, MxTBatch) ) MxTBatch = -1 ! Determined by memory space
C
      If ( .NOT. RTDB_CGet( RTDB, 'mp2:energy mem minimize',
     $   1, Mem_Opt) ) Mem_Opt = 'I'
C
      If ( .NOT. RTDB_Get( RTDB, 'mp2:energy batch size',
     $   MT_Int, 2, MxEBatch) ) then
         MxEBatch(1) = -1          ! Determined by memory space
         MxEBatch(2) = -1          ! Determined by memory space
      EndIf
C
C     ********************
C     * RI Approximation *
C     ********************
C
      If ( .NOT. RTDB_CGet( RTDB, 'mp2:ri approximation',
     $   1, Approx) ) Approx = 'V'
      If ( .NOT. RTDB_Get( RTDB, 'mp2:fit min eval',
     $   MT_Dbl, 1, MinEVal) ) MinEVal = 1.0d-8
C
C     ************************
C     * Check for bad inputs *
C     ************************
C     These are no longer supported, and we consider it a favor to draw
C     the attention of users to this fact, lest they get wrong results!
C
      Fail = .False.
      If ( RTDB_Get_Info( RTDB, 'mp2:freeze by atoms', Type, Junk,
     $     Date) ) then
         Fail = .True.
         Write (LuOut, '(A,A/A)')
     $        'ERROR: ''mp2:freeze by atoms'' no longer supported ',
     $        'by RI-MP2 module.',
     $        'ERROR: Use FREEZE directive instead.'
      EndIf
      If ( RTDB_Get_Info( RTDB, 'mp2:freeze orbitals', Type, Junk,
     $     Date) ) then
         Fail = .True.
         Write (LuOut, '(A,A/A)')
     $        'ERROR: ''mp2:freeze orbitals'' no longer supported ',
     $        'by RI-MP2 module.',
     $        'ERROR: Use FREEZE and VECTORS SWAP directives instead.'
      EndIf
      If ( Fail ) Call ErrQuit('rimp2: aborting', 0, UNKNOWN_ERR)
C
      Return
      End
