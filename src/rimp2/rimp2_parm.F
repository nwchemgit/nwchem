      Subroutine RIMP2_Get_Param( RTDB, Thresh, MxTBatch, Mem_Opt,
     $   MxEBatch, Approx, MinEVal)
C$Id: rimp2_parm.F,v 1.3 1998-03-17 00:14:36 bernhold Exp $
      Implicit NONE
      Integer RTDB, MxTBatch, MxEBatch(2)
      Double precision Thresh, MinEVal
      Character*(*) Approx
      Character*(*) Mem_Opt
C
      Integer Type, I
      Character*(26) Date
C
#include "mafdecls.fh"
#include "rtdb.fh"
C
C     **********************
C     * Integral Threshold *
C     **********************
C
      If ( .NOT. RTDB_Get( RTDB, 'mp2:thresh',
     $  MT_Dbl, 1, Thresh) ) then
         If ( .NOT. RTDB_Get( RTDB, 'scf:thresh',
     $      MT_Dbl, 1, Thresh) ) Thresh = 1.0d-12
      EndIf
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
      Return
      End
