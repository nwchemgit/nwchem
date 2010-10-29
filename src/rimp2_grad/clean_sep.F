C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_clean_sep -- Clean up after separable gradient terms
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_Clean_Sep(TopSpin)
      Implicit NONE
#include "errquit.fh"
C
      Integer TopSpin ![in]
C
C DESCRIPTION
C PRINT CONTROLS
C EXTERNAL ROUTINES
C INCLUDE FILES
#include "stdio.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "global.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MyMaxSpin
      Parameter (MyMaxSpin = 2)
C
      Character*(6) Def_FN_P2(MyMaxSpin), Def_FN_W2(MyMaxSpin)
      Data Def_fn_p2 /'p_a_ao', 'p_b_ao'/
      Data Def_fn_w2 /'w_a_ao', 'w_b_ao'/
C
C LOCAL VARIABLES
      Character*(NW_Max_Path_Len) FN_P2_ao, FN_W2_ao
      Integer Spin, Type, NAO, NMO
      Integer g_Tmp, g_AOXF_tr, g_ao
C
C STATEMENT FUNCTIONS

C
      If ( TopSpin .gt. MyMaxSpin ) Call ErrQuit( 
     $     'RIMP2_Clean_Sep: MyMaxSpin inconsistent with TopSpin',
     $     MyMaxSpin , INPUT_ERR)
C
      Do Spin = 1, TopSpin
         Call Util_File_Name(Def_FN_W2(Spin), .TRUE., .TRUE., FN_W2_ao)
         If ( GA_NodeID() .eq. 0) Call Util_File_Unlink(FN_W2_ao)
         Call Util_File_Name(Def_FN_P2(Spin), .TRUE., .TRUE., FN_P2_ao)
         If ( GA_NodeID() .eq. 0) Call Util_File_Unlink(FN_P2_ao)
      EndDo
C
      Return
      End
