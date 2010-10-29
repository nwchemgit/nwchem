C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_prep_sep -- Prepare for computation of separable gradient terms
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_Prep_Sep(TopSpin, g_P2, g_W2, g_AOXF)
      Implicit NONE
#include "errquit.fh"
C
      Integer TopSpin ![in]
      Integer g_P2(TopSpin) ![in]
      Integer g_W2(TopSpin) ![in]
      Integer g_AOXF(TopSpin) ![in]
C
C DESCRIPTION
C PRINT CONTROLS
C EXTERNAL ROUTINES
      Logical  File_Write_GA
      External File_Write_GA
C
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
     $     'RIMP2_Prep_Sep: MyMaxSpin inconsistent with TopSpin',
     $     MyMaxSpin , INPUT_ERR)
C
C     **************************************************
C     * Transform densities to AO basis & dump to disk *
C     **************************************************
C
      Do Spin = 1, TopSpin
C
C        Need transposed MO vectors for this "back transformation"
C        Inside loop because AOXF may change size with spin
C
         Call GA_Inquire(g_AOXF(Spin), Type, NAO, NMO)
         If ( .NOT. GA_Create(Type, NMO, NAO, 'Transposed MO vectors',
     $        0, 0, g_AOXF_tr) ) Call ErrQuit(
     $        'RIMP2_Prep_Sep: can''t allocate transposed MO vectors',
     $        NAO*NMO, GA_ERR)
         Call GA_Transpose( g_AOXF(Spin), g_AOXF_tr )
C
C        Make scratch array
C
         If ( .NOT. GA_Create(MT_Dbl, NMO, NAO, 'transf scr',
     $        0, 0, g_Tmp) ) Call ErrQuit(
     $        'RIMP2_Prep_Sep: can''t allocate transformation scratch',
     $        NMO*NAO, GA_ERR)
C
C        Make temporary AO density array
C
         If ( .NOT. GA_Create(MT_Dbl, NAO, NAO, 'AO basis P(2) or W(2)',
     $        0, 0, g_ao) ) Call ErrQuit(
     $        'RIMP2_Prep_Sep: can''t allocate ao basis W(2)', NAO*NAO,
     &       GA_ERR)
C
C        Transform W(2)
C
         Call Two_Index_Transf( g_W2(Spin), g_AOXF_tr, g_AOXF_tr,
     $        g_Tmp, g_ao)
         Call GA_Symmetrize( g_ao )
C
         If ( Util_Print('ao w2', Print_Debug) ) then
            Write (LuOut, *) 'AO basis W(2)...'
            Call GA_Print( g_ao )
         EndIf
C
C        Write W(2)
C
         Call Util_File_Name(Def_FN_W2(Spin), .TRUE., .TRUE., FN_W2_ao)
         If ( .NOT. File_Write_GA(FN_W2_ao, g_ao) ) Call ErrQuit(
     $        'RIMP2_Prep_Sep: Failed to write AO W2 to disk', Spin,
     &       GA_ERR)
C
C        Transform P(2)
C
         Call Two_Index_Transf( g_P2(Spin), g_AOXF_tr, g_AOXF_tr,
     $        g_Tmp, g_ao)
         Call GA_Symmetrize( g_ao )
C
         If ( Util_Print('ao p2', Print_Debug) ) then
            Write (LuOut, *) 'AO basis P(2)...'
            Call GA_Print( g_ao )
         EndIf
C
C        Write P(2)
C
         Call Util_File_Name(Def_FN_P2(Spin), .TRUE., .TRUE., FN_P2_ao)
         If ( .NOT. File_Write_GA(FN_P2_ao, g_ao) ) Call ErrQuit(
     $        'RIMP2_Prep_Sep: Failed to write AO P2 to disk', Spin,
     &       DISK_ERR)
C
C        Clean up temporaries
C
         If ( .NOT. GA_Destroy( g_ao ) ) Call ErrQuit(
     $        'RIMP2_Prep_Sep: can''t destroy ao basis density',
     $        Spin, GA_ERR)
         If ( .NOT. GA_Destroy( g_Tmp ) ) Call ErrQuit(
     $        'RIMP2_Prep_Sep: can''t destroy temporary GA',
     $        Spin, GA_ERR)
         If ( .NOT. GA_Destroy( g_AOXF_tr ) ) Call ErrQuit(
     $        'RIMP2_Prep_Sep: can''t destroy transposed MO vectors',
     $        Spin, GA_ERR)
      EndDo
C
      Return
      End
