C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCDen_Driver -- Top-level routine for coupled cluster density
C
C REVISION
C     $Id: ccden_driver.F,v 2.1 1997-03-24 02:01:47 gg502 Exp $
C
C SYNOPSIS
      Subroutine CCDen_Driver(NOcc, NVir, T1, LDT1, g_T2, Z1, LDZ1,
     $   g_Z2, g_1PDM)
      Implicit NONE
      Integer NOcc, NVir, LDT1, LDZ1
      Integer g_T2, g_Z2, g_1PDM
      Double Precision T1(LDT1, NVir), Z1(LDZ1, NVir)
C
C ARGUMENTS
C DESCRIPTION
C INCLUDE FILES
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
C
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer H_SmOmeg1, H_SmOmeg2, H_SmOmeg3, H_SmOmeg4
      Integer I_SmOmeg1, I_SmOmeg2, I_SmOmeg3, I_SmOmeg4
C
C
C
C     ******************
C     * Initialization *
C     ******************
C
C     ********************************
C     * Create 2-index intermediates *
C     ********************************
C     These are the "small omega" quantities of APR & TJL Eq. 52--55
C
      If (.NOT. MA_Push_Get(MT_Dbl, NOcc*NVir, 'CCDen SmOmeg1',
     $   H_SmOmeg1, I_SmOmeg1) ) Call ErrQuit(
     $   'CCDen_Driver: unable to allocate SmOmeg1', NOcc*NVir)
C
      If (.NOT. MA_Push_Get(MT_Dbl, NOcc*NVir, 'CCDen SmOmeg2',
     $   H_SmOmeg2, I_SmOmeg2) ) Call ErrQuit(
     $   'CCDen_Driver: unable to allocate SmOmeg2', NOcc*NVir)
C
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir, 'CCDen SmOmeg3',
     $   H_SmOmeg3, I_SmOmeg3) ) Call ErrQuit(
     $   'CCDen_Driver: unable to allocate SmOmeg3', NVir*NVir)
C
      If (.NOT. MA_Push_Get(MT_Dbl, NOcc*NOcc, 'CCDen SmOmeg4',
     $   H_SmOmeg4, I_SmOmeg4) ) Call ErrQuit(
     $   'CCDen_Driver: unable to allocate SmOmeg4', NOcc*NOcc)
C
      Call CCDen_Interm2( NOcc, NVir, T1, LDT1, G_T2, Z1, LDZ1, g_Z2,
     $   Dbl_MB(I_SmOmeg1), NOcc, Dbl_MB(I_SmOmeg2), NOcc,
     $   Dbl_MB(I_SmOmeg3), NVir, Dbl_MB(I_SmOmeg4), NOcc )
C
C     Print these things if requested
C
      If (Util_Print('density 2-index intermediates', Print_Debug)
     $   .AND. GA_NodeID() .eq. 0) then
         Call DGeWr(NOcc, NVir, Dbl_MB(I_SmOmeg1), NOcc, 0, 0, LuOut,
     $      'Density Intermediate (small) Omega 1', 80, 'F20.15')
         Call DGeWr(NOcc, NVir, Dbl_MB(I_SmOmeg2), NOcc, 0, 0, LuOut,
     $      'Density Intermediate (small) Omega 2', 80, 'F20.15')
         Call DGeWr(NVir, NVir, Dbl_MB(I_SmOmeg3), NVir, 0, 0, LuOut,
     $      'Density Intermediate (small) Omega 3', 80, 'F20.15')
         Call DGeWr(NOcc, NOcc, Dbl_MB(I_SmOmeg4), NOcc, 0, 0, LuOut,
     $      'Density Intermediate (small) Omega 4', 80, 'F20.15')
      EndIf
C
C     Now we can compute the 1-PDM
C
      Call CCDen_1PDM( NOcc, NVir, T1, LDT1, G_T2, Z1, LDZ1, g_Z2,
     $   Dbl_MB(I_SmOmeg1), NOcc, Dbl_MB(I_SmOmeg2), NOcc,
     $   Dbl_MB(I_SmOmeg3), NVir, Dbl_MB(I_SmOmeg4), NOcc, G_1PDM )
C
C     Print this thing if requested
C
      If (Util_Print('1pdm', Print_Debug) ) Call GA_Print(G_1PDM)
C
C     Get rid of temporary arrays
C
      If (.NOT. MA_Pop_Stack(H_SmOmeg4) ) Call ErrQuit(
     $   'CCDen_Driver: unable to free SmOmeg4', 0)
      If (.NOT. MA_Pop_Stack(H_SmOmeg3) ) Call ErrQuit(
     $   'CCDen_Driver: unable to free SmOmeg3', 0)
      If (.NOT. MA_Pop_Stack(H_SmOmeg2) ) Call ErrQuit(
     $   'CCDen_Driver: unable to free SmOmeg2', 0)
      If (.NOT. MA_Pop_Stack(H_SmOmeg1) ) Call ErrQuit(
     $   'CCDen_Driver: unable to free SmOmeg1', 0)
C
      Return
      End

