C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCDen_Driver -- Top-level routine for coupled cluster density
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine CCDen_Driver(NOcc, NVir, T1, LDT1, g_T2, Z1, LDZ1,
     $   g_Z2, g_1PDM)
      Implicit NONE
#include "errquit.fh"
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
      Integer g_LgOmeg1, g_LgOmeg2
      Integer g_2PDMoooo
      Integer g_2PDMvvoo
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
     $   'CCDen_Driver: unable to allocate SmOmeg1', NOcc*NVir,
     &       MA_ERR)
C
      If (.NOT. MA_Push_Get(MT_Dbl, NOcc*NVir, 'CCDen SmOmeg2',
     $   H_SmOmeg2, I_SmOmeg2) ) Call ErrQuit(
     $   'CCDen_Driver: unable to allocate SmOmeg2', NOcc*NVir,
     &       MA_ERR)
C
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir, 'CCDen SmOmeg3',
     $   H_SmOmeg3, I_SmOmeg3) ) Call ErrQuit(
     $   'CCDen_Driver: unable to allocate SmOmeg3', NVir*NVir,
     &       MA_ERR)
C
      If (.NOT. MA_Push_Get(MT_Dbl, NOcc*NOcc, 'CCDen SmOmeg4',
     $   H_SmOmeg4, I_SmOmeg4) ) Call ErrQuit(
     $   'CCDen_Driver: unable to allocate SmOmeg4', NOcc*NOcc,
     &       MA_ERR)
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
C     *****************
C     * Compute 1-PDM *
C     *****************
C
      Call CCDen_1PDM( NOcc, NVir, T1, LDT1, G_T2, Z1, LDZ1, g_Z2,
     $   Dbl_MB(I_SmOmeg1), NOcc, Dbl_MB(I_SmOmeg2), NOcc,
     $   Dbl_MB(I_SmOmeg3), NVir, Dbl_MB(I_SmOmeg4), NOcc, G_1PDM )
C
C     Print this thing if requested
C
      If (Util_Print('1pdm', Print_Debug) ) Call GA_Print(G_1PDM)
C
C     -----WRITE TO DISK HERE-----
C
C     ****************************
C     * Compute 2-PDM oooo terms *
C     ****************************
C     No multi-passing -- assumes v > o and holds G(i,j,k,l) in core
C
C     To save memory, we do not save Tau = T2 + T1 * T1, or 
C     [TauZ](i,j,k,l), both of which are computed in CCDen_2PDMoooo.
C     Formation of [TauZ] costs O(o**4 v**2), and it appears also in
C     G(a,i,b,j) and G(a,i,j,k).
C
      If ( .NOT. GA_Create(MT_Dbl, NOcc*NOcc, NOcc*NOcc,
     $     'CCDen 2-PDM oooo', NOcc, NOcc, G_2PDMoooo) ) Call ErrQuit(
     $     'CCDen_Driver: unable to allocate 2-PDM oooo',
     $     NOcc*NOcc*NOcc*NOcc, GA_ERR)
C
      Call CCDen_2PDMoooo(NOcc, NVir, T1, LDT1, g_T2, g_Z2,
     $     g_1PDM, g_2PDMoooo)
C
      If (Util_Print('2pdm', Print_Debug) ) Call GA_Print(G_2PDMoooo)
C
C     -----WRITE TO DISK HERE-----
C
      If (.NOT. GA_Destroy(g_2PDMoooo) ) Call ErrQuit(
     $     'CCDen_Driver: unable to free 2-PDM oooo', 0, GA_ERR)
C
C     ****************************
C     * Compute 2-PDM vvvv terms *
C     ****************************
C
C     ********************************
C     * Create 4-index intermediates *
C     ********************************
C     These are the "large omega" quantities of APR & TJL Eq. 60--62
C     Not required for oooo and vvvv terms of 2-PDM.
C     Costs are O(o**3 v**3) for each term.
C
C     Since LgOmeg3 is easily computed as the scaled sum of LgOmeg1 and
C     LgOmeg2, we do it on the fly as needed.
C
C
      If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $     'CCDen LgOmeg1', NVir, NVir, G_LgOmeg1) ) Call ErrQuit(
     $     'CCDen_Driver: unable to allocate LgOmeg1',
     $     NVir*NOcc*NVir*NOcc, GA_ERR)
C
      If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $     'CCDen LgOmeg2', NVir, NVir, G_LgOmeg2) ) Call ErrQuit(
     $     'CCDen_Driver: unable to allocate LgOmeg2',
     $     NVir*NOcc*NVir*NOcc, GA_ERR)
C
      Call CCDen_Interm4( NOcc, NVir, G_T2, g_Z2, g_LgOmeg1, g_LgOmeg2)
C
C     Print these things if requested
C
      If (Util_Print('density 4-index intermediates',
     $     Print_Debug) ) then
         Call GA_Print(g_LgOmeg1)
         Call GA_Print(g_LgOmeg2)
      EndIf
C
C     ****************************
C     * Compute 2-PDM vvoo terms *
C     ****************************
C     No multi-passing -- assumes v > o and holds G(a,b,i,j) in core
C     NOTE: stored like T2 and Z2, despite "vvoo" labeling!
C
      If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $     'CCDen 2-PDM vvoo', NVir, NVir, G_2PDMvvoo) ) Call ErrQuit(
     $     'CCDen_Driver: unable to allocate 2-PDM vvoo',
     $     NVir*NOcc*NVir*NOcc, GA_ERR)
C
      Call CCDen_2PDMvvoo( NOcc, NVir, T1, LDT1, Z1, LDZ1, g_Z2,
     $     g_LgOmeg1, g_LgOmeg2, g_1PDM, g_2PDMvvoo)
C
C     Print these things if requested
C
      If (Util_Print('2pdm', Print_Debug) )Call GA_Print(g_2PDMvvoo)
C
C     -----WRITE TO DISK HERE-----
C
      If (.NOT. GA_Destroy(g_2PDMvvoo) ) Call ErrQuit(
     $     'CCDen_Driver: unable to free 2-PDM vvoo', 0, GA_ERR)








C
C     Get rid of temporary arrays
C
      If (.NOT. GA_Destroy(g_LgOmeg2) ) Call ErrQuit(
     $     'CCDen_Driver: unable to free LgOmeg2', 0, GA_ERR)
      If (.NOT. GA_Destroy(g_LgOmeg1) ) Call ErrQuit(
     $     'CCDen_Driver: unable to free LgOmeg1', 0, GA_ERR)
C
      If (.NOT. MA_Pop_Stack(H_SmOmeg4) ) Call ErrQuit(
     $   'CCDen_Driver: unable to free SmOmeg4', 0, MA_ERR)
      If (.NOT. MA_Pop_Stack(H_SmOmeg3) ) Call ErrQuit(
     $   'CCDen_Driver: unable to free SmOmeg3', 0, MA_ERR)
      If (.NOT. MA_Pop_Stack(H_SmOmeg2) ) Call ErrQuit(
     $   'CCDen_Driver: unable to free SmOmeg2', 0, MA_ERR)
      If (.NOT. MA_Pop_Stack(H_SmOmeg1) ) Call ErrQuit(
     $   'CCDen_Driver: unable to free SmOmeg1', 0, MA_ERR)
C
      Return
      End

