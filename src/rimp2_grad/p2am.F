C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_p2am -- Solve CPHF eqns for relaxed part of density
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_P2am(RTDB, TopSpin, NFrzO, NAct, NVir, g_L, g_P2)
      Implicit NONE
#include "errquit.fh"
      Integer RTDB ![in]
      Integer TopSpin ![in]
      Integer NFrzO(TopSpin) ![in]
      Integer NAct(TopSpin) ![in]
      Integer NVir(TopSpin) ![in]
      Integer g_L(TopSpin) ![in]
      Integer g_P2(TopSpin) ![in]
C
C DESCRIPTION
C PRINT CONTROLS
C EXTERNAL ROUTINES
      Logical  File_Write_GA, File_Read_GA, CPHF
      External File_Write_GA, File_Read_GA, CPHF
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
C LOCAL VARIABLES
      Integer ISpin, Length, C, A, V
      Integer g_RHS
      Character*(NW_Max_Path_Len) FN_CPHF_RHS, FN_CPHF_Sol
C
C STATEMENT FUNCTIONS
C

C
C     Create array to hold RHS
C
      Length = 0
      Do ISpin = 1, TopSpin
         Length = Length + (NVir(ISpin) * (NFrzO(ISpin)+NAct(ISpin) ) )
      EndDo
C
*ga:1:-1
      If ( .NOT. GA_Create(MT_Dbl, Length, 1, 'RI-MP2 CPHF RHS/Soln', 
     $     -1, -1, g_RHS ) ) Call ErrQuit(
     $     'RIMP2_P2am: can''t allocate RHS', 0, GA_ERR)
C
C     Copy appropriate bits of L into RHS
C
      Length = 0
      Do ISpin = 1, TopSpin
         C = NFrzO(ISpin)
         A = NAct(ISpin)
         V = NVir(ISpin)
C
         Call GA_Copy_Patch('N', g_L(ISpin), C+A+1, C+A+V, 1, C+A,
     $        g_RHS, 1+Length, (V*(C+A))+Length, 1, 1)
C
C        Need length of alpha spin for next spin case
C
         Length = Length + (V * (C+A) )
      EndDo
C
      If ( Util_Print('partial p2', Print_Debug) .OR.
     $     Util_Print('cphf', Print_Debug) ) Call GA_Print(g_RHS)
C
C     CPHF solver expects RHS on a file.  Name is fixed by CPHF module.
C
      Call Util_File_name('cphf_rhs', .TRUE., .TRUE., FN_CPHF_RHS)
      If (.NOT. File_Write_GA(FN_CPHF_RHS, g_RHS) ) Call ErrQuit(
     $     'RIMP2_P2am: can''t write CPHF RHS to file', 0, GA_ERR)
C
C     Invoke the CPHF solver
C
      If ( .NOT. CPHF( RTDB ) ) Call ErrQuit(
     $     'RIMP2_P2am: cphf solution failed', 0, CALC_ERR)
C
C     Once again, solution is returned in a file, name fixed by CPHF module
C
      Call Util_File_Name('cphf_sol', .TRUE., .TRUE., FN_CPHF_Sol)
      If ( .NOT. File_Read_GA(FN_CPHF_Sol, g_rhs) ) Call ErrQuit(
     $     'RIMP2_P2am: can''t read CPHF sol''n file', 0, DISK_ERR)
C
      If ( Util_Print('partial p2', Print_Debug) .OR.
     $     Util_Print('cphf', Print_Debug) ) Call GA_Print(g_RHS)
C
C     Copy appropriate bits of CPHF solution into P2
C
      Call GA_Scale(g_RHS, 2.0d0) ! CHEATING!
      Length = 0
      Do ISpin = 1, TopSpin
         C = NFrzO(ISpin)
         A = NAct(ISpin)
         V = NVir(ISpin)
C
         Call GA_Copy_Patch('N', g_RHS, 1+Length, (V*(C+A))+Length,
     $        1, 1, g_P2(ISpin), C+A+1, C+A+V, 1, C+A)
         Call GA_Copy_Patch('T', g_P2(ISpin), C+A+1, C+A+V, 1, C+A, 
     $        g_P2(ISpin), 1, C+A, C+A+1, C+A+V )
C
C        Need length of alpha spin for next spin case
C
         Length = Length + (V * (C+A) )
C
         If ( Util_Print('partial p2', Print_Debug) ) then
            Write (LuOut, *) 'After RIMP2_P2am...'
            Call GA_Print(g_P2(ISpin) )
         EndIf
      EndDo
C
C     THIS IS A TEMPORARY HACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      Write (LuOut, *)
     $     'rimp2_p2am: WARNING: copying P(2) alpha to P(2) beta!!!!'
      Call GA_Copy(g_P2(1), g_P2(2) ) ! CHEATING!!!
C
C     Free temporary memory
C
      If ( .NOT. GA_Destroy( g_RHS) ) Call ErrQuit(
     $     'RIMP2_P2am: can''t destroy CPHF RHS GA', 0, GA_ERR)
C
C     Clean up files
C
      If ( GA_NodeID() .eq. 0) Call Util_File_Unlink(FN_CPHF_RHS)
      If ( GA_NodeID() .eq. 0) Call Util_File_Unlink(FN_CPHF_sol)
C
      Return
      End
