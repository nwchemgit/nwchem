C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_NonSep2 -- Non-separable gradient contrib from 2c integrals
C
C REVISION
C     $Id: nonsep2.F,v 1.1 1999-05-28 16:03:55 bernhold Exp $
C
C SYNOPSIS
      Subroutine RIMP2_NonSep2(BraKetSame, TopSpin, NFrzO, NAct, NVir,
     $     FitBas, D_Int, D_Gam, Grad, LDGrad)
      Implicit NONE
      Logical BraKetSame        ![in]
      Integer TopSpin           ![in]
      Integer NFrzO(TopSpin)    ![in]
      Integer NAct(TopSpin)     ![in]
      Integer NVir(TopSpin)     ![in]
      Integer FitBas            ![in]
      Integer D_Int(TopSpin)    ![in]
      Integer D_Gam(TopSpin)    ![in]
      Integer LDGrad            ![in]
      Double Precision Grad(LDGrad) ![inout]
C
C DESCRIPTION
C
C     2 (ia|D) V-1(D,D') V'(D',T) G(T,ia)
C
C NOTES
C     Equations have been derived only for the 'V' approximation
C     Assumes that the ket integrals are actually (ia|D) V-1(D,D')
C
C PRINT CONTROLS
C EXTERNAL ROUTINES
c$$$      Logical  GA_Create
c$$$      External GA_Create
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "bas.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MinChunk
      Parameter (MinChunk = -64)
C
C LOCAL VARIABLES
      Integer Me
      Integer NFit, Spin
      Integer g_V
      Integer NTotal, IBlock, GBlock
      Double Precision Initial, Final, Block
C
C STATEMENT FUNCTIONS
C
C
      Me = GA_NodeID()
C
      If ( BraKetSame ) then
         If ( Me .eq. 0) Write (LuOut, 9000)
         Call ErrQuit(
     $     'RIMP2_NonSep2: invocation inconsistent with implementation',
     $        0)
      EndIf
 9000 Format(/1X, 'ERROR: RIMP2_NonSep2 implementation requires an ',
     $            'asymmetric integral',
     $       /1X, 'ERROR: approximation with the V^{-1} subsumed ',
     $            'into the ket fitting',
     $       /1X, 'ERROR: basis index (BraKetSame must be .FALSE.).')
C
      If ( .NOT. Bas_NumBF( FitBas, NFit) ) Call ErrQuit(
     $     'RIMP2_NonSep2: can''t get fitting basis size', FitBas)
C
C     *******************************************************
C     * Create target GA for 3c integral-gamma3 contraction *
C     *******************************************************
C     NOTE: Should probably distribute this on shell boundaries
C
      If ( .NOT. GA_Create(MT_Dbl, NFit, NFit, 'V-like temp',
     $     MinChunk, MinChunk, g_V) ) Call ErrQuit(
     $     'RIMP2_NonSep2: can''t allocate V-like temporary', 0)
      Call GA_Zero(g_V)
C
C     ************************************
C     * Contract 3c integrals and Gamma3 *
C     ************************************
C
      Do Spin = 1, TopSpin
         Call RIMP2_NS2_I3G3(NFrzO(Spin), NAct(Spin), NVir(Spin),
     $        NFit, D_Int(Spin), D_Gam(Spin), g_V)
      EndDo
c$$$      Call GA_Print( g_V )
C
C     **************************************
C     * Contract with derivative integrals *
C     **************************************
C
c$$$      call int_mem_2e2c(max1e, mem1)
c$$$      if (.not. MA_push_get(MT_DBL,max1e,'int_2c_ga:buf',l_buf,k_buf))
c$$$     $     call errquit('int_2c_ga: ma failed', max1e)
c$$$      if (.not. MA_push_get(MT_DBL, mem1,'int_2c_ga:scr',l_scr,k_scr))
c$$$     $     call errquit('int_2c_ga: ma failed', mem1)
C
C     Free remaining GAs
C
      If ( .NOT. GA_Destroy( g_V ) ) Call ErrQuit(
     $     'RIMP2_NonSep2: can''t destroy V-like GA', 0)
C
      Return
      End
