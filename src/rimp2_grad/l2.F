C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_l2 -- Compute RI-MP2 L2 Lagrangian term
C
C REVISION
C     $Id: l2.F,v 1.2 1998-07-16 19:59:49 bernhold Exp $
C
C SYNOPSIS
      Subroutine RIMP2_L2(NFrzO, NAct, NVir, NFit, D_Int_ai,
     $     d_Int_ij, D_Gam, g_L2)
      Implicit NONE
      Integer NFrzO  ![in]
      Integer NAct   ![in]
      Integer NVir   ![in]
      Integer NFit   ![in]
      Integer D_Int_ai  ![in]
      Integer D_Int_ij  ![in]
      Integer D_Gam  ![in]
      Integer g_L2   ![in]
C
C DESCRIPTION
C     Compute L^{2}_{ap} = -<ij||[b>_{RI} t_{ij}^{ab} as
C     L^{2}_{ap} = -2 (ip|D) V^{-1}_{DT} (T|jb) t_{ij}^{ab}
C                = -2 (pi|D) Gamma3(ai,D)
C
C NOTES
C     Current implementation is not scalable because ga_matmul_patch
C     is collective.
C
C     No attempt is made to overlap computation with disk i/o, though
C     this could easily be done using double buffering.
C
C EXTERNAL ROUTINES
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
#include "msgids.fh"
#include "util.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
C LOCAL VARIABLES
      Integer G_Int, G_Gam
      Integer NTotal, IBlock, PBlock
      Integer ILo, IHi, NI, PLo, PHi, NP, P, I, NAll
      Integer ReqI, ReqP
      Double Precision Initial, Final, Block, L2

C
      NAll = NFrzO + NAct + NVir
C
      Call GA_Zero(g_L2)
C
C     Determine size of a single integral/gamma block and total
C     available memory for integrals and gammas.
C
      Initial = GA_Memory_Avail()
C
      If ( .NOT. GA_Create( MT_Dbl, NAll+NVir, NFit,
     $     'test for memory allocation', NAll+NVir, -1, G_Int) )
     $     Call ErrQuit(
     $     'RIMP2_L2: can''t allocate minimal GA', (NAll+NVir)*NFit)
C
      Final = GA_Memory_Avail()
      Block = Final - Initial
C
      If( .NOT. GA_Destroy( G_Int) ) Call ErrQuit(
     $     'RIMP2_L2: can''t destroy GA!', 0)
C
C     If this node has no data, then it can hold the full load of
C     integrals and gammas.  We really want the minimum any node can
C     handle.
C
      If ( Block .gt. 0) then
         IBlock = Min( Int( Initial / Block ), NAct)
      Else
         IBlock = NAct
      EndIf
      Call GA_IGOp(Msg_RIMP2_Mem, IBlock, 1, 'min')
C
C     Create GAs appropriate for blocking
C
      If ( .NOT. GA_Create( MT_Dbl, IBlock*NAll, NFit,
     $     '3c2e integrals', IBlock*NAll, -1, G_Int) ) Call ErrQuit(
     $     'RIMP2_L2: can''t allocate integral GA', IBlock*NAll*NFit)
C     
      If ( .NOT. GA_Create( MT_Dbl, IBlock*NVir, NFit,
     $     '3-index gamma', IBlock*NVir, -1, G_Gam) ) Call ErrQuit(
     $     'RIMP2_L2: can''t allocate gamma GA', IBlock*NVir*NFit)
C
C     Loop over batches of I
C
      Do ILo = 1, NAct, IBlock
         IHi = Min( NAct, ILo + IBlock - 1)
         NI = IHi - ILo + 1
C
C        Assemble Int(pi,D) from Int(ai,D) and Int(ji,D) using
C        the g_Gam GA as scratch space.
C
C        First, read Int(mi,D) into g_Int as one patch
C
         If ( DRA_Read_Section(.FALSE.,
     $        G_Int, 1,                  NI*(NAct+NFrzO), 1, NFit,
     $        D_Int_ij, (ILo-1)*(NAct+NFrzO)+1, IHi*(NAct+NFrzO),
     $        1, NFit, ReqI) .ne. 0)
     $        Call ErrQuit('RIMP2_L2: Error reading integrals', ILo)
         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
     $        'RIMP2_L2: Error from DRA_Wait on integrals', 1)
C
C        Now copy blockwise into proper locations in g_Int.  Start
C        at the bottom and work up to avoid overwriting.  Data layout
C        guarantees that ga_copy_patch_dp can be used in this context.
C
         Do I = IHi - ILo + 1, 1, -1
            Call GA_Copy_Patch_DP('N',
     $           g_Int, (I-1)*(NAct+NFrzO)+1, I*(NAct+NFrzO), 1, NFit,
     $           g_Int, (I-1)*NAll+1, (I-1)*NAll+NFrzO+NAct,  1, NFit)
         EndDo                  ! I
         Call GA_Sync
C
C        Next, read Int(ai,D) into g_Gam 
C
         If ( DRA_Read_Section(.FALSE.,
     $        G_Gam, 1,                  NI*NVir, 1, NFit,
     $        D_Int_ai, (ILo-1)*NVir+1, IHi*NVir, 1, NFit, ReqI) .ne. 0)
     $        Call ErrQuit('RIMP2_L2: Error reading integrals', ILo)
         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
     $        'RIMP2_L2: Error from DRA_Wait on integrals', 1)
C
C        Now copy blockwise into proper locations in g_Int
C
         Do I = 1, IHi - ILo + 1
            Call GA_Copy_Patch_DP('N',
     $           g_Gam, (I-1)*NVir+1,            I*NVir, 1, NFit,
     $           g_Int, (I-1)*NAll+NFrzO+NAct+1, I*NAll, 1, NFit)
         EndDo                  ! I
C
C        Now we can finally read in the actual Gamma block
C
         Call GA_Zero(g_Gam)
         If ( DRA_Read_Section(.FALSE.,
     $        G_Gam, 1,               NI*NVir, 1, NFit,
     $        D_Gam, (ILo-1)*NVir+1, IHi*NVir, 1, NFit, ReqI) .ne. 0)
     $        Call ErrQuit('RIMP2_L2: Error reading gamma', ILo)
         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
     $        'RIMP2_L2: Error from DRA_Wait on gamma', 1)
C
C     Now process all gamma and integrals in core
C
         Do I = 1, IHi - ILo + 1
C
            Call GA_MatMul_Patch('N', 'T', -2.0d0, 1.0d0,
     $           g_Gam, (I-1)*NVir+1, I*NVir, 1, NFit,
     $           g_Int, 1, NFit, (I-1)*NAll+1, I*NAll, 
     $           g_L2, 1, NVir, 1, NAll)
         EndDo                  ! I
C
      EndDo                     ! I blocks
C
C     Free memory
C
      If ( .NOT. GA_Destroy( g_Gam ) ) Call ErrQuit(
     $     'RIMP2_L2: can''t destroy gamma GA', g_Gam)
      If ( .NOT. GA_Destroy( g_Int ) ) Call ErrQuit(
     $     'RIMP2_L2: can''t destroy integral GA', g_Gam)
C
C     Print results if requested
C
      If ( Util_Print('partial l', Print_Debug) ) Call GA_Print( g_L2 )
C
      Return
      End
