C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_l2 -- Compute RI-MP2 L2 Lagrangian term
C
C REVISION
C     $Id: l2.F,v 1.3 1998-10-06 17:19:16 bernhold Exp $
C
C SYNOPSIS
      Subroutine RIMP2_L2(C, A, V, NFit, D_Int_ai, d_Int_ij, D_Gam,
     $     g_L2)
      Implicit NONE
      Integer C  ![in]
      Integer A  ![in]
      Integer V  ![in]
      Integer NFit   ![in]
      Integer D_Int_ai  ![in]
      Integer D_Int_ij  ![in]
      Integer D_Gam  ![in]
      Integer g_L2   ![in]
C
C DESCRIPTION
C     Compute L^{2}_{ap} = -<ij||pb>_{RI} t_{ij}^{ab} as
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
      Integer ILo, IHi, NI, PLo, PHi, NP, P, I, CAV
      Integer JLo, JHi, NJ, J, JBlock
      Integer ReqI, ReqP
      Double Precision Initial, Final, Block, L2

C
      CAV = C + A + V
C
      Call GA_Zero(g_L2)
C
C     Determine size of a single integral/gamma block and total
C     available memory for integrals and gammas.
C
      Initial = GA_Memory_Avail()
C
      If ( .NOT. GA_Create( MT_Dbl, CAV+V, NFit,
     $     'test for memory allocation', CAV+V, -1, G_Int) )
     $     Call ErrQuit(
     $     'RIMP2_L2: can''t allocate minimal GA', (CAV+V)*NFit)
C
      Final = GA_Memory_Avail()
      Block = Final - Initial
C
      If( .NOT. GA_Destroy( G_Int) ) Call ErrQuit(
     $     'RIMP2_L2: can''t destroy GA!', 0)
C
C     Figure out a safe blocking for all nodes
C
      If ( Block .gt. 0) then
         IBlock = Min( Int( Initial / Block ), A)
      Else
         IBlock = A
      EndIf
      Call GA_IGOp(Msg_RIMP2_Mem, IBlock, 1, 'min')
C
C     Create GAs appropriate for blocking
C
      If ( .NOT. GA_Create( MT_Dbl, IBlock*CAV, NFit,
     $     '3c2e integrals', IBlock*CAV, -1, G_Int) ) Call ErrQuit(
     $     'RIMP2_L2: can''t allocate integral GA', IBlock*CAV*NFit)
C     
      If ( .NOT. GA_Create( MT_Dbl, IBlock*V, NFit,
     $     '3-index gamma', IBlock*V, -1, G_Gam) ) Call ErrQuit(
     $     'RIMP2_L2: can''t allocate gamma GA', IBlock*V*NFit)
C
C     Loop over batches of I
C
      Do ILo = 1, A, IBlock
         IHi = Min( A, ILo + IBlock - 1)
         NI = IHi - ILo + 1
C
C        Assemble Int(pi,D) from Int(ai,D) and Int(ji,D) using
C        the g_Gam GA as scratch space.
C
C        First, read Int(ai,D) into g_Gam
C
         If ( DRA_Read_Section(.FALSE.,
     $        G_Gam, 1,                  NI*V, 1, NFit,
     $        D_Int_ai, (ILo-1)*V+1, IHi*V, 1, NFit, ReqI) .ne. 0)
     $        Call ErrQuit('RIMP2_L2: Error reading integrals', ILo)
         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
     $        'RIMP2_L2: Error from DRA_Wait on integrals', 1)
C
C        Now copy blockwise into proper locations in g_Int
C
         Do I = 1, IHi - ILo + 1
            Call GA_Copy_Patch_DP('N', g_Gam, (I-1)*V+1, I*V, 1, NFit,
     $           g_Int, (I-1)*CAV+C+A+1, I*CAV, 1, NFit)
         EndDo                  ! I
         Call GA_Sync
C
C        Next, read Int(mi,D) into g_Gam.
C
C        Note: These integral blocks are of dimension IBlock*(C+A) while
C        GA is sized for IBlock*V.  If V < (C+A), we must actually do
C        this in several passes.
C
         JBlock = Int( IBlock*V / (C+A) )
         Do JLo = ILo, IHi, JBlock
            JHi = Min( IHi, JLo + JBlock - 1)
            NJ = JHi - JLo + 1
C
            If ( DRA_Read_Section(.FALSE., G_Gam, 1, NJ*(C+A), 1, NFit,
     $           D_Int_ij, (JLo-1)*(C+A)+1, JHi*(C+A), 1, NFit, ReqI)
     $           .ne. 0)
     $           Call ErrQuit('RIMP2_L2: Error reading integrals', JLo)
            If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
     $           'RIMP2_L2: Error from DRA_Wait on integrals', 1)
C
C           Now copy blockwise into proper locations in g_Int.
C
            Do J = 1, NJ
               Call GA_Copy_Patch_DP('N',
     $              g_Gam, (J-1)*(C+A)+1, J*(C+A),         1, NFit,
     $              g_Int, (J+JLo-ILo-1)*CAV+1, (J+JLo-ILo-1)*CAV+(C+A),
     $              1, NFit)
            EndDo               ! J
            Call GA_Sync
         EndDo                  ! J Blocks
C
C        Now we can finally read in the actual Gamma block
C
         If ( DRA_Read_Section(.FALSE., G_Gam, 1, NI*V, 1, NFit,
     $        D_Gam, (ILo-1)*V+1, IHi*V, 1, NFit, ReqI) .ne. 0)
     $        Call ErrQuit('RIMP2_L2: Error reading gamma', ILo)
         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
     $        'RIMP2_L2: Error from DRA_Wait on gamma', 1)
C
C     Now process all gamma and integrals in core
C
         Do I = 1, NI
C
            Call GA_MatMul_Patch('N', 'T', -2.0d0, 1.0d0,
     $           g_Gam, (I-1)*V+1, I*V, 1, NFit,
     $           g_Int, 1, NFit, (I-1)*CAV+1, I*CAV, 
     $           g_L2, 1, V, 1, CAV)
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
