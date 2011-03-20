C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_l1 -- Compute RI-MP2 L1 Lagrangian term
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_L1(NFrzO, NAct, NVir, NFit, D_Int, D_Gam,
     $     g_L1)
      Implicit NONE
#include "errquit.fh"
      Integer NFrzO ![in]
      Integer NAct  ![in]
      Integer NVir  ![in]
      Integer NFit  ![in]
      Integer D_Int ![in]
      Integer D_Gam ![in]
      Integer g_L1  ![in]
C
C DESCRIPTION
C     Compute L^{1}_{pi} = <pj||ab>_{RI} t_{ij}^{ab} as
C     L^{1}_{pi} = 2 (pa|D) V^{-1}_{DT} (T|jb) t_{ij}^{ab}
C                = 2 Int(ap|D) Gam(ai,D)
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
      Integer ILo, IHi, NI, PLo, PHi, NP, P, I
      Integer ReqI, ReqP
      Double Precision Initial, Final, Block, L1
C
      Call GA_Zero( g_L1 )
C
C     Determine size of a single integral/gamma block and total
C     available memory for integrals and gammas.
C
      Initial = GA_Memory_Avail()
C
      If ( .NOT. GA_Create( MT_Dbl, NVir, NFit,
     $     'test for memory allocation', NVir, -1, G_Int) )
     $     Call ErrQuit(
     $     'RIMP2_L1: can''t allocate minimal GA', NVir*NFit, GA_ERR)
C
      Final = GA_Memory_Avail()
      Block = Initial - Final
C
      If( .NOT. GA_Destroy( G_Int) ) Call ErrQuit(
     $     'RIMP2_L1: can''t destroy GA!', 0, GA_ERR)
C
C     If this node has no data, then it can hold the full load of
C     integrals and gammas.  We really want the minimum any node can
C     handle.
C
      If ( Block .gt. 0) then
         NTotal = Int( Initial / Block )
      Else
         NTotal = NFrzO+NAct+NVir + NAct
      EndIf
      Call GA_IGOp(Msg_RIMP2_Mem, NTotal, 1, 'min')
c$$$      Write ( 6, *) 'Initial, Final, Block, NTotal = ',
c$$$     $     Initial, Final, Block, NTotal
C
C     Allocate the total number of blocks we have space for so as to
C     minimize the total number of passes through the disk files.
C
      Do IBlock = NAct, 1, -1
         PBlock = NTotal-IBlock
c$$$         Write ( 6, *) 'Testing IBlock, PBlock = ', IBlock, PBlock
         If ( PBlock .gt. 0 ) Goto 1000
      EndDo
      Call ErrQuit('RIMP2_L1: insufficient memory for memory buffers',
     $     NTotal, MEM_ERR)
 1000 Continue
      PBlock = Min( PBlock, NFrzO+NAct+NVir)
c$$$      Write ( 6, *) 'Chosen PBlock = ', PBlock
C
C     Create GAs appropriate for blocking
C
      If ( .NOT. GA_Create( MT_Dbl, PBlock*NVir, NFit,
     $     '3c2e integrals', PBlock*NVir, -1, G_Int) ) Call ErrQuit(
     $     'RIMP2_L1: can''t allocate integral GA', PBlock*NVir*NFit,
     &       GA_ERR)
C     
      If ( .NOT. GA_Create( MT_Dbl, IBlock*NVir, NFit,
     $     '3-index gamma', IBlock*NVir, -1, G_Gam) ) Call ErrQuit(
     $     'RIMP2_L1: can''t allocate gamma GA', IBlock*NVir*NFit,
     &       GA_ERR)
C
C     Loop over batches of P and I
C
      Do ILo = 1, NAct, IBlock
         IHi = Min( NAct, ILo + IBlock - 1)
         NI = IHi - ILo + 1
C
c$$$         Write (6, *) 'L1: ILo, IHi, NAct, IBlock = ', ILo, IHi,
c$$$     $        NAct, IBlock
c$$$C
C        Read gamma block
C
         If ( DRA_Read_Section(.FALSE.,
     $        G_Gam, 1,               NI*NVir, 1, NFit,
     $        D_Gam, (ILo-1)*NVir+1, IHi*NVir, 1, NFit, ReqI) .ne. 0)
     $        Call ErrQuit('RIMP2_L1: Error reading gamma', ILo,
     &       DISK_ERR)
         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
     $        'RIMP2_L1: Error from DRA_Wait on gamma', 1, CALC_ERR)
C
         Do PLo = 1, NFrzO+NAct+NVir, PBlock
            PHi = Min( NFrzO+NAct+NVir, PLo + PBlock - 1)
            NP = PHi - PLo + 1
C
c$$$            Write (6, *) 'L1: PLo, PHi, (limit), PBlock = ', PLo, PHi,
c$$$     $           NFrzO+NAct+NVir, PBlock
c$$$C
C        Read integral block
C
            If ( DRA_Read_Section(.FALSE.,
     $           G_Int, 1,              NP*NVir,  1, NFit,
     $           D_Int, (PLo-1)*NVir+1, PHi*NVir, 1, NFit, ReqP) .ne. 0)
     $           Call ErrQuit('RIMP2_L1: Error reading integrals', PLo,
     &       DISK_ERR)
            If ( DRA_Wait(ReqP) .ne. 0) Call ErrQuit(
     $           'RIMP2_L1: Error from DRA_Wait on integrals', 1,
     &       CALC_ERR)
C
C           Now process all gamma and integrals in core
C
            Call GA_Sync
            Do P = 1, PHi - PLo + 1
               Do I = 1, IHi - ILo + 1
C
                  L1 = GA_DDot_Patch(
     $                 G_Int, 'N', (P-1)*NVir+1, P*NVir, 1, NFit,
     $                 G_Gam, 'N', (I-1)*NVir+1, I*NVir, 1, NFit)
C
C                 Caution:  This may be too much traffic at once.  Check!
C
                  If ( L1 .ne. 0.0d0 .and. ga_nodeid().eq.0 ) then
                     Call GA_Acc( g_L1, P+PLo-1, P+PLo-1, I+ILo-1,
     $                    I+ILo-1, L1, 1, 2.0d0)
                  EndIf
               EndDo            ! I
            EndDo               ! P
            Call GA_Sync
C
         EndDo                  ! P blocks
      EndDo                     ! I blocks
C
C     Free memory
C
      If ( .NOT. GA_Destroy( g_Gam ) ) Call ErrQuit(
     $     'RIMP2_L1: can''t destroy gamma GA', g_Gam, GA_ERR)
      If ( .NOT. GA_Destroy( g_Int ) ) Call ErrQuit(
     $     'RIMP2_L1: can''t destroy integral GA', g_Gam, GA_ERR)
C
C     Print results if requested
C
      If ( Util_Print('partial l', Print_Debug) ) Call GA_Print( g_L1 )
C
      Return
      End
