C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_NS2_I3G3 -- Contraction of 3-c integrals with Gamma3 for
C     non-separable 2-c contribution to gradient
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_NS2_I3G3(NFrzO, NAct, NVir, NFit,
     $     D_Int, D_Gam, g_V)

      Implicit NONE
#include "errquit.fh"
      Integer NFrzO
      Integer NAct
      Integer NVir
      Integer NFit
      Integer D_Int
      Integer D_Gam
      Integer g_V
C
C DESCRIPTION
C
C     Performs contraction: {(ia|D) V-1(D,D')} G(T,ia)
C
C NOTES
C     Equations have been derived only for the 'V' approximation
C     Assumes that the ket integrals are actually (ia|D) V-1(D,D')
C
C PRINT CONTROLS
C EXTERNAL ROUTINES
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
C LOCAL VARIABLES
      Integer g_Int, g_Gam, ReqI, ReqG
      Integer NTotal, IBlock
      Double Precision Initial, Final, Block
      Integer ILo, IHi, NI
C
C STATEMENT FUNCTIONS
C
C
C     ***************************************************
C     * Determine memory available for integrals/gammas *
C     ***************************************************
C
      Initial = GA_Memory_Avail()
C
      If ( .NOT. GA_Create( MT_Dbl, NVir, NFit,
     $     'test for memory allocation', -1, -1, G_Int) )
     $     Call ErrQuit(
     $     'RIMP2_NS2_I3G3: can''t allocate minimal GA', NVir*NFit)
C
      Final = GA_Memory_Avail()
      Block = Initial - Final
C
      If( .NOT. GA_Destroy( G_Int) ) Call ErrQuit(
     $     'RIMP2_NS2_I3G3: can''t destroy GA!', 0)
C
C     If this node has no data, then it can hold the full load of
C     integrals and gammas.  We really want the minimum any node can
C     handle.
C
      If ( Block .gt. 0) then
         NTotal = Int( Initial / Block )
      Else
         NTotal = 2 * NAct
      EndIf
      Call GA_IGOp(Msg_RIMP2_Mem, NTotal, 1, 'min')
      IBlock = NTotal/2
C
C     Allocate GAs for integrals and gammas.  Both must use identical
C     blocking
C
      If ( .NOT. GA_Create( MT_Dbl, IBlock*NVir, NFit,
     $     'integrals', -1, -1, G_Int) ) Call ErrQuit(
     $     'RIMP2_NS2_I3G3: can''t allocate integral GA',
     $     IBlock*NVir*NFit)
      If ( .NOT. GA_Create( MT_Dbl, IBlock*NVir, NFit,
     $     'gammas', -1, -1, G_Gam) ) Call ErrQuit(
     $     'RIMP2_NS2_I3G3: can''t allocate integral GA',
     $     IBlock*NVir*NFit)
C
C     Loop over chunks that fit into memory
C
      Do ILo = 1, NAct, IBlock
         IHi = Min( NAct, ILo + IBlock - 1)
         NI = IHi - ILo + 1
C
C        Read integral block.  Account for frozen orbitals.
C
         If ( DRA_Read_Section(.FALSE.,
     $        G_Int, 1,                    NI*NVir,  1, NFit,
     $        D_Int, (NFrzO+ILo-1)*NVir+1, (NFrzO+IHi)*NVir,
     $        1, NFit, ReqI) .ne. 0) Call ErrQuit(
     $        'RIMP2_NS2_I3G3: Error reading integrals', ILo)
C
C        Read gamma block
C
         If ( DRA_Read_Section(.FALSE.,
     $        G_Gam, 1,               NI*NVir, 1, NFit,
     $        D_Gam, (ILo-1)*NVir+1, IHi*NVir, 1, NFit, ReqG) .ne. 0)
     $        Call ErrQuit(
     $        'RIMP2_NS2_I3G3: Error reading gamma', ILo)
C
C        Insure data is ready
C
         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
     $        'RIMP2_NS2_I3G3: Error from DRA_Wait on integrals', 1)
         If ( DRA_Wait(ReqG) .ne. 0) Call ErrQuit(
     $        'RIMP2_NS2_I3G3: Error from DRA_Wait on gamma', 1)
C
C        Multiply integrals and gamma
C
         Call GA_DGEMM('T', 'N', NFit, NFit, NI*NVir,
     $        1.0d0, G_Int, G_Gam, 1.0d0, g_V)
C
      EndDo                     ! ILo
C     
C     Free integral and gamma GAs
C     
      If ( .NOT. GA_Destroy(g_Gam) ) Call ErrQuit(
     $     'RIMP2_NS2_I3G3: can''t free gamma GA', 0)
      If ( .NOT. GA_Destroy(g_Int) ) Call ErrQuit(
     $     'RIMP2_NS2_I3G3: can''t free integral GA', 0)
C
      Return
      End
