C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_e_mem -- Determine optimal memory usage for energy calc.
C
C REVISION
C     $Id: rimp2_e_mem.F,v 1.1 1995-03-29 06:05:47 gg502 Exp $
C
C SYNOPSIS
      Subroutine RIMP2_E_Mem(TopSpin, NOcc, NVir, FitBas,
     $   IBlock, JBlock)
      Implicit NONE
      Integer TopSpin        ! [input]  Number of spin cases
      Integer NOcc(TopSpin)  ! [input]  Size of occupied spaces
      Integer NVir(TopSpin)  ! [input]  Size of virtual spaces
      Integer FitBas         ! [input]  Fitting basis set handle
      Integer IBlock         ! [output] Blocksize for outer loop index
      Integer JBlock         ! [output] Blocksize for inner loop index
C
C DESCRIPTION
C
C EXTERNAL ROUTINES
      Integer iamax
      External iamax
      Logical GA_Create_One_By_Shells
      External GA_Create_One_By_Shells
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.h"
#include "msgids.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer BigOcc, BigVir, I, J, Need1, Need2, G_I, NProc
      Integer Mem_Avail, Final_Avail, IBlock_Size, ABlock_Size
C
C     For simplicity, we'll figure out a distribution that works for
C     the largest spin case.  Note that nrs. of occupied and virtual
C     orbitals may change, but fitting basis doesn't.
C
      BigOcc = IAMax(TopSpin, NOcc, 1)
      BigVir = IAMax(TopSpin, NVir, 1)
C
C     Find out how much space each occupied block of the integral
C     arrays take up.
C
      Mem_Avail = MA_SizeOf(MT_Byte, GA_Memory_Avail(), MT_Dbl)
C
      If ( .NOT. GA_Create_One_By_Shells( MT_Dbl, 1 * BigVir,
     $      FitBas, 'Test Integral matrix', G_I) ) Call ErrQuit(
     $      'RIMP2_E_Mem: Unable to allocate minimal GA', 0)
C
      Final_Avail = MA_SizeOf(MT_Byte, GA_Memory_Avail(), MT_Dbl)
      IBlock_Size = Mem_Avail - Final_Avail
C
      If (.NOT. GA_Destroy( G_I ) ) Call ErrQuit(
     $   'RIMP2_E_Mem: Unable to destroy a GA!', 0)
C
C     The space in a block of the (approximate) four-center matrix
C     is easy to figure out. NOTE:  This is NOT a per-node figure
C     like IBlock_Size because the likelihood of rounding errors
C     being magnified.
C
      ABlock_Size = BigVir * BigVir
C
C     Now we have to take these two values and try to find block
C     sizes that optimize the memory usage.
C
C     Give ourselves a safety margin to cover uneven distributions
C     and the fact that we might compete with later operations for
C     a small amount of memory space.
C
      Mem_Avail = 0.9 * Mem_Avail
C
C     Start from the largest possible memory usage and work back
C     until we get something that fits.  The larger I is (outer loop
C     of energy evaluation) the less I/O we have to do.
C
      NProc = GA_NNodes()
C
      Do I = BigOcc, 1, -1
         Do J = BigOcc, 1, -1
C
C           Per-node memory needs for approx-4-index construction:
C              I block of 3-center integrals
C              J block of 3-center integrals
C              one set of four-center integrals
C
            Need1 = I * IBlock_Size + J * IBlock_Size
     $         + (I * J * ABlock_Size / NProc )
C
C           Per-node memory needs for pair energy evaluation:
C              I block of 3-center integrals (to reduce I/O)
C              two sets of four-center integrals
C
            Need2 = I * IBlock_Size + (2 * I * J * ABlock_Size / NProc )
C
C           See if these will fit in memory
C
            If ( Max( Need1, Need2) .lt. Mem_Avail ) Goto 1000
         EndDo
      EndDo
C
C     Not good -- it looks like we can't do the calculation
C
      Call ErrQuit(
     $   'RIMP2_E_Mem: Insufficient memory for energy evaluation',
     $   Mem_Avail)
C
      Return
C
C     We have acceptable block sizes, now must use the smallest set.
C
 1000 Continue
      IBlock = I
      Call GA_IGOp(Msg_RIMP2_Mem, IBlock, 1, 'min')
      JBlock = J
      Call GA_IGOp(Msg_RIMP2_Mem, JBlock, 1, 'min')
C
      Return
      End
