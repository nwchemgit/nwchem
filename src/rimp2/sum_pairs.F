C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     sum_pairs -- Sum up a matrix of pair energies to give an MBPT(2)
C     energy.
C
C REVISION
C     $Id: sum_pairs.F,v 1.2 1995-10-17 08:55:56 d3g681 Exp $
C
C SYNOPSIS
      Double Precision Function Sum_Pairs( G_PairE, LowTri)
      Implicit NONE
      Integer G_PairE
      Logical LowTri
C
C ARGUMENTS
C     G_PairE  GA handle for pair energy matrix [IN]
C     LowTri   Use only lower triangle of matrix [IN]
C
C RETURNS
C     The sum of the elements of the matrix.
C
C DESCRIPTION
C     Returns either the simple sum of all elements of the matrix 
C     (LowTri = .FALSE.), or the total sum derived from the lower
C     triangle of the matrix only, assuming the matrix is symmetric.
C     In this case, the off diagonals contribute with a factor of
C     two to account for their twin the upper triangle.
C
C INCLUDE FILES
#include "numerical_constants.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer ILo, IHi, JLo, JHi, LD, I, J, I_Pair, Off
C
      Sum_Pairs = FP_Zero
C
C     Find out what patch we have and if there is anything, grab it.
C
      Call GA_Distribution(G_PairE, GA_NodeID(), ILo, IHi, JLo, JHi)
      If ( (IHi-ILo+1) * (JHi-JLo+1) .gt. 0) then
         Call GA_Access(G_PairE, ILo, IHi, JLo, JHi, I_Pair, LD)
C
C        For simplicity, we just loop over the whole patch and pick
C        up the bits we need.
C
         Off = -1
         Do J = JLo, JHi
            Do I = ILo, IHi
               Off = Off + 1
C
C              If we're using only the lower triangle, the diagonal
C              comes in with a factor of one in front, and the lower
C              triangle with a factor of two.  Nothing comes from the
C              upper triangle.
C
               If ( LowTri .AND. I .eq. J) then
                  Sum_Pairs = Sum_Pairs + Dbl_MB(I_Pair + Off)
               ElseIf (LowTri .AND. I .gt. J) then
                  Sum_Pairs = Sum_Pairs + FP_Two * Dbl_MB(I_Pair + Off)
C
C                 If there is no permutational symmetry, everyone
C                 contributes equally.
C
               ElseIf ( .NOT. LowTri) then
                  Sum_Pairs = Sum_Pairs + Dbl_MB(I_Pair + Off)
               Endif
C
            EndDo
            Off = Off + LD - (IHi-ILo+1)
         EndDo
C
         Call GA_Release(G_PairE, ILo, IHi, JLo, JHi)
      EndIf
C
C     Collect up everyone's partial sums so we all return the right 
C     answer.
C
      Call GA_DGOp(Msg_PrSum, Sum_Pairs, 1, '+')
      Call GA_Sync
C
      Return
      End
