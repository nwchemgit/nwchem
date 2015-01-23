C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     sum_pairs -- Sum up a matrix of pair energies to give an MBPT(2)
C     energy.
C
C REVISION
C     $Id$
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
      Integer ILo, IHi, JLo, JHi, LD,  I_Pair
#ifdef BAD_GACCESS
      integer l_local
#endif

C
      Sum_Pairs = FP_Zero
C
C     Find out what patch we have and if there is anything, grab it.
C
      Call GA_Distribution(G_PairE, GA_NodeID(), ILo, IHi, JLo, JHi)
      If ( (IHi-ILo+1) * (JHi-JLo+1) .gt. 0) then
#ifdef  BAD_GACCESS
      if(.not.ma_push_get(MT_DBL,(ihi-ilo+1)*(jhi-jlo+1),
     $  'scratch buff', l_local, i_pair)) call
     $  errquit('sumpairs: pushget failed',0,0)
      ld=ihi-ilo+1
      call ga_get(g_paire,ilo,ihi,jlo,jhi,dbl_mb(i_pair),ld)
#else
         Call GA_Access(G_PairE, ILo, IHi, JLo, JHi, I_Pair, LD)
#endif
C
C        For simplicity, we just loop over the whole patch and pick
C        up the bits we need.
C
      call sum_pairs_0(G_PairE, ILo, IHi, JLo, JHi, 
     ,dbl_mb(I_Pair), LD, Sum_Pairs, LowTri)
#ifdef  BAD_GACCESS
      if(.not.ma_pop_stack(l_local)) call 
     $  errquit('sumpairs: popstack failed',0,0)
#else
         Call GA_Release(G_PairE, ILo, IHi, JLo, JHi)
#endif
      EndIf
C
C     Collect up everyone-s partial sums so we all return the right 
C     answer.
C
      Call GA_DGOp(Msg_PrSum, Sum_Pairs, 1, '+')
      Call GA_Sync

C
      Return
      End
      subroutine sum_pairs_0(G_PairE, ILo, IHi, JLo, JHi, 
     ,Pair, LD, Sum_Pairs, LowTri)
      implicit none
#include "numerical_constants.fh"
      integer G_PairE, ILo, IHi, JLo, JHi, LD
      double precision Pair(*)
      logical LowTri
      double precision  Sum_Pairs
c
      integer Off,J,I
c
       
         Off = 0
         Do J = JLo, JHi
            Do I = ILo, IHi
               Off = Off + 1
C
C              If we are using only the lower triangle, the diagonal
C              comes in with a factor of one in front, and the lower
C              triangle with a factor of two.  Nothing comes from the
C              upper triangle.
C
               If ( LowTri .AND. I .eq. J) then
                  Sum_Pairs = Sum_Pairs + Pair(Off)
               ElseIf (LowTri .AND. I .gt. J) then
                  Sum_Pairs = Sum_Pairs + FP_Two * Pair(Off)
C
C                 If there is no permutational symmetry, everyone
C                 contributes equally.
C
               ElseIf ( .NOT. LowTri) then
                  Sum_Pairs = Sum_Pairs + Pair(Off)
               Endif
C
            EndDo
            Off = Off + LD - (IHi-ILo+1)
         EndDo
C
        return
	end
