C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     GA_MatPow -- Form V^{-1/2}
C
C REVISION
C     $Id: ga_matpow.F,v 1.2 1995-10-17 08:56:31 d3g681 Exp $
C
C SYNOPSIS
      Subroutine GA_MatPow(G_V, Pow, MinEval)
      Implicit NONE
      Integer G_V
      Double Precision Pow, MinEval
C
C ARGUMENTS
C     G_V      Global array handle to the matrix of interest [INOUT]
C     Pow      Power to which matrix is to be taken [IN]
C     MinEval  Minimum eigenvalue of input matrix to be retained in 
C              the case of inversion [IN]
C
C DESCRIPTION
C     Forms V^{1/2}, V^{-1/2} or V^{-1} depending on the value of Pow.
C     Uses the eigenvalue decomposition approach.  In the case of
C     Pow < 0, MinEval allows filtering of small, potentially
C     problemmatic eigenvalues.
C
C     The input GA is overwritten with the exponentiated result.  It is
C     NOT guaranteed that the same handle will be returned -- if it is
C     most efficient, the original V may be destroyed and a new GA created
C     to hold the result.
C
C MEMORY USE
C     Uses a GA the size of V, and a local array the size of the number
C     of rows of V.  The eigensolver requires additional memory.
C
C     Due to the use of a generalized eigensolver, an additional GA the
C     size of V is also used.
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "numerical_constants.fh"
#include "stdio.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MinChunk
      Parameter (MinChunk = 64)
C
C LOCAL VARIABLES
      Logical Status, Have_Diag
      Integer H_Eval, Eval_Idx
      Integer D, DiagLo, DiagHi
      Integer G_EV, G_Tmp
      Integer I, G_Junk, VType, Rows, Cols, ILo, IHi, JLo, JHi
      Double precision Elem
      Integer Junk
C
C     ***************************
C     * Check inputs for sanity *
C     ***************************
C
      Call GA_Check_Handle( G_V, 'GA_MatPow: Input matrix')
C
      If ( Pow .ne. FP_Half .AND. Pow .ne. -FP_Half
     $   .AND. Pow .ne. -FP_One) Call ErrQuit(
     $   'GA_MatPow: Don''t know how to deal with this exponent',
     $   Int( 100 * Pow ) )
C
C     *********************************
C     * Get basic info about input GA *
C     *********************************
C
      Call GA_Inquire( G_V, VType, Rows, Cols)
      If ( Rows .ne. Cols )
     $   Call Errquit('GA_MatPow: Matrix must be square', Rows-Cols)
C
C     ************************
C     * Find eigenstuff of V *
C     ************************
C     Create a global array for the eigenvectors...
C     And a local array for the eigenvalues (that's the way the routine
C     does it!)
C
      Status = .TRUE.
      Status = Status .and. GA_Create( VType, Rows, Cols,
     $   'GA_MatPow eigenvectors', -MinChunk, -MinChunk, G_EV)
      Status = Status .and. MA_Push_Get( VType, Rows,
     $   'GA_MatPow eigenvalues', H_EVal, EVal_Idx)
C
      If ( .NOT. Status )
     $   Call ErrQuit('GA_MatPow: Unable to allocate temporaries',
     $   (Rows+1)*Cols)
C
C     Now we can solve the stupid thing...
C
#if defined(USE_GEN_DIAG)
C     Create an array for the metric and turn it into the unit matrix
C
      If ( .NOT. GA_Create( VType, Rows, Cols, 'Metric',
     $   -MinChunk, -MinChunk, G_Junk) )
     $   Call ErrQuit('GA_MatPow: Unable to allocate metric',
     $   Rows*Cols)
C
C     This is clearly serial, but hopefully soon we should be able to
C     go back to the standard eigensolver instead of the generalized one.
C
      Call GA_Zero( G_Junk)
      If ( GA_NodeID() .eq. 0) then
         Do I = 1, Rows
            Call GA_Put(G_Junk, I, I, I, I, FP_One, 1)
         EndDo
      EndIf
C
#if defined(PARALLEL_DIAG)
      Call GA_Diag( G_V, G_Junk, G_EV, Dbl_MB(Eval_Idx))
#else
      Call GA_Diag_Seq( G_V, G_Junk, G_EV, Dbl_MB(Eval_Idx))
#endif
C
      If ( .NOT. GA_Destroy( G_Junk) )
     $   Call ErrQuit('GA_MatPow: Unable to free metric', 0)
C
#else /* defined(USE_GEN_DIAG) */
#if defined(PARALLEL_DIAG)
      Call GA_Diag_Std( G_V, G_EV, Dbl_MB(Eval_Idx))
#else
      Call GA_Diag_Std_Seq( G_V, G_EV, Dbl_MB(Eval_Idx))
#endif
#endif /* defined(USE_GEN_DIAG) */
C
c$$$      If ( GA_NodeID() .eq. 0)
c$$$     $   Call DGeWr(1, Rows, Dbl_MB(Eval_Idx), 1, 0, 0, LuOut,
c$$$     $   'GA_MatPow eigenvalues', 80, 'E15.6', Junk)
c$$$      Call GA_Print(G_EV)
C
C     **************************************************
C     * Filter out potentially problematic eigenvalues *
C     **************************************************
C     Everybody find out what panel of V they hold, and what its
C     diagonal may be.  Then we can have each node work on
C     what it "needs"
C
      Call GA_Distribution( G_V, GA_NodeID(), Ilo, Ihi,
     $   Jlo, Jhi)
C
C     If our panel includes some of the matrix diagonal, this will
C     be true.
C
      Have_Diag = ( IHi-JLo ) * ( ILo-JHi ) .lt. 0
C
      DiagLo = Max( ILo, JLo )
      DiagHi = Min( IHi, JHi )
C
      If ( Have_Diag ) then
         Do D = DiagLo, DiagHi
            Elem = Dbl_MB( Eval_Idx + D - 1) 
c$$$C
c$$$            If ( Elem .lt. FP_Zero ) then
c$$$               Write (LuErr, *) 'V is not positive definite! ',
c$$$     $            D, Elem
c$$$               GA_MatPow = .FALSE.
c$$$               Return
c$$$            EndIf
C
            If ( Pow .lt. FP_Zero .AND. Elem .lt. MinEval ) then
               Write (LuErr, *)
     $            'Setting small eigenvalue to infinity ', D, Elem
               Elem = FP_Zero
            ElseIf (Pow .eq. -FP_Half) then
               Elem = FP_One / Sqrt( Elem )
            ElseIf (Pow .eq. -FP_One) then
               Elem = FP_One / Elem
            ElseIf (Pow .eq. FP_Half) then
               Elem = Sqrt( Elem )
            EndIf
C
            Dbl_MB( Eval_Idx + D - 1) = Elem
C            
         EndDo
      EndIf
C
c$$$      If ( GA_NodeID() .eq. 0)
c$$$     $   Call DGeWr(1, Rows, Dbl_MB(Eval_Idx), 1, 0, 0, LuOut,
c$$$     $   'V eigenvalues^{-1/2}', 80, 'E15.6', Junk)
C
C     Turn V into a diagonal matrix by copying in our bit of the
C     diagonal.
C
      Call GA_Zero( G_V)
C
      If ( Have_Diag ) then
         Do D = DiagLo, DiagHi
            Call GA_Put( G_V, D, D, D, D,
     $         Dbl_MB( Eval_Idx + D - 1), 1)
         EndDo
      EndIf
      Call GA_Sync
c$$$      Write (LuOut, *) 'Diagonal matrix '
c$$$      Call GA_Print(G_V)
C
C     Now we are done with the local eigenvalue array
C
      If ( .NOT. MA_Pop_Stack( H_Eval ) )
     $   Call ErrQuit('GA_MatPow: Unable to free local eigenvalues', 0)
C
C     ****************************
C     * V^{-1/2} = X' v^{-1/2} X *
C     ****************************
C     Note that the eigenvalues, v^{-1/2} is actually stored in the
C     form of a full matrix.  Currently GA offers no way to go in
C     and scale rows/columns, which would be more efficient.
C
      If ( .NOT. GA_Create( VType, Rows, Cols, 'Temporary matrix',
     $   -MinChunk, -MinChunk, G_Tmp) )
     $   Call ErrQuit('GA_MatPow: Unable to allocate temporary matrix',
     $   Rows*Cols)
C
      Call GA_DGEMM( 'N', 'T', Rows, Cols, Cols, FP_One, G_V, G_EV,
     $   FP_Zero, G_Tmp)
      Call GA_DGEMM( 'N', 'N', Rows, Cols, Cols, FP_One, G_EV, G_Tmp,
     $   FP_Zero, G_V)
C
C     Get rid of the global arrays that are not persistent
C
      Status = .TRUE.
      Status = Status .and. GA_Destroy( G_Tmp)
      Status = Status .and. GA_Destroy( G_EV)
C
      If ( .NOT. Status )
     $   Call ErrQuit('GA_MatPow: Unable to free global arrays', 0)
c$$$C
c$$$      Write (Luout, *) 'V^{-1/2} MATRIX'
c$$$      CALL GA_PRINT( G_V)
C
      Return
      End
