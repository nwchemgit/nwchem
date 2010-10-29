C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     GA_MatPow -- Form V^{-1/2}
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine GA_MatPow(G_V, Pow, MinEval)
      Implicit NONE
#include "errquit.fh"
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
C     This routine responds to the print directive "ga_matpow details",
C     which is at the Print_Debug level by default.  The print control
C     context is that of the caller.
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
#include "util.fh"
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
      Integer I, J, VType, Rows, Cols, ILo, IHi, JLo, JHi
#if defined(USE_GEN_DIAG)
      Integer G_Junk
#endif
      Double precision Elem
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
     $   Int( 100 * Pow ) , GA_ERR)
C
      If ( Util_Print('ga_matpow details', Print_Debug) ) then
         If ( GA_NodeID() .eq. 0 )
     $        Call Util_Print_Centered(LuOut, 'GA_MatPow Input', 40,
     $        .TRUE.)
         Call GA_Sync
         Call GA_Print( G_V )
      EndIf
C
C     *********************************
C     * Get basic info about input GA *
C     *********************************
C
      Call GA_Inquire( G_V, VType, Rows, Cols)
      If ( Rows .ne. Cols )
     $   Call Errquit('GA_MatPow: Matrix must be square', Rows-Cols,
     &       GA_ERR)
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
     $   (Rows+1)*Cols, GA_ERR)
C
C     Now we can solve the stupid thing...
C
#if defined(USE_GEN_DIAG)
C     Create an array for the metric and turn it into the unit matrix
C
      If ( .NOT. GA_Create( VType, Rows, Cols, 'Metric',
     $   -MinChunk, -MinChunk, G_Junk) )
     $   Call ErrQuit('GA_MatPow: Unable to allocate metric',
     $   Rows*Cols, GA_ERR)
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
#ifdef SCALAPACK
         call ga_pdsygv(g_v, g_junk, g_ev, dbl_mb(eval_idx),0)
#else
      Call GA_Diag( G_V, G_Junk, G_EV, Dbl_MB(Eval_Idx))
#endif
#else
      Call GA_Diag_Seq( G_V, G_Junk, G_EV, Dbl_MB(Eval_Idx))
#endif
C
      If ( .NOT. GA_Destroy( G_Junk) )
     $   Call ErrQuit('GA_MatPow: Unable to free metric', 0, GA_ERR)
C
#else /* defined(USE_GEN_DIAG) */
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
      call ga_pdsyev(g_v, g_ev, dbl_mb(eval_idx),0)
#else
      Call GA_Diag_Std( G_V, G_EV, Dbl_MB(Eval_Idx))
#endif
#else
      Call GA_Diag_Std_Seq( G_V, G_EV, Dbl_MB(Eval_Idx))
#endif
#endif /* defined(USE_GEN_DIAG) */
C
C     Sometimes people may just want a look at the eigenvalues.
C
      If ( Util_Print('ga_matpow eigenvalues', Print_Debug)
     $   .OR. Util_Print('ga_matpow details', Print_Debug) ) then
         If ( GA_NodeID() .eq. 0 ) then
            Call Util_Print_Centered(LuOut, 'GA_MatPow Eigenvalues', 40,
     $           .TRUE.)
            Do I = 1, Rows, 4
               Write (LuOut, '(1X, I5, 4(1X, E15.6) )') I,
     $            (Dbl_MB(Eval_Idx + J - 1), J=I,Min(I+3, Rows) )
            EndDo
         EndIf
         Call GA_Sync
      EndIf
      If ( Util_Print('ga_matpow details', Print_Debug) ) then
         If ( GA_NodeID() .eq. 0 ) then
            Call Util_Print_Centered(LuOut, 'GA_MatPow Vectors', 40,
     $           .TRUE.)
         EndIf
         Call GA_Sync
         Call GA_Print( G_EV )
      EndIf
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
C     Consider the diagonal of a matrix.  It divides the upper triangle,
C     where I (row) < J (column) from the lower triangle, where I > J.
C
C     If IHi-JLo > 0 the lower-left corner of the patch is below the diagonal.
C     If ILo-JHi < 0 the upper-right corner of the patch is above the diagonal.
C     If the diagonal passes through the patch, these two corners must
C     be on opposite sides of the diagonal, so the two differences must be
C     of opposite sign, or their product must be negative.
C
C     If the diagonal hits precisely on the ll or ur corner, the
C     corresponding difference will be zero, and the product will
C     be zero.  This we must also catch.
C
      Have_Diag = ( IHi-JLo ) * ( ILo-JHi ) .le. 0
C
      DiagLo = Max( ILo, JLo )
      DiagHi = Min( IHi, JHi )
C
      If ( Have_Diag ) then
         Do D = DiagLo, DiagHi
            Elem = Dbl_MB( Eval_Idx + D - 1) 
C
C           It seems that in the RI-MP2 application, very large
C           fitting basis sets can lead to a fair number of eigenvalues
C           which are negative, but small (< ~1e-5).  We suppose that
C           these arise from small numerical inaccuracies in the
C           integrals, and perhaps in the eigensolver.  In any case,
C           requiring strict positive definiteness seems too harsh
C           in practice.
C
c$$$C
c$$$            If ( Elem .lt. FP_Zero ) then
c$$$               Write (LuErr, *) 'V is not positive definite! ',
c$$$     $            D, Elem
c$$$               GA_MatPow = .FALSE.
c$$$               Return
c$$$            EndIf
C
            If ( Pow .lt. FP_Zero .AND. Elem .lt. MinEval
     $           .OR. Elem .lt. FP_Zero) then
               Write (LuErr, *)
     $            'Dropping small eigenvalue ', D, Elem
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
C
      If ( Util_Print('ga_matpow details', Print_Debug) ) then
         If ( GA_NodeID() .eq. 0 )
     $        Call Util_Print_Centered(LuOut, 'GA_MatPow Diagonal', 40,
     $           .TRUE.)
         Call GA_Sync
         Call GA_Print( G_V )
      EndIf
C
C     Now we are done with the local eigenvalue array
C
      If ( .NOT. MA_Pop_Stack( H_Eval ) )
     $   Call ErrQuit('GA_MatPow: Unable to free local eigenvalues', 0,
     &       MA_ERR)
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
     $   Rows*Cols, GA_ERR)
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
     $   Call ErrQuit('GA_MatPow: Unable to free global arrays', 0,
     &       GA_ERR)
C
      If ( Util_Print('ga_matpow details', Print_Debug) ) then
         If ( GA_NodeID() .eq. 0 )
     $        Call Util_Print_Centered(LuOut, 'GA_MatPow Output', 40,
     $           .TRUE.)
         Call GA_Sync
         Call GA_Print( G_V )
      EndIf
C
      Return
      End
