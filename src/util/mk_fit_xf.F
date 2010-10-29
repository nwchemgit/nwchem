C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     mk_fit_xf -- Form fitting basis transformation
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function Mk_Fit_Xf(Approx, Split, Basis, MinEval, G_V)
      Implicit NONE
#include "errquit.fh"
      Character*(*) Approx, Split
      Integer Basis, G_V
      Double Precision MinEval
C
C ARGUMENTS
C     Approx   RI approximation used (SVS, S, or V) [IN]
C     Split    Whether or not to return the square root of the matrix
C              so that it can be used to transform both sets of 3c ints.
C              Y or N [IN]
C     Basis    Handle for basis set in which to evaluate V [IN]
C     MinEval  Minimum eigenvalue of V matrix to be retained in 
C              the inversion [IN]
C     G_V      Global array handle to the V^{-1/2} matrix [OUT]
C
C RETURN VALUE
C     .TRUE.  if successful, even if some eigenvalues fell below MinEval.
C     .FALSE. if errors occured in dynamic memory (MA or GA) operations,
C             inquiries about Basis, or in obtaining the required integrals.
C
C DESCRIPTION
C     Forms V^{-1/2} for the specified basis set and returns it in 
C     the newly created global array G_V.
C
C NOTES
C     The integral package must be initialized before calling this routine.
C
C MEMORY USE
C     Creates and returns a global array (G_V) the size of 
C     Bas_NumBF(Basis)^2.
C
C     Additional temporary usage consists of the largest of:
C     1) Integral requirements, reported by Int_Mem_2e2c.
C     2) Bas_NumBF(Basis)^2 + Bas_NumBF(Basis) + whatever additional
C        space is required by GA_Diag_Std.
C     3) 2 * Bas_NumBF(Basis)^2
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "numerical_constants.fh"
#include "stdio.fh"
#include "bas.fh"
#include "inp.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MinChunk, MinTasks
      Parameter (MinChunk = 10, MinTasks = 10)
C
C LOCAL VARIABLES
      Logical Status, DoSplit
      Integer NBF, Type, G_Tmp, G_S
C
      If ( Inp_Compare(.FALSE., Approx, 'SVS') ) then
         Type = 1
      ElseIf ( Inp_Compare(.FALSE., Approx, 'S') ) then
         Type = 2
      ElseIf ( Inp_Compare(.FALSE., Approx, 'V') ) then
         Type = 3
      Else
         Write (LuErr, 9000) Approx
         Call ErrQuit('Mk_Fit_XF: Bad argument', 1, INPUT_ERR)
      EndIf
      If ( Inp_Compare(.FALSE., Split, 'Y') ) then
         DoSplit = .TRUE.
      ElseIf ( Inp_Compare(.FALSE., Split, 'N') ) then
         DoSplit = .FALSE.
      Else
         Write (LuErr, 9010) Split
         Call ErrQuit('Mk_Fit_XF: Bad argument', 2, INPUT_ERR)
      EndIf
 9000 Format(1X, 'Mk_Fit_XF: Unknown approximation requested: ', A)
 9010 Format(1X, 'Mk_Fit_XF: ''split'' argument is Y/N: ', A)
C
C     ****************************
C     * Get basis set parameters *
C     ****************************
C
      Status = .TRUE.
      Status = Status .AND. Bas_NumBF( Basis, NBF)
      If (.NOT. Status )
     $   Call ErrQuit('Mk_Fit_XF: Failed to get expansion basis info',
     $   0, BASIS_ERR)
C
      If (Type .eq. 1) then
C
C        **************************************************
C        * Produce S^{-1}VS^{-1} or (S^{-1}VS^{-1})^{1/2} *
C        **************************************************
C        This is the matrix that will eventually be returned to the caller.
C
C        Do not worry about shell or atom-based blocking at the moment,
C        since this should not be a time-consuming routine.
C
         If ( .NOT. GA_Create( MT_Dbl, NBF, NBF, 'S integral matrix',
     $      -MinChunk, -MinChunk, G_S) )
     $      Call ErrQuit('Mk_Fit_XF: Unable to allocate S temporary',
     $      0, GA_ERR)
C
         Call GA_Zero(G_S)
C
C        The GA_Zero is SIMD, so it should sync us all
C
         Call Int_2c_GA( Basis, Basis, G_S, 'overlap', FP_One)
         Call GA_Sync
C
C        Form S^-1
C
         Call GA_MatPow( G_S, -FP_One, MinEval)
C
C        Form V
C
         If ( .NOT. GA_Create( MT_Dbl, NBF, NBF, 'SVS integral matrix',
     $      -MinChunk, -MinChunk, G_V) ) 
     $      Call ErrQuit('Mk_Fit_XF: Unable to allocate SVS',
     $      0, GA_ERR)
C
         Call GA_Zero(G_V)
C
C        The GA_Zero is SIMD, so it should sync us all
C
         Call Int_2c_GA( Basis, Basis, G_V, 'eri', FP_One)
         Call GA_Sync
C
C        Multiply everyone together
C
         If ( .NOT. GA_Create( MT_Dbl, NBF, NBF, 'tmp integral matrix',
     $      -MinChunk, -MinChunk, G_Tmp) )
     $      Call ErrQuit('Mk_Fit_XF: Unable to allocate temporary',
     $      0, GA_ERR)
C
         Call GA_DGEMM('N', 'N', NBF, NBF, NBF,
     $      FP_One, G_S, G_V, FP_Zero, G_Tmp)
C
         Call GA_DGEMM('N', 'N', NBF, NBF, NBF,
     $      FP_One, G_Tmp, G_S, FP_Zero, G_V)
C
         Status = .TRUE.
         Status = Status .AND. GA_Destroy(G_Tmp)
         Status = Status .AND. GA_Destroy(G_S)
C
         If (.NOT. Status)
     $      Call ErrQuit('Mk_Fit_XF: Unable to deallocate temporaries',
     $      0, GA_ERR)

      ElseIf (Type .eq. 2) then
C
C        ******************************
C        * Produce S^{-1} or S^{-1/2} *
C        ******************************
C        This is the matrix that will eventually be returned to the caller.
C
C        Do not worry about shell or atom-based blocking at the moment,
C        since this should not be a time-consuming routine.
C
         If ( .NOT. GA_Create( MT_Dbl, NBF, NBF, 'S integral matrix',
     $      -MinChunk, -MinChunk, G_V) ) 
     $      Call ErrQuit('Mk_Fit_XF: Unable to allocate S',
     $      0, GA_ERR)
C
         Call GA_Zero(G_V)
C
C        The GA_Zero is SIMD, so it should sync us all
C
         Call Int_2c_GA( Basis, Basis, G_V, 'overlap', FP_One)
         Call GA_Sync
C
      ElseIf (Type .eq. 3) then
C
C        ******************************
C        * Produce V^{-1} or V^{-1/2} *
C        ******************************
C        This is the matrix that will eventually be returned to the caller.
C
C        Do not worry about shell or atom-based blocking at the moment,
C        since this should not be a time-consuming routine.
C
         If ( .NOT. GA_Create( MT_Dbl, NBF, NBF, 'V integral matrix',
     $      -MinChunk, -MinChunk, G_V) )
     $      Call ErrQuit('Mk_Fit_XF: Unable to allocate V',
     $      0, GA_ERR)
C
         Call GA_Zero(G_V)
C
C        The GA_Zero is SIMD, so it should sync us all
C
         Call Int_2c_GA( Basis, Basis, G_V, 'eri', FP_One)
         Call GA_Sync
      EndIf
C
C     ****************************************
C     * Compute M^{1/2}, M^{-1/2}, or M^{-1} *
C     ****************************************
C     As appropriate to the approximation and "split" requested
C
      If (DoSplit .AND. Type .eq. 1) then
         Call GA_MatPow( G_V, FP_Half, MinEval)
      ElseIf ( DoSplit .AND. (Type .eq. 2 .OR. Type .eq. 3) ) then
         Call GA_MatPow( G_V, -FP_Half, MinEval)
      ElseIf ( .NOT. DoSplit .AND. (Type .eq. 2 .OR. Type .eq. 3) ) then
         Call GA_MatPow( G_V, -FP_One, MinEval)
      EndIf
C
      Mk_Fit_Xf = .TRUE.
C
      Return
      End
