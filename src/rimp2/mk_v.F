C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     mk_v -- Form V^{-1/2}
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function Mk_V(Basis, MinEval, G_V)
      Implicit NONE
      Integer Basis, G_V
      Double Precision MinEval
C
C ARGUMENTS
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
#include "tcgmsg.fh"
C
C EXTERNAL ROUTINES
      Integer NxTask
      External NxTask     ! From TCGMSG
      Logical Panel_Diagonal
      External Panel_Diagonal
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MinChunk, MinTasks
      Parameter (MinChunk = 10, MinTasks = 10)
C
C LOCAL VARIABLES
      Logical Status, Have_Diag
      Integer NSh, NBF, ISh, JSh, IBFlo, IBFhi, JBFlo, JBFhi
      Integer LEri, H_Eri, ERI_Idx, LScr, H_Scr, Scr_Idx
      Integer H_Eval, Eval_Idx
      Integer IJ, Next
      Integer D, DiagLo, DiagHi
      Integer G_EV, G_Tmp
      Integer I, G_Junk
      Double precision Elem
      Integer Junk
C
C     ****************************
C     * Get basis set parameters *
C     ****************************
C
      Status = .TRUE.
      Status = Status .AND. Bas_NumCont( Basis, NSh)
      Status = Status .AND. Bas_NumBF( Basis, NBF)
      If (.NOT. Status ) then
         Write (LuErr, *) 'Failed to get expansion basis info.'
         Mk_V = .FALSE.
         Return
      EndIf
C
C     ******************************
C     * Generate V integral matrix *
C     ******************************
C     This is the matrix that will eventually be returned to the caller.
C
C     Don't worry about shell or atom-based blocking at the moment,
C     since this shouldn't be a time-consuming routine.
C
      If ( .NOT. GA_Create( MT_Dbl, NBF, NBF, 'V integral matrix',
     $   MinChunk, MinChunk, G_V) ) then
         Write (LuErr, *) 'Unable to allocate V integral matrix'
         Mk_V = .FALSE.
         Return
      EndIf
C
      Call GA_Zero(G_V)
C
C     Figure out memory needs for integral evaluation and allocate the
C     necessary arrays.
C
      Call Int_Mem_2e2c( LEri, LScr)
C
      Status = .TRUE.
      Status = Status .AND.
     $   MA_Push_Get( MT_Dbl, LEri, 'Integral values', H_ERI, ERI_Idx)
      Status = Status .AND.
     $   MA_Push_Get( MT_Dbl, LScr, 'Integral scratch', H_Scr, Scr_Idx)
C
      If (.NOT. Status ) then
         Write (LuErr, *) 'Failed to get temporary space for integrals'
         Mk_V = .FALSE.
         Return
      EndIf
C
C     This algorithm is dynamically load balanced, but may result in
C     as much as N^2 data movement.  It is not worth going with anything
C     more sophisticated unless timings indicate problems, and changes
C     here should probably be coupled with a different distribution for
C     V, as declared above.
C
      IJ = 0
      Next = NxTask( GA_NNodes(), MinTasks)
C
      Do JSh = 1, NSh
         If ( .NOT. Bas_Cn2Bfr( Basis, JSh, JBFlo, JBFhi)) then
            Write (LuErr, *) 'Unable to get info on Jsh ', Jsh
            Mk_V = .FALSE.
            Return
         EndIf
C
         Do ISh = 1, NSh
            If ( .NOT. Bas_Cn2Bfr( Basis, ISh, IBFlo, IBFhi)) then
               Write (LuErr, *) 'Unable to get info on Ish ', Ish
               Mk_V = .FALSE.
               Return
            EndIf
C
c$$$            Write (LuOut, *) 'ISh, JSh, NSh, IJ, Next ',
c$$$     $         ISh, JSh, NSh, IJ, Next
            If ( IJ .eq. Next ) then
               Call Int_2e2c( Basis, ISh, Basis, JSh,
     $            LScr, Dbl_MB( Scr_Idx), LEri, Dbl_MB( Eri_Idx) )
C
C              Plunk this block into the proper location
C
               Call GA_Put(G_V, IBFlo, IBFhi, JBFlo, JBFhi, 
     $            Dbl_MB( ERI_Idx), IBFhi-IBFlo+1)
c$$$C
c$$$               If ( GA_NodeID() .eq. 0)
c$$$     $            Call DGeWr(IBFHi-IBFLo+1, JBFHi-JBFLo+1,
c$$$     $            Dbl_MB(ERI_Idx), IBFHi-IBFLo+1, IBFLo, JBFLo, LuOut,
c$$$     $            'integrals', 80, 'E15.6', Junk)
C
               Next = NxTask( GA_NNodes(), MinTasks)
            EndIf
            IJ = IJ + 1
         EndDo
      EndDo
      Next = NxTask( -GA_NNodes(), MinTasks)
      Call GA_Sync
C
      Status = .TRUE.
      Status = Status .AND. MA_Pop_Stack( H_Scr )
      Status = Status .AND. MA_Pop_Stack( H_ERI )
C
      If ( .NOT. Status ) then
         Write (LuErr, *) 'Unable to free integral temporaries.'
         Mk_V = .FALSE.
         Return
      EndIf
C
      Call GA_Print( G_V)
C
      Mk_V = .TRUE.
C
      Return
      End
