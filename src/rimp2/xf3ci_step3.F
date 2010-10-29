C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     XF3cI_Step3 -- Do the last step of the three-center integral 
C     transformation
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function XF3cI_Step3(G_In, NAIn, NIIn,
     $   AInlo, AInhi, IInlo, IInhi, SInlo, SInHi,
     $   G_C, SOutLo, SOutHi,
     $   G_Out, NAOut, NIOut, AOut, IOut, SOut)
      Implicit NONE
      Integer NAIn, NIIn, AInLo, AInHi, IInLo, IInHi, SInlo, SInHi
      Integer SOutLo, SOutHi, NAOut, NIOut, AOut, IOut, SOut
      Integer G_In, G_C, G_Out
C
C ARGUMENTS
C     G_In     GA handle for incoming (ia|s) integrals [IN]
C     NAIn     Size of A range in G_In [IN]
C     NIIn     Size of I range in G_In [IN]
C     AInLo    Lower limit of A range in G_In to transform [IN]
C     AInHi    Upper limit of A range in G_In to transform [IN]
C     IInLo    Lower limit of I range in G_In to transform [IN]
C     IInHi    Upper limit of I range in G_In to transform [IN]
C     SInLo    Lower limit of S range in G_In to transform [IN]
C     SInHi    Upper limit of S range in G_In to transform [IN]
C     G_C      GA handle for fitting basis transformation matrix [IN]
C     SOutLo   Lower limit of S in columns of G_C or in G_Out [IN]
C     SOutHi   Upper limit of S in columns of G_C or in G_Out [IN]
C     G_Out    GA handle for resulting (ia|s') [IN]
C     NAOut    Size of A range in G_Out [IN]
C     NIOut    Size of I range in G_Out [IN]
C     AOut     Destination A in G_Out for transformed integrals [IN]
C     IOut     Destination I in G_Out for transformed integrals [IN]
C     SOut     Destination S in G_Out for transformed integrals [IN]
C
C RETURN VALUE
C     .TRUE.   Successful completion
C     .FALSE.  Error condition.  Could be in input, basis info queries,
C              or memory operations (MA or GA).
C DESCRIPTION
C     Performs the transformation on the fitting basis index of the
C     three center integrals.  Both input and output integral matrices
C     are organized with the combined (ia| index in the first dimension
C     and the fitting basis index |s) as the second dimension.  Selected
C     regions of the input matrix can be specified.
C
C     Note that because of the index collapse in the first dimension
C     of the integral matrices, it is necessary to ask for the individual
C     sizes of A and I in order to do the calculations.  But the size
C     of S can be determined from the global arrays themselves, to
C     they do not need to be input.
C
C LIMITATIONS
C     Currently the transformed region must be the entire matrix.  This
C     restriction can be lifted when GA_DGEMM is capable of working
C     with submatrices.
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "numerical_constants.fh"
#include "stdio.fh"
#include "rimp2_timers.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer NAXF, NIXF, NSXFOut, Dim1In, Dim1Out, Type,
     $  NCRow, NCCol, NSIn, NSOut
      Logical Status
C
C     This routine is essentially synchronous anyway because GA_DGEMM
C     is SIMD, so to be fair, we should include only the time after
C     everyone meets up here.
C
      Call GA_Sync
      Call PStat_On(PS_Step3)
C
C     For convenience, determine the size of the region to be transformed
C
      NAXF = AInHi - AInLo + 1
      NIXF = IInHi - IInLo + 1
      NSXFOut = SOutHi - SOutLo + 1
c$$$C
c$$$      Write (LuOut, *) 'Step3 inputs ', GA_NodeID(), NAIn, NIIn,
c$$$     $   AInLo, AInHi, IInLo, IInHi, SInLo, SInHi, SOutLo, SOutHi,
c$$$     $   NAOut, NIOut, AOut, IOut, SOut
C
C     *****************
C     * Sanity checks *
C     *****************
C     Verify that we've got good GA & basis set handles
C
      Call GA_Check_Handle ( G_In, 'XF3cI_Step3: Input integrals')
      Call GA_Check_Handle ( G_Out, 'XF3cI_Step3: Output integrals')
      Call GA_Check_Handle ( G_C, 'XF3cI_Step3: Transformation')
C
C     Verify that the dimensions declared are consistent with arguments
C     and obtain the size of the S dimension for everyone.
C
      Call GA_Inquire( G_In,  Type, Dim1In, NSIn)
      Call GA_Inquire( G_Out, Type, Dim1Out, NSOut)
      Call GA_Inquire( G_C,   Type, NCRow, NCCol)
C
      If ( Dim1In .lt. NAIn*NIIn .OR. Dim1Out .lt. NAOut*NIOut ) then
         Write (LuErr, *) 'XF3cI_Step3: Declared I and A ',
     $      'dimensions inconsistent with GA for G_In or G_Out.'
         XF3cI_Step3 = .FALSE.
         Return
      EndIf
C
C     Verify that the specified regions fall within the matrices
C
      Status = .TRUE.
      Status = Status .AND. AInLo .ge. 1 .AND. AInLo .le. NAIn
      Status = Status .AND. AInHi .ge. 1 .AND. AInHi .le. NAIn
      Status = Status .AND. IInLo .ge. 1 .AND. IInLo .le. NIIn
      Status = Status .AND. IInHi .ge. 1 .AND. IInHi .le. NIIn
      Status = Status .AND. SInLo .ge. 1 .AND. SInLo .le. NSIn
      Status = Status .AND. SInHi .ge. 1 .AND. SInHi .le. NSIn
C
      Status = Status .AND. SInHi .le. NCRow
      Status = Status .AND. SOutHi .ge. 1 .AND. SOutHi .le. NCCol
C
      Status = Status .AND. AOut .ge. 1 .AND. AOut .le. NAOut
      Status = Status .AND. IOut .ge. 1 .AND. IOut .le. NIOut
      Status = Status .AND. SOut .ge. 1 .AND. SOut .le. NSOut
      Status = Status .AND. AOut + NAXF - 1 .le. NAOut
      Status = Status .AND. IOut + NIXF - 1 .le. NIOut
      Status = Status .AND. SOut + NSXFOut - 1 .le. NSOut
C
      If ( .NOT. Status ) then
         Write (LuErr, *) 'XF3cI_Step3: Trans. region ',
     $      'inconsistent with dimensions: G_In, G_C, or G_Out.'
         XF3cI_Step3 = .FALSE.
         Return
      EndIf
C
C     For the time being, make sure we're not asked to do a submatrix
C     that isn't rooted at (1,1)
C
      Status = .TRUE.
      Status = Status .AND. AInLo .eq. 1 .AND. AInHi .eq. NAIn
      Status = Status .AND. IInLo .eq. 1 
      Status = Status .AND. SInLo .eq. 1 
C
      Status = Status .AND. SOutLo .eq. 1
C
      Status = Status .AND. AOut .eq. 1 .AND. IOut .eq. 1
     $   .AND. SOut .eq. 1
C
      If ( .NOT. Status ) then
         Write (LuErr, *) 'XF3cI_Step3: Cannot transform ',
     $      ' submatrices yet!'
         XF3cI_Step3 = .FALSE.
         Return
      EndIf
C
C     *****************************************
C     * A very straightforward operation here *
C     *****************************************
C
      Call GA_DGEMM('N', 'N', NAXF*NIXF, NSOut, NSIn, FP_One, G_In, G_C,
     $   FP_Zero, G_Out)
C
      Call PStat_Off(PS_Step3)
C
      XF3cI_Step3 = .TRUE.
C
      Return
      End

