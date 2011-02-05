C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Asym_I -- Antisymmetrize IJ integral batch in memory
C     
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function RIMP2_Asym_I(TopSpin, SameSpin, ILo, IHi, NVirI,
     $   JLo, JHi, NVirJ, IJSym, IJTri, G_Int, SpinLabel)
      Implicit NONE
#include "errquit.fh"
      Integer TopSpin ![in] Indicates RHF (1) or UHF (2) calculation
      Logical SameSpin![in] Indicates if I and J have the same spin case
      Integer ILo     ![in] Lowest occ index in batch (row)
      Integer IHi     ![in] Highest occ index in batch (row)
      Integer NVirI   ![in] Size of blocks (row)
      Integer JLo     ![in] Lowest occ index in batch (col)
      Integer JHi     ![in] Highest occ index in batch (col)
      Integer NVirJ   ![in] Size of blocks (col)
      Logical IJSym   ![in] Use I <--> J permutational symmetry if possible
      Logical IJTri   ![in] G_I/G_T hold only triangular I > J
      Integer G_Int   ![in] GA handle for integral batch
      Character*(*) SpinLabel ![in] text label for current spin case
C
C DESCRIPTION
C     Antisymmetrizes a batch of integrals "in place" (as far as the
C     caller can tell).  Exact operation depends on type of calculation
C     (RHF/UHF) and whether the spin cases are the same.
C     
C                 Input               Output
C     RHF        (IA|JB)     3 (IA|JB) - 2 (IB|JA)
C     UHF AA     (IA|JB)       (IA|JB) -   (IB|JA)
C     UHF Ab     (IA|jb)       (IA|jb) 
C     UHF bb     (ia|jb)       (ia|jb) -   (ib|ja)
C     
C     Note that for UHF Ab, this is a no-op!
C
C NOTES
C     This routine happily uses a large temporary GA because we assume
C     the output of this routine is going into something like a
C     t amplitude calculation, which requires as much space.
C
C PRINT CONTROLS
C     'asym 4c ints' prints the batch of antisymmetrized integrals
C     at the end of this routine.
C
C INCLUDED FILES
#include "global.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I, J, JTop, NI, NJ, LT1, LT2, Type
      Double Precision FactA, FactB
      Integer G_AsymI
      Character*255 String1
C     
      Call GA_Sync
      If ( SameSpin .OR. TopSpin .eq. 1 ) then
C
C        Prepare a GA for the amplitudes which matches the integrals
C        in shape and distribution.  That way we can capture anything
C        special the caller has done.
C     
         String1 = 'Antisym 4c2e ints ' // SpinLabel
         RIMP2_Asym_I = GA_Duplicate(G_Int, G_AsymI,
     $      String1(:Inp_StrLen(String1)) )
         If ( .NOT. RIMP2_Asym_I ) Return ! Maybe caller can handle failure
C
C        Determine the size of the matrix.  Dimensions must be multiples
C        of NVirI and NVirJ respectively, but beyond that, we don't
C        really care.
C
         Call GA_Inquire( G_Int, Type, LT1, LT2)
         NI = LT1 / NVirI
         NJ = LT2 / NVirJ

         Call GA_Zero(G_AsymI)
C        
C        Use I <--> J symmetry if requested
C
         Do I = 1, NI
            If (IJSym .AND. ILo .eq. JLo .AND. IHi. eq. JHi
     $         .AND. .not. IJTri) then
               JTop = I
            Else
               JTop = NJ
            EndIf
            Do J = 1, JTop
C              
C              So in every case except UHF Ab, we can initialize
C              G_AsymI with G_Int transpose, but on a block-wise basis.
C              
C              Copy from an IJ block of G_Int to the 
C              coresponding IJ block of G_AsymI, but with
C              transposition.  This routine is data parallel
C              in that nodes which don't have any part of
C              the patch to be copies will just fall through
C              and immediately go on to the next loop
C              iteration.
C              
               Call GA_Copy_Patch('T',
     $            G_Int,   (I-1)*NVirI+1, I*NVirI,
     $                     (J-1)*NVirJ+1, J*NVirJ,
     $            G_AsymI, (I-1)*NVirI+1, I*NVirI,
     $                     (J-1)*NVirJ+1, J*NVirJ )
            EndDo
         EndDo
C        
C        Now we must make sure everyone has a consistent
C        copy of G_AsymI before moving on.
C        
         Call GA_Sync
C        
C        Combine appropriate multiples of original and permuted integrals
C        
         FactA =  Dble( 5 - 2 * TopSpin)  !  3 for RHF,  1 for UHF
         FactB = -Dble( 3 - TopSpin)      ! -2 for RHF, -1 for UHF
         Call GA_DAdd( FactA, G_Int, FactB, G_AsymI, G_AsymI)
C
C        Now destroy original integral array and return new one
C
         If ( .NOT. GA_Destroy(G_Int) ) Call ErrQuit(
     $      'RIMP2_Asym_I: can''t destroy original integrals', G_Int,
     &       GA_ERR)
C
         G_Int = G_AsymI
      EndIf
      Call GA_Sync
C     
      If ( Util_Print('asym 4c ints', Print_Debug) ) then
         Write (LuOut, *) 'After RIMP2_Asym_I...'
         Call GA_Print( G_Int )
      EndIf
      Return
      End


