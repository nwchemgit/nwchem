C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_prep_dra -- Setup DRA sizes and names
C
C REVISION
C     $Id: prep_dra.F,v 1.1 1998-07-06 21:10:23 bernhold Exp $
C
C SYNOPSIS
      Subroutine RIMP2_Prep_DRA(TopSpin, ILo, IHi, ALo, AHi, NFit,
     $     DoGrad, BraKetSame, NrDRA, MaxSiz, TotSiz,
     $     BaseNameI, BaseNameG, NameI, NameG)
      Implicit NONE
      Integer TopSpin ![in]
      Integer ILo(TopSpin) ![in]
      Integer IHi(TopSpin) ![in]
      Integer ALo(TopSpin) ![in]
      Integer AHi(TopSpin) ![in]
      Integer NFit ![in]
      Logical DoGrad ![in]
      Logical BraKetSame ![in]
      Integer NrDRA ![out]
      Integer MaxSiz ![out]
      Integer TotSiz ![out]
      Character*(*) BaseNameI ![in]
      Character*(*) BaseNameG ![in]
      Character*(*) NameI(TopSpin, 2, 2) ![out]
      Character*(*) NameG(TopSpin) ![out]
C
C DESCRIPTION
C
C EXTERNAL ROUTINES
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "inp.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer Size, SizeB
      Integer NrTypes, I
      Integer ISpin, IType
      Character*80 Tmp
C
C STATEMENT FUNCTIONS
      Character*1 SpinItoA
      SpinItoA(ISpin) = Char( ISpin + IChar('a') - 1 )
C
C     *************************************
C     * Integral transformation (to disk) *
C     *************************************
C     Compute largest and total DRA sizes (in bytes)
C
      NrDRA = 0
      MaxSiz = 0
      TotSiz = 0
      Do I = 1, TopSpin
         Size =  (IHi(i)-ILo(i)+1) * (AHi(i)-ALo(i)+1) * NFit
         SizeB = MA_SizeOf(MT_Dbl, Size, MT_Byte)
         MaxSiz = Max( MaxSiz, SizeB)
C
         TotSiz = TotSiz + SizeB
         NrDRA = NrDRA + 1
C
         Tmp = BaseNameI(:Inp_StrLen(BaseNameI)) // '-ab'
     $        // SpinItoA(I)
         Call Util_File_Name(Tmp, .TRUE., .FALSE., NameI(I, 1, 1) )
C
C        If the bra and ket integrals are not equivalent, we'll
C        have two instead of one (same size, though).
C
         If ( BraKetSame ) then
            NameI(I, 2, 1) = NameI(I, 1, 1)
         Else
            TotSiz = TotSiz + SizeB
            NrDRA = NrDRA + 1
C
            Tmp = BaseNameI(:Inp_StrLen(BaseNameI)) // '-ak'
     $           // SpinItoA(I)
            Call Util_File_Name(Tmp, .TRUE., .FALSE., NameI(I, 2, 1) )
         EndIf
C
         If ( DoGrad ) then
C
C           If we're doing gradients, we'll also be storing the 
C           "3-index MO Gamma", which is the same size as above, but
C           comes in only one "type".
C
            TotSiz = TotSiz + SizeB
            NrDRA = NrDRA + 1
C
            Tmp = BaseNameG(:Inp_StrLen(BaseNameG)) // '-'
     $           // SpinItoA(I)
            Call Util_File_Name(Tmp, .TRUE., .FALSE., NameG(I) )
C
C           We also need to do the (ij|D) integrals.  Also only one
C           "type" but a different size from above.
C
            Size =  (IHi(i)-ILo(i)+1) * (IHi(i)-ILo(i)+1) * NFit
            SizeB = MA_SizeOf(MT_Dbl, Size, MT_Byte)
            MaxSiz = Max( MaxSiz, SizeB)
            TotSiz = TotSiz + SizeB
            NrDRA = NrDRA + 1
C
            Tmp = BaseNameI(:Inp_StrLen(BaseNameI)) // '-i'
     $           // SpinItoA(I)
            Call Util_File_Name(Tmp, .TRUE., .FALSE., NameI(I, 1, 2) )
         EndIf
      EndDo
C
      Return
      End
