C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_prep_dra -- Setup DRA sizes and names
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_Prep_DRA(MaxSpin, TopSpin, NFrzO, NAct, NVir,
     $     NFrzV, NFit, DoGrad, BraKetSame, NrDRA, MaxSiz, TotSiz,
     $     BaseNameI, BaseNameG, NameI, NameG, IntRange, GamRange)
      Implicit NONE
#include "errquit.fh"
C
      Integer MaxSpin ![in]
      Integer TopSpin ![in]
      Integer NFrzO(TopSpin) ![in]
      Integer NAct(TopSpin) ![in]
      Integer NVir(TopSpin) ![in]
      Integer NFrzV(TopSpin) ![in]
      Integer NFit ![in]
      Logical DoGrad ![in]
      Logical BraKetSame ![in]
      Integer NrDRA ![out]
      Double precision MaxSiz ![out]
      Double precision TotSiz ![out]
      Character*(*) BaseNameI ![in]
      Character*(*) BaseNameG ![in]
      Character*(*) NameI(MaxSpin, 2, 2) ![out]
      Character*(*) NameG(TopSpin) ![out]
      Integer IntRange(2, 2, MaxSpin, 2, 2) ![out]
      Integer GamRange(2, 2, MaxSpin) ![out]
C
C DESCRIPTION
C     Prepare names, sizes, and transformation limits for DRAs
C
C                    Value of BraKetSame/DoGrad
C     D_Int()     T/F       F/F       T/T       F/T
C
C     *, 1, 1   (ai|D')   (ai|D)    (ap|D')   (ap|D) 
C     *, 2, 1      -      (ai|D")      -      (ap|D")    
C     *, 1, 2      -         -      (mi|D')   (mi|D) 
C     *, 2, 2      -         -         -      (mi|D")    
C
C     D_Gam()   Gamma
C     *         (ai|D)
C
C     Key:
C     D  = fitting basis, untransformed
C     D' = fitting basis with symmetric transformation
C     D" = fitting basis with asymmetric transformation
C     i  = A (active)
C     a  = V (virtual)
C     m  = C+A (core + active)
C     p  = C+A+V
C
C EXTERNAL ROUTINES
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer MyMaxSpin
      Parameter (MyMaxSpin = 2)
C
C LOCAL VARIABLES
      Double precision Size, SizeB
      Integer NrTypes, I, J, K, C, A, V, N1, N2
      Integer ISpin, IType
      Character*80 Tmp
      Logical PrFilInfo
C
C STATEMENT FUNCTIONS
      Character*1 SpinItoA
      SpinItoA(ISpin) = Char( ISpin + IChar('a') - 1 )
      Character*1 ItoA
      ItoA(I) = Char( I + IChar('0') )
C
      If ( TopSpin .gt. MyMaxSpin ) Call ErrQuit(
     $     'RIMP2_Prep_DRA: fatal program error: TopSpin > MyMaxSpin',
     $     MyMaxSpin, INPUT_ERR)
C
C     Initialize transformation ranges to 0:-1 (zero trip loop)
C
      Call IFill(1*2*MaxSpin*2*2,  0, IntRange(1, 1, 1, 1, 1), 2)
      Call IFill(1*2*MaxSpin*2*2, -1, IntRange(2, 1, 1, 1, 1), 2)
      Call IFill(1*2*MaxSpin,      0, GamRange(1, 1, 1), 2)
      Call IFill(1*2*MaxSpin,     -1, GamRange(2, 1, 1), 2)
C
      PrFilInfo = Util_Print('file details', Print_Debug)
     $     .AND. (GA_NodeID() .eq. 0)
C
      If ( PrFilInfo ) Write (LuOut, 9000)
C
 9000 Format(/1X, 'Disk File Information',
     $     /1X, 'Spin', 1X, 'Bra/Ket', 1X, 'Class', 13X, 'Size (B)',
     $     3X, 'Name')
 9010 Format(1X, I4, 1X, I7, 1X, I5, 1X, F20.0, 3X, A)
 9020 Format(1X, I4, 7X,  A, 5X,  A, 1X, F20.0, 3X, A)
 9030 Format(1X, I4, 7X,  A, 1X, I5, 1X, F20.0, 3X, A)
 9040 Format(1X, 4X, 1X, 7X, 1X, 5X, 1X, F20.0, 3X, A)
C
C     Drive everything off of the transformation ranges
C
C     First two indices correspond to lower/upper limit of range
C     for index 1/2 of the 3c-integrals.  Remaining indices are
C     exactly like D_Int/D_Gam (spin, bra/ket, class)
C
      Do ISpin = 1, TopSpin
         C = NFrzO(ISpin)
         A = NAct(ISpin)
         V = NVir(ISpin)
C
         If ( DoGrad) then
C
C           (ap|D)
C
            IntRange(1, 1, ISpin, 1, 1) = C+A+1
            IntRange(2, 1, ISpin, 1, 1) = C+A+V
C
            IntRange(1, 2, ISpin, 1, 1) =     1
            IntRange(2, 2, ISpin, 1, 1) = C+A+V
C
C           (mi|D)
C
            IntRange(1, 1, ISpin, 1, 2) =   1
            IntRange(2, 1, ISpin, 1, 2) = C+A
C
            IntRange(1, 2, ISpin, 1, 2) = C+1
            IntRange(2, 2, ISpin, 1, 2) = C+A
C
C           Gamma: (ai|D)
C
            GamRange(1, 1, ISpin) = C+A+1
            GamRange(2, 1, ISpin) = C+A+V
C
            GamRange(1, 2, ISpin) = C+1
            GamRange(2, 2, ISpin) = C+A
         Else                   ! DoGrad
C
C           (ai|D)
C
            IntRange(1, 1, ISpin, 1, 1) = C+A+1
            IntRange(2, 1, ISpin, 1, 1) = C+A+V
C
            IntRange(1, 2, ISpin, 1, 1) = C+1
            IntRange(2, 2, ISpin, 1, 1) = C+A
         EndIf                  ! DoGrad
C
C        Integrals may not be split symmetrically.  Sizes & ranges
c        are identical, of course.  We _don't_ need the ket integrals
C        for class=2
C
         If ( .NOT. BraKetSame ) then
            Do J = 1, 2
               Do K = 1, 2
                  IntRange(J, K, ISpin, 2, 1)
     $                 = IntRange(J, K, ISpin, 1, 1)
               EndDo            ! K
            EndDo               ! J
         EndIf                  ! .NOT. BraKetSame
C
      EndDo                     ! ISpin
C
C     Generate sizes and names
C
      NrDRA = 0
      MaxSiz = 0
      TotSiz = 0
      Do ISpin = 1, TopSpin
         Do J = 1, 2            ! Bra/Ket
            Do K = 1, 2         ! Class
               N1 = IntRange(2, 1, ISpin, J, K)
     $              - IntRange(1, 1, ISpin, J, K) + 1
               N2 = IntRange(2, 2, ISpin, J, K)
     $              - IntRange(1, 2, ISpin, J, K) + 1
               Size = Dble(N1) * Dble(N2) * Dble(NFit)
               SizeB = MA_SizeOf(MT_Dbl, 1, MT_Byte) * Size
C
               If ( Size .gt. 0 ) then
                  MaxSiz = Max( MaxSiz, SizeB)
                  TotSiz = TotSiz + SizeB
                  NrDRA = NrDRA + 1
C
                  Tmp = BaseNameI(:Inp_StrLen(BaseNameI))
     $                 // SpinItoA(ISpin) // ItoA(J) // ItoA(K)
                  Call Util_File_Name(Tmp, .TRUE., .FALSE.,
     $                 NameI(ISpin, J, K) )
C
                  If (PrFilInfo) Write (LuOut, 9010) ISpin, J, K, SizeB,
     $               NameI(ISpin, J, K)(:Inp_StrLen(NameI(ISpin, J, K)))
               Else
                  NameI(ISpin, J, K) = ' '
               EndIf
            EndDo               ! Class (K)
         EndDo                  ! Bra/Ket (J)
C
C        Handle Gammas too, which have only spin cases
C
         N1 = GamRange(2, 1, ISpin)
     $        - GamRange(1, 1, ISpin) + 1
         N2 = GamRange(2, 2, ISpin)
     $        - GamRange(1, 2, ISpin) + 1
         Size = Dble(N1) * Dble(N2) * Dble(NFit)
         SizeB = MA_SizeOf(MT_Dbl, 1, MT_Byte) * Size
C
         If ( Size .gt. 0 ) then
            MaxSiz = Max( MaxSiz, SizeB)
            TotSiz = TotSiz + SizeB
            NrDRA = NrDRA + 1
C
            Tmp = BaseNameG(:Inp_StrLen(BaseNameG))
     $           // SpinItoA(ISpin)
            Call Util_File_Name(Tmp, .TRUE., .FALSE., NameG(ISpin) )
C
            If (PrFilInfo) Write (LuOut, 9020) ISpin, '-', '-', SizeB,
     $           NameG(ISpin)(:Inp_StrLen(NameG(ISpin)))
         Else
            NameG(ISpin) = ' '
         EndIf
      EndDo                     ! ISpin
C
      If ( PrFilInfo ) then
         Write (LuOut, 9040) MaxSiz, 'Maximum File Size'
         Write (LuOut, 9040) TotSiz, 'Total File Size'
      EndIf
      Call Util_Flush(LuOut)
C
      Return
      End
