C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     XF3cI_Step12a -- Do the first two steps of the three-center
C     integral transformation
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine XF3cI_Step12a(UseERI, AOBas, PDim, QDim, FitBas,
     $   SDim, Neglect, Int, LInt, Scr, LScr, IBuf,
     $   CP, PXFLo, PXFHi, CQ, QXFLo, QXFHi, G_C,
     $   Interm, Out, LDOut, Local, LDLocal, NodeSLo, NodeSHi, FLOPs)
      Implicit NONE
#include "errquit.fh"
      Logical UseERI
      Integer AOBas, PDim, QDim, SDim
      Integer FitBas, LInt, LScr
      Integer LDOut, PXFLo, PXFHi, QXFLo, QXFHi, G_C
      Integer LDLocal, NodeSLo, NodeSHi
      Double precision Neglect, FLOPs
      Double precision CP(PDim, PXFHi-PXFLo+1), CQ(QDim, QXFHi-QXFLo+1),
     $   Int(LInt), Scr(LScr),
     $   Out(LDOut, QXFHi-QXFLo+1, SDim),
     $   IBuf(QDim, PDim, SDim),
     $   Interm(QXFHi-QXFLo+1, PDim, SDim),
     $   Local(LDLocal, *)
C
C ARGUMENTS
C     UseERI   Compute 2e3c ERIs (.TRUE.) or 1e3c overlaps [IN]
C     AOBas    Handle for AO basis. [IN]
C     PDim     Space allocated to P index of IBuf [IN]
C     QDim     Space allocated to Q index of IBuf [IN]
C     FitBas   Handle for fitting basis. [IN]
C     SDim     Space allocated to S index of IBuf [IN]
C     Neglect  Neglect raw integrals smaller than this value [IN]
C     Int      Integral evaluation array for Eval_3ci_Block [scratch]
C     LInt     Length of Int
C     Scr      Integral evaluation scratch array for Eval_3ci_Block [scr]
C     LScr     Length of Scr
C     IBuf     Target integral buffer for Eval_3ci_Block [scratch]
C     CP       Space for P index transformation matrix [scratch]
C     PXFLo    Lower limit of MO functions to transform P into [IN]
C     PXFHi    Upper limit of MO functions to transform P into [IN]
C     CQ       Space for Q index transformation matrix [scratch]
C     QXFLo    Lower limit of MO functions to transform Q into [IN]
C     QXFHi    Upper limit of MO functions to transform Q into [IN]
C     G_C      GA handle for transformation matrix [IN]
C     Interm   Used to hold 1/3 transformed integrals before transposition
C              to (pi|s) organization [scratch]
C     Out      Matrix of 1/3 transformed integrals (pi|s) [OUT]
C     LDOut    Leading dimension of Out [IN]
C     Local    Output array containing 2/3 transformed integrals [INOUT]
C     LDLocal  Leading dimension of Local [IN]
C     NodeSLo  Lower limit of this node's range of the fitting basis [IN]
C     NodeSLo  Upper limit of this node's range of the fitting basis [IN]
C     FLOPs    Number of floating point operations performed [INOUT]
C
C DESCRIPTION
C     Performs first two index transformations of the three-center
C     integrals.  The output is a matrix of integrals (ai|s) where
C     the fitting basis s is distributed across the machine. Perm-
C     utational symmetry of the AO basis indices is NOT used in this
C     algorithm.
C
C     Each index is covered in multiple passes governed by the PDim,
C     QDim, and SDim provided by the caller.  Sheel boundaries are
C     respected so it is not possible to translate the dimensions
C     into a number of passes without considering the shell blocking
C     of the basis.
C
C INCLUDE FILES
#include "global.fh"
#include "numerical_constants.fh"
#include "bas.fh"
#include "rimp2_timers.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Integer PBFLo, PBFHi, PShLo, PShHi,
     $   NPXF, NQXF, SBFLo, SBFHi, SShLo, SShHi, MaxP, High, Low, NP,
     $   NS
      Logical Status
C
      NPXF = PXFHi-PXFLo+1
      NQXF = QXFHi-QXFLo+1
      Status = .TRUE.
      Status = Status .AND. Bas_NumBF( AOBas, MaxP)
C
      If ( .NOT. Status) Call ErrQuit(
     $   'XF3cI_Step12a: Problems with fitting basis inquiries', AOBas,
     &       BASIS_ERR)
C
      If ( LDLocal .ne. NPXF*NQXF) Call ErrQuit(
     $   'XF3cI_Step12a: Dimensions of local patch incorrect', LDLocal,
     &       BASIS_ERR)
C
C     *******************
C     * Start of P loop *
C     *******************
C     Set the range for the next block of P to be done.  We have been
c     given a fixed overall size (number of basis functions) we can
c     accomodate in PDim.  We want to use as much of that as possible
c     insuring that we end up on a shell block boundary. Procedure:
C     1) Guess end of block
C     2) Find out what shells the beginning and ending functions are
C     3) Insure that we start and end on shell block boundaries at most
c     PDim functions apart.
C
      PBFLo = 1
C
 1000 PBFHi = Min( PBFLo + PDim - 1, MaxP)
C
      Status = Bas_BF2Cn( AOBas, PBFLo, PShLo)
      Status = Bas_BF2Cn( AOBas, PBFHi, PShHi)
C
      Status = Bas_Cn2BFR( AOBas, PShLo, PBFLo, High)
C
      Status = Bas_Cn2BFR( AOBas, PShHi, Low, High)
      If ( High .gt. PBFHi) then
         PShHi = PShHi - 1
         Status = Bas_Cn2BFR( AOBas, PShHi, Low, PBFHi)
      EndIf
C
      NP = PBFHi - PBFLo + 1
C
C     **********************************************************
C     * Get the appropriate patch of the transformation matrix *
C     **********************************************************
C
      Call GA_Get(G_C, PBFlo, PBFhi, PXFLo, PXFhi,
     $   CP, PDim)
C
C        *******************
C        * Start of S loop *
C        *******************
C        Same procedure as for S
C
         SBFLo = NodeSLo
C
 1100    SBFHi = Min( SBFLo + SDim - 1, NodeSHi)
C
         Status = Bas_BF2Cn( FitBas, SBFLo, SShLo)
         Status = Bas_BF2Cn( FitBas, SBFHi, SShHi)
C
         Status = Bas_Cn2BFR( FitBas, SShLo, SBFLo, High)
C
         Status = Bas_Cn2BFR( FitBas, SShHi, Low, High)
         If ( High .gt. SBFHi) then
            SShHi = SShHi - 1
            Status = Bas_Cn2BFR( FitBas, SShHi, Low, SBFHi)
         EndIf
C
         NS = SBFHi - SBFLo + 1
C
C           *************************************************
C           * Accumulate a 1/3 transformed result for all Q *
C           *************************************************
C
         Call dlaset('Full', LDOut, NQXF*SDim, FP_Zero, FP_Zero,
     $        Out, LDOut)
C
            Call XF3cI_Step12b(UseERI, AOBas, PShLo, PShHi, NP,
     $         QDim, FitBas, SShLo, SShhi, SDim, Neglect,
     $         Int, LInt, Scr, LScr, IBuf,
     $         CQ, QXFLo, QXFHi, G_C,
     $         Interm, Out, LDOut, FLOPs )
C
C           ******************************
c           * C_{pa}^T (pi|s) --> (ai|s) *
C           ******************************
C           Notice that we are thinking of the result we generate
C           as an NPXF by NQXF*NS matrix, while in reality, it is
C           an NPXF*NQXF by NS matrix.  Therefore, of the leading
c           dimension is not exactly NPXF*NQXF, we are in trouble!  
c           This is checked when we get access to the local patch,
c           above.
C
            Call PStat_On( PS_Step12)
C
            FLOPs = FLOPs + 2 * NPXF * NQXF*NS * NP
            Call dgemm('T', 'N', NPXF, NQXF*NS, NP,
     $         FP_One, CP(1, 1), PDim, Out, LDOut,
     $         FP_One, Local(1, SBFLo-NodeSLo+1), NPXF)
C
            Call PStat_Off( PS_Step12)
C        
C        **********************
C        * Continue the loops *
C        **********************
C
         If ( SBFHi .lt. NodeSHi) then
            SBFLo = SBFHi+1
            Goto 1100
         EndIf
C
      If ( PBFHi .lt. MaxP) then
         PBFLo = PBFHi+1
         Goto 1000
      EndIf
C
      Return
      End




