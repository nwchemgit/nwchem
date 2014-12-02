C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     XF3cI_Step12b -- Perform the first 1/3 transformation,
c     accumulating the result for a complete pass through the innermost
c     index.
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine XF3cI_Step12b(UseERI, AOBas, PShLo, PShHi, PDim,
     $   QDim, FitBas, SShLo, SShHi, SDim, Neglect,
     $   Int, LInt, Scr, LScr, IBuf, CQ, QXFLo, QXFHi, G_C,
     $   Interm, Out, LDOut, FLOPs)
      Implicit NONE
#include "errquit.fh"
      Logical UseERI
      Integer AOBas, PShLo, PShHi, PDim, QDim, SDim
      Integer FitBas, SShLo, SShHi, LInt, LScr
      Integer LDOut, QXFLo, QXFHi, G_C
      Double precision Neglect, FLOPs
      Double precision Int(LInt), Scr(LScr), CQ(QDim, QXFHi-QXFLo+1),
     $   IBuf(QDim, PDim, SDim), Interm(QXFHi-QXFLo+1, PDim, SDim),
     $   Out(LDOut, QXFHi-QXFLo+1, SDim)
C
C ARGUMENTS
C     UseERI   Compute 2e3c ERIs (.TRUE.) or 1e3c overlaps [IN]
C     AOBas    Handle for AO basis. [IN]
C     PShLo    Lower limit of shell range of P index [IN]
C     PShHi    Upper limit of shell range of P index [IN]
C     PDim     Space allocated to P index of IBuf (must be at least the
C              number of functions in the range PShLo:PShHi) [IN]
C     QDim     Space allocated to Q index of IBuf
C     FitBas   Handle for fitting basis. [IN]
C     SShLo    Lower limit of shell range of fitting (S) index [IN]
C     SShHi    Upper limit of shell range of fitting (S) index [IN]
C     SDim     Space allocated to S index of IBuf (must be at least the
C              number of functions in the range SShLo:SShHi) [IN]
C     Neglect  Neglect raw integrals smaller than this value [IN]
C     Int      Integral evaluation array for Eval_3ci_Block [scratch]
C     LInt     Length of Int
C     Scr      Integral evaluation scratch array for Eval_3ci_Block [scr]
C     LScr     Length of Scr
C     IBuf     Target integral buffer for Eval_3ci_Block [scratch]
C     CQ       Space for Q index transformation matrix [scratch]
C     QXFLo    Lower limit of MO functions to transform Q into [IN]
C     QXFHi    Upper limit of MO functions to transform Q into [IN]
C     G_C      GA handle for Q index transformation matrix [IN]
C     Interm   Used to hold 1/3 transformed integrals before transposition
C              to (pi|s) organization [scratch]
C     Out      Matrix of 1/3 transformed integrals (pi|s) [OUT]
C     LDOut    Leading dimension of Out [IN]
C     FLOPs    Number of floating point operations performed [INOUT]
C
C DESCRIPTION
C     Performs the first index transformation for three-center integrals,
C     making a complete pass through the Q index for a specified range of
C     P and S functions.  The output is a matrix of integrals (pi|s),
C     where q --> i is the transformation done here.
C
C     The Q index is covered in multiple passes according to the QDim
C     provided by the caller and respecting shell boundaries.  The size
C     QDim governs the contraction length in a matrix multiplication
C     so ideally it should be made large enough to obtain asymptotic
C     performance from the matrix multiplication.
C
C     The structure of the integrals that is to be output would not allow
C     a single large DGEMM in general, so this routine evaluates the 
C     integrals in such a way as to produce the large single DGEMM
C     desired, then transposes the result into the form expected by the
C     parent.  This will slow things a little compared to asymptotic 
C     DGEMM speeds, but unless QDim is very small, it is worth it.

C MEMORY USE
C     All work arrays are provided by parent.
C
C INCLUDE FILES
#include "numerical_constants.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "global.fh"
#include "rimp2_timers.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer P, Q, S, QShLo, QShHi, QBFLo, QBFHi, NQXF, Low, High, MaxQ
      Integer PCenLo, PCenHi, QCenLo, QCenHi
      Logical Status, MustDo
C
      NQXF = QXFHi - QXFLo + 1
      Status = Bas_NumBF( AOBas, MaxQ)
C
      Call dlaset( 'Full', NQXF, PDim*SDim, FP_Zero, FP_Zero,
     $     Interm, NQXF)
C
C     This is the starting place for our romp through the Q space
C
      QBFLo = 1
C
C     This is where we want end up -- careful not to run off the end of
c     the basis, since the inquiry routines do not appreciate that.
C
 1000 QBFHi = Min( QBFLo + QDim - 1, MaxQ)
C
C     Find out what contractions these belong to
C
      Status = Bas_BF2Cn( AOBas, QBFLo, QShLo)
      Status = Bas_BF2Cn( AOBas, QBFHi, QShHi)
C
C     We must be at the beginning of QShLo, so reset QBFLo to guarantee
c     this.
C
      Status = Bas_Cn2BFR( AOBas, QShLo, QBFLo, High)
C
C     And we need to be at the first shell boundary below QBFHi.  This
c     may necessitate backing off one shell if QBFHi is not already on a
c     boundary.
C
      Status = Bas_Cn2BFR( AOBas, QShHi, Low, High)
      If ( High .gt. QBFHi) then
         QShHi = QShHi - 1
         Status = Bas_Cn2BFR( AOBas, QShHi, Low, QBFHi)
      EndIf
C
C     Now (finally) we have a good range to use: QShLo:QShHi, or
c     QBFLo:QBFHi for shells or functions!
C
C     ***************************************
C     * See if we can skip the entire batch *
C     ***************************************
C     Use the atom-level Schwarz info to check that we have work to do
c     at this level.
C
      Status = .TRUE.
      Status = Status .AND. Bas_Cn2Ce( AOBas,  PShLo, PCenLo)
      Status = Status .AND. Bas_Cn2Ce( AOBas,  PShHi, PCenHi)
      Status = Status .AND. Bas_Cn2Ce( AOBas,  QShLo, QCenLo)
      Status = Status .AND. Bas_Cn2Ce( AOBas,  QShHi, QCenHi)
C
      If ( .NOT. Status) Call ErrQuit(
     $   'XF3cI_Step12b: Problems with center inquiries', PSHLo,
     &       BASIS_ERR)
C
      MustDo = .FALSE.
      Do P = PCenLo, PCenHi
         Do Q = QCenLo, QCenHi
            MustDo = MustDo
     $         .OR. (Schwarz_Atom( P, Q) * Schwarz_Max() .gt. Neglect)
         EndDo                  ! Q
      EndDo                     ! P
C
      If ( MustDo ) then
C        
C        ************************************************
C        * Grab the transformation matrix for this pass *
C        ************************************************
C
         Call GA_Get(G_C, QBFlo, QBFhi, QXFlo, QXFhi, CQ, QDim)
C
         Call PStat_On( PS_Int_Gen)
C
         Call Eval_3ci_Block( UseERI, AOBas, QShLo, QShHi, PShLo, PShHi,
     $      FitBas, SShLo, SShHi, Neglect, Int, LInt, Scr, LScr, IBuf,
     $      QDim, PDim, SDim)
C
         Call PStat_Off( PS_Int_Gen)
         Call PStat_On( PS_Step12)
C
C        Transform this block: C_{qi}^T (qp|s) --> (ip|s)
C        
         FLOPs = FLOps + 2 * NQXF * PDim*SDim * (QBFHi-QBFLo+1)
         Call dgemm('T', 'N', NQXF, PDim*SDim, QBFHi-QBFLo+1,
     $      FP_One, CQ, QDim, IBuf, QDim,
     $      FP_One, Interm, NQXF)
C
      Call PStat_Off( PS_Step12)
C     
      EndIf                     ! Schwarz test
C        
C     Pick the next starting point
C
      If ( QBFHi .lt. MaxQ) then
         QBFLo = QBFHi+1
         Goto 1000
      EndIf
C
C     (ip|s) --> (pi|s) (necessary to be able to do a matrix
c     multiplication for the second index too, unfortunately)
C
      Call PStat_On( PS_Step12)
C     
      Do S = 1, SDim
         Call DGECop('T', PDim, NQXF, Interm(1, 1, S), NQXF,
     $      Out(1, 1, S), LDOut)
      EndDo
C     
      Call PStat_Off( PS_Step12)
C
      Return
      End

