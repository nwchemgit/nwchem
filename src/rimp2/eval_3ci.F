C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     Eval_3ci_block -- Produce a block of three-center integrals for
c     the shell block ranges specified.
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine Eval_3ci_block(UseERI, AOBas, PShLo, PShHi, QShLo,
     $   QSHHi, FitBas, SShLo, SShHi, Neglect, Int, LInt, Scr, LScr,
     $   IBuf, PDim, QDim, SDim)
      Implicit NONE
#include "errquit.fh"
      Logical UseERI
      Integer AOBas, PShLo, PShHi, QShLo, QShHi, FitBas, SShLo, SShHi
      Integer LInt, LScr, PDim, QDim, SDim
      Double precision Neglect
      Double precision Int(LInt), Scr(LScr),IBuf(PDim, QDim, SDim)
C
C ARGUMENTS
C     UseERI   Compute 2e3c ERIs (.TRUE.) or 1e3c overlaps [IN]
C     AOBas    Handle for AO basis. [IN]
C     PShLo    Lower limit for shell range of first AO index [IN]
C     PShHi    Upper limit for shell range of first AO index [IN]
C     QShLo    Lower limit for shell range of second AO index [IN]
C     QShHi    Upper limit for shell range of second AO index [IN]
C     FitBas   Handle for fitting basis. [IN]
C     SShLo    Lower limit for shell range of fitting index [IN]
C     SShHi    Upper limit for shell range of fitting index [IN]
C     Neglect  Neglect raw integrals smaller than this value [IN]
C     Int      Target array for integral evaluation [scratch]
C     LInt     Length of Int
C     Scr      Scratch arrya for integral evaluation [scratch]
C     LScr     Length of Scr
C     IBuf     Output integral buffer [OUT]
C     PDim     First dimension of IBuf (first AO index) [IN]
C     QDim     Second dimension of IBuf (second AO index) [IN]
C     SDim     Third dimension of IBuf (fitting index) [IN]
C
C DESCRIPTION
C     Produce a single three-dimensional array containing the three
c     center integrals for a range of shell blocks.  This cannot be done
c     in general by simply putting successive integral batches next to
c     each other in an array because of the three-dimensional nature of
c     the data.  This routine loops over the requested shell blocks and
c     copies the integrals into appropriate locations in the output
c     array.
C
C MEMORY USE
C     All work arrays are provided by parent.
C
C INCLUDE FILES
#include "numerical_constants.fh"
#include "bas.fh"
#include "schwarz.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer SLo, SHi, Qlo, Qhi, Plo, Phi
      Integer P, Q, S, NP, NQ, NS
      Integer SSh, PSh, QSh
      Integer PBase, QBase, SBase, Junk, IntOff
      Logical Status
C
C     **************************
C     * Clear out target array *
C     **************************
C     Must do this so we can safely skip blocks later
C
      Call dlaset('Full', PDim*QDim, SDim, FP_Zero, FP_Zero,
     $     IBuf, PDim*QDim)
C
C     ****************
C     * Do the work! *
C     ****************
C     First, we need to know where we're starting from to get the
c     indexing right in IBuf.
C
      Status = .TRUE.
      Status = Status .AND. Bas_Cn2BFR( AOBas,  PShLo, PBase, Junk)
      Status = Status .AND. Bas_Cn2BFR( AOBas,  QShLo, QBase, Junk)
      Status = Status .AND. Bas_Cn2BFR( FitBas, SShLo, SBase, Junk)
C
      If ( .NOT. Status) Call ErrQuit(
     $   'Eval_3ci_Block: Problems with ''base'' inquiries', PSHLo,
     &       CALC_ERR)
C
      Do PSh = PShLo, PShHi
         Status = Status .AND. Bas_Cn2BFR( AOBas, PSh, Plo, Phi)
         NP = Phi - Plo +1
C
         Do Qsh = QShLo, QShHi
            Status = Status .AND. Bas_Cn2BFR( AOBas, QSh, Qlo, Qhi)
            NQ = Qhi - Qlo +1
C
            If ( Schwarz_Shell( PSh, QSh) * Schwarz_Max()
     $         .gt. Neglect) then
C
               Do SSh = SShlo, SShhi
                  Status = Status
     $               .AND. Bas_Cn2BFR( FitBas, SSh, Slo, Shi)
                  NS = Shi - Slo +1
C
C                 Evaluate the integrals as we've been asked
C
                  If ( UseERI ) then
                     Call Int_2e3c( FitBas, SSh, AOBas, QSh, PSh,
     $                  LScr, Scr, LInt, Int)
                  Else
                     Call Int_1e3ov( FitBas, SSh, AOBas, QSh, AOBAs,
     $                  PSh, LScr, Scr, LInt, Int)
                  EndIf
C
C                 Transfer these integrals to the buffer.  Note that the
c                 integrals are returned to us organized as an array of
c                 dimensions NPxNQxNS, with P running fastest.
C
                  IntOff = 1
                  Do S = 1, NS
                     Do Q = 1, NQ
                        Do P = 1, NP
                           IBuf(P+PLo-PBase, Q+QLo-QBase, S+SLo-SBase)
     $                        = Int( IntOff )
                           IntOff = IntOff + 1
                        EndDo   ! P
                     EndDo      ! Q
                  EndDo         ! S
C
               EndDo            ! SSh
            EndIf               ! Schwarz screening
         EndDo                  ! QSh
      EndDo                     ! PSh
C
      If ( .NOT. Status) Call ErrQuit(
     $   'Eval_3ci_Block: Problems with contraction inquiries', PSHLo,
     &       CALC_ERR)
C
      Return
      End
