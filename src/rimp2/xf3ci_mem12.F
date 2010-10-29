C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     xf3ci_mem12 -- Compute best local memory usage for transformation
c     steps 1 and 2
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine XF3cI_Mem12( MAvail, PMin, QMin, SMin, PMax, QMax,
     $   SMax, NPXF, NQXF, PDim, QDim, SDim, MUsed, LCP, LCQ,
     $   LInterm, LOut, LIBuf)
      Implicit NONE
      Integer MAvail, PMin, QMin, SMin, PMax, QMax,
     $   SMax, NPXF, NQXF, PDim, QDim, SDim, MUsed, LCP, LCQ,
     $   LInterm, LOut, LIBuf
C
C ARGUMENTS
C     MAvail   Local memory available for use in transformation
C              (in double precision words) [IN]
C     PMin     Minimum size allowed for P blocking [IN]
C     QMin     Minimum size allowed for Q blocking [IN]
C     SMin     Minimum size allowed for S blocking [IN]
C     PMax     Maximum size allowed for P blocking [IN]
C     QMax     Maximum size allowed for Q blocking [IN]
C     SMax     Maximum size allowed for S blocking [IN]
C     NPXF     Number of functions P index is transformed into [IN]
C     NQXF     Number of functions Q index is transformed into [IN]
C     PDim     Actual size allowed for P blocking [OUT]
C     QDim     Actual size allowed for Q blocking [OUT]
C     SDim     Actual size allowed for S blocking [OUT]
C     MUsed    Memory required by local arrays [OUT]
C     LCP      Length of CP array (double precision words) [OUT]
C     LCQ      Length of CQ array (double precision words) [OUT]
C     LInterm  Length of Interm array (double precision words) [OUT]
C     LOut     Length of Out array (double precision words) [OUT]
C     LIBuf    Length of IBuf array (double precision words) [OUT]
C
C DESCRIPTION
C     Attempts to compute best use of local memory to improve
C     performance of first two transformation steps.
C
C     The minimum limits for P, Q, and S should be at least the
C     size of the largest shell block in the appropriate basis to 
C     insure things will work.  The maximums would generally be
C     the size of the basis (P, Q) or the size of the local piece
C     of the basis (S) but may be lower to controll memory usage.
C
C     If there is insufficient memory to even meet the minimums, 
C     PDim, QDim and SDim are returned 0 and MUsed holds the ammount
C     of memory that would be required for the minimum size blocks.
C
C     All array lengths are computed here to simplify keeping
C     memory allocation consistent with the computational code.
C     This way all our parent routine has to do is allocate whatever
C     size we tell it.
C
C     Uses a simplistic algorithm right now which simply maximizes
C     Q, then P, then S.  S is last because it always forms a combined
C     index with one of the others, so it is usually large enough.
C     But more emphasis should be placed on the size of P since
C     about 1/3 of the work has PDim as the contraction length.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer MStep
C
C STATEMENT FUNCTIONS
      Integer P, Q, S
      Integer MCP, MCQ, MInterm, MOut, MIBuf, MemReq
      MCP(P,Q,S)     = P * NPXF
      MCQ(P,Q,S)     = Q * NQXF
      MInterm(P,Q,S) = NQXF * P * S
      MOut(P,Q,S)    = NQXF * P * S
      MIBuf(P,Q,S)   = P * Q * S
C
      MemReq(P,Q,S) = MCP(P,Q,S) + MCQ(P,Q,S) + MInterm(P,Q,S)
     $   + MOut(P,Q,S) + MIBuf(P,Q,S)
C
C     Start with the minimums
C
      PDim = PMin
      QDim = QMin
      SDim = SMin
C
C     First, make sure we can do the minimum.  If not, the signal is to
c     set all dimensions to zero.
C
      MUsed = MemReq(PDim, QDim, SDim)
      If ( MUsed .gt. MAvail) then
         PDim = 0
         QDim = 0
         SDim = 0
         Goto 8000
      EndIf
C
C     The following procedure depends on the MemReq function being
c     linear in each of the arguments.
C
C     Make Q as large as we can given the available memory since that is
c     most important to performance.  To do this, take the "derivative"
c     of the memory requirement function w.r.t. Q and use that memory
c     step size to fill in the rest of the memory.
C
      MStep = MemReq(PDim, 1, SDim)
      QDim = Min( QDim + (MAvail - MUsed) / MStep, QMax)
C
C     Now maximize P in the same way.
C
      MUsed = MemReq(PDim, QDim, SDim)
      MStep =  MemReq( 1, QDim, SDim)
      PDim = Min( PDim + (MAvail - MUsed) / MStep, PMax)
C
C     And finally, S...
C
      MUsed = MemReq(PDim, QDim, SDim)
      MStep =  MemReq( PDim, QDim, 1)
      SDim = Min( SDim + (MAvail - MUsed) / MStep, SMax)
C
C     We've got something workable.  Set the sizes and return
C
 8000 Continue
      MUsed   = MemReq(PDim, QDim, SDim)
      LCP     = MCP(PDim, QDim, SDim)
      LCQ     = MCQ(PDim, QDim, SDim)
      LInterm = MInterm(PDim, QDim, SDim)
      LOut    = MOut(PDim, QDim, SDim)
      LIBuf   = MIBuf(PDim, QDim, SDim)
C
      Return
      End
