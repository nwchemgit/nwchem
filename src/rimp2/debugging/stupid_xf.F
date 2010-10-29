      Subroutine Stupid_XF(G_C, G_V, AOBas, FitBas)
C$Id$
      Implicit NONE
C
      Integer G_C, G_V, AOBas, FitBas
C
      Integer NBF, NOcc, NVir, NFit
C     -------------------- biphenyl ------------------------------
      Parameter (NBF = 204, NOcc = 73, NVir = NBF-NOcc, NFit = 204)
c$$$C     -------------------------- He16 ----------------------------
c$$$      Parameter (NBF = 32, NOcc = 16, NVir = NBF-NOcc, NFit = 272)
C
      Double precision AOInts(NBF, NFit, NBF), First(NVir, NFit, NBF),
     $   Sec(NVir, NFit, NOcc), Third(NVir, NFit, NOcc), C(NBF, NBF),
     $   V(NFit, NFit)
C
      Logical Status
      Integer NAOSh, NFitSh, PSh, QSh, SSh, PLo, PHi, QLo, QHi,
     $   SLo, SHi, P, Q, S, Offset, LInt, LScr, H_Int, Int_Idx,
     $   H_Scr, Scr_Idx, I
C
#include "numerical_constants.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
C
C     Step 0: Generate AO integrals
C
      Write (6, *) '%%%%% Starting integral generation'
      Call Util_Flush(6)
C
      Call Int_Mem_2e3c(LInt, LScr)
      Status = .TRUE.
      Status = Status .AND. MA_Push_Get( MT_Dbl, LInt,
     $   'Integral eval', H_Int, Int_Idx)
      Status = Status .AND. MA_Push_Get( MT_Dbl, LScr,
     $   'Integral scr', H_Scr, Scr_Idx)
C
      Status = Status .AND. Bas_NumCont(AOBas, NAOSh)
      Status = Status .AND. Bas_NumCont(FitBas, NFitSh)
C
C
C
      Do Psh = 1, NAOSh
         Status = Bas_Cn2BFR( AOBas, PSh, Plo, Phi)
C
         Do Qsh = 1, NAOSh
            Status = Bas_Cn2BFR( AOBas, QSh, Qlo, Qhi)
C
            Do SSh = 1, NFitSh
               Status = Bas_Cn2BFR( FitBas, SSh, Slo, Shi)
C
C              Note that these guys come out (pq|s), despite the order
c              of the arguments.
C
               Call Int_2e3c( FitBas, SSh, AOBas, PSh, QSh,
     $            LScr, Dbl_MB(Scr_Idx), LInt, Dbl_MB(Int_Idx) )
C
               Offset = 0
               Do S = SLo, SHi
                  Do Q = QLo, QHi
                     Do P = PLo, PHi
                        AOInts(P, S, Q) = Dbl_MB(Int_Idx + Offset)
                        Offset = Offset + 1
                     EndDo
                  EndDo
               EndDo
C
            EndDo
         EndDo
      EndDo
C
      Status = MA_Pop_Stack(H_Scr)
      Status = MA_Pop_Stack(H_Int)
C
C     Step 1: Transform first index
C
      Write (6, *) '%%%%% Starting first transformation'
      Call Util_Flush(6)
C
      Call GA_Get(G_C, 1, NBF, 1, NBF, C, NBF)
      Call SGEMM('T', 'N', NVir, NFit*NBF, NBF,
     $  FP_One, C(1, NOcc+1), NBF, AOInts, NBF,
     $   FP_Zero, First, NVir)
C
C     Step 2: Transform second index
C
      Write (6, *) '%%%%% Starting second transformation'
      Call Util_Flush(6)
C
      Call SGEMM('N', 'N', NVir*NFit, NOcc, NBF,
     $   FP_One, First, NVir*NFit, C(1, 1), NBF, 
     $   FP_Zero, Sec, NVir*NFit)
C
C     Step 2a: write out partially transformed integrals
C
      Write (6, *) '%%%%% Starting first output'
      Call Util_Flush(6)
C
      Call OpenRec(50, NVir*NFit)
      Do I = 1, NOcc
         Call PutRec(50, I, NVir*NFit, Sec(1, 1, I) )
      EndDo
C
C     Step 3: Transform fitting index
C
      Write (6, *) '%%%%% Starting third transformation'
      Call Util_Flush(6)
C
      Call GA_Get(G_V, 1, NFit, 1, NFit, V, NFit)
      Do I = 1, NOcc
         Call SGEMM('N', 'N', NVir, NFit, NFit,
     $      FP_One, Sec(1, 1, I), NVir, V, NFit,
     $      FP_Zero, Third(1, 1, I), NVir)
      EndDo
C
C     Step 3a: write out fully transformed integrals
C
      Write (6, *) '%%%%% Starting second output'
      Call Util_Flush(6)
C
      Call OpenRec(51, NVir*NFit)
      Do I = 1, NOcc
         Call PutRec(51, I, NVir*NFit, Third(1, 1, I) )
      EndDo
C
C
C
      Return
      End
C
C
C
      Subroutine OpenRec(LU, NElem)
      Implicit NONE
      Integer LU, NElem
C
      Open (UNIT=LU, ACCESS='DIRECT', RECL=NElem*8)
      Return
      End
C
C
C
      Subroutine PutRec(LU, RecNr, NElem, A)
      Implicit NONE
      Integer LU, RecNr, NElem
      Double Precision A(NElem)
C
      Integer Junk
      Double precision ChkSum
C
      External SASUM, ISAMAX
      Double precision SASUM
      Integer IASMAX
C
      Write (LU, REC=RecNr) A
C
      ChkSum = SASUM(NElem, A, 1)
      Junk = ISAMAX(NElem, A, 1)
      Write (6, 9000) RecNr, ChkSum, A(Junk)
 9000 Format(1X, 'Record ', I5, ' checksum = ', F30.15,
     $   ' max = ', F30.15)
      Return
      End
