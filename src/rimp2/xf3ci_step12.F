C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     XF3cI_Step12 -- Do the first two steps of the three-center
C     integral transformation
C
C REVISION
C     $Id: xf3ci_step12.F,v 1.1 1994-06-14 21:54:37 gg502 Exp $
C
C SYNOPSIS
      Logical Function XF3cI_Step12(IntTyp, Alo, Ahi, Ilo, Ihi,
     $   AOBas, G_C, FitBas, G_I)
      Implicit NONE
      Character*1 IntTyp
      Integer Alo, Ahi, Ilo, Ihi, FitBas, AOBas, G_I, G_C
C
C ARGUMENTS
C     IntTyp   Type of integral to be transformed. 'E' for three-center
C              two-electron ERI, (pq|s), or 'O' for three-center
C              one-electron overlap, (pqs). [IN]
C     Alo      Lower limit of A to be generated. [IN]
C     Ahi      Upper limit of A to be generated. [IN]
C     Ilo      Lower limit of I to be generated. [IN]
C     Ihi      Upper limit of I to be generated. [IN]
C     AOBas    Handle for AO basis. [IN]
C     G_C      GA handle for SCF eigenvectors. [IN]
C     FitBas   Handle for fitting basis. [IN]
C     G_I      Handle for global array to hold (ia|s). [IN]
C
C RETURN VALUE
C     .TRUE.   Successful completion
C     .FALSE.  Error condition.  Could be in input, basis info queries,
C              or memory operations (MA or GA).
C DESCRIPTION
C     Produce a block of partially transformed three-center integrals
C     in the caller-provided array G_I.  
C
C     Either 2-electron ERIs or 1-center overlaps may be used by setting
C     IntTyp appropriately.  The result is (Alo:Ahi, Ilo:Ihi | S) or
C     ( Alo:Ahi, Ilo:Ihi, S) according to the type selected.  The global
C     array holds the combined index IA in the first dimension, which is
C     assumed to be undistributed, and S in the second dimension, which
C     is distributed.
C
C     This routine accesses C( : , Alo:Ahi) and C( : , Ilo:Ihi), so
C     as long as everything matches, you can do whatever you want
C     w.r.t. dropping MOs, etc.  If you have multiple spin cases,
C     you simply have to feed the appropriate C in to get what you
C     want.  Three-center integrals approximating four-center two-
C     electron ERIs do not require mixing spins within the three-center
C     pieces.
C
C NOTES
C     In order to minimize caller overhead (and therefore possibilities
C     for errors), we try to make this as self-contained as possible.
C     Consequently, we may end up doing more MA requests or inquiries
C     than is strictly necessary.  If this appears to be a performance
C     issue, this routine could easily be made lean and mean by insisting
C     that the caller provide the two scratch arrays, and some of the
C     basis set info.  But please don't do this unless its necessary!
C
C MEMORY USE
C     Two local arrays are used for multiple purposes:
C     1) Integral value and scratch arrays as sized by Int_Mem_{2e3c,3ov}
C     2) Scratch array must be at least as large as the value array
C        (This will normally be the case just from Imt_Mem_*)
C     3) The integral value array must also be able to hold 
C        (Alo:Ahi, MxBFSh), where MxBFSh is the largest shell in the
C        AOBas basis.
C     Two additional arrays are used to localize panels of C:
C     1) NA * MxBFSh
C     2) NI * MxBFSh
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.h"
#include "numerical_constants.fh"
#include "stdio.fh"
#include "bas.fh"
C
C EXTERNAL ROUTINES
      Logical LSame      ! From BLAS-[23]
      External LSame
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer IAlo, IAhi, Slo, Shi, Plo, Phi, Qlo, Qhi, NI, NA
      Integer NAOSh, MxBFSh, LInt, LScr, MySLo, MySHi
      Integer S, SSh, SShlo, SShhi, PSh, QSh, NBFSSh, NBFPSh, NBFQSh
      Logical Status, UseERI
      Integer Local_Idx, LDLocal, H_Int, H_Scr, Int_Idx, Scr_Idx
      Integer H_CP, H_CQ, CP_Idx, CQ_Idx, LCP, LCQ
      Integer Junk
C
C     Dimensions of A and I ranges
C
      NA = Ahi - Alo + 1
      NI = IHi - ILo + 1
C
C     *****************
C     * Sanity checks *
C     *****************
C     We can handle either 2-electron 3-center ERIs or 1-electron 3-center
C     overlaps here.
C
      UseERI = LSame( IntTyp, 'E')
      If ( .NOT. ( UseERI .or. LSame( IntTyp, 'O') ) )then
         Write (LuErr, *) 'XF3cI_Step12: Unknown integral type ',
     $      IntTyp
         XF3cI_Step12 = .FALSE.
         Return
      EndIf
C
C     Verify that we've got good GA & basis set handles
C
      Call GA_Check_Handle ( G_I, 'XF3cI_Step12: Integral matrix')
      Call GA_Check_Handle ( G_C, 'XF3cI_Step12: SCF eigenvectors')
C
      Status = .TRUE.
      Status = Status .AND.
     $   Bas_Check_Handle( AOBas,  'XF3cI_Step12: AO Basis')
      Status = Status .AND.
     $   Bas_Check_Handle( FitBas, 'XF3cI_Step12: Fitting basis')
C
C     Message already printed, so no need for one here
C
      If (.NOT. Status ) then
         XF3cI_Step12 = .FALSE.
         Return
      EndIf
C
C     Zero out target array
C
      Call GA_Zero( G_I )
C
C     ********************************************
C     * Determine our range of the fitting basis *
C     ********************************************
C     Each node will fill the range that lives in its memory.  So we
C     get the basis functions that define the local panel and translate
C     them to shell blocks.
C
      Call GA_Distribution( G_I, GA_NodeID(), IAlo, IAhi, MySLo, MySHi)
c$$$      Write (LuOut, *) 'Step12 distribution ', IALo,IAhi, MySLo, mySHi
C
      Status = .TRUE.
      Status = Status .AND. Bas_BF2Cn( FitBas, MySLo, SShlo)
      Status = Status .AND. Bas_BF2Cn( FitBas, MySHi, SShhi)
C
C     We also need some info about the AO basis
C
      Status = Status .AND. Bas_NumCont( AOBas, NAOSh)
      Status = Status .AND. Bas_NBF_Cn_Max( AOBas, MxBFSh)
C
      If (.NOT. Status ) then
         Write (LuErr, *) 'XF3cI_Step12: Problems with basis set ',
     $      'queries.'
         XF3cI_Step12 = .FALSE.
         Return
      EndIf
C
C     **********************************
C     * Allocate local "scratch" space *
C     **********************************
C     These guys are both fed into the integral package, and reused
C     during the transformation, so we must be sure to size for both
C     uses.
C
C     This is what the integral package needs
C
      If ( UseERI ) then
         Call Int_Mem_2e3c( LInt, LScr)
      Else
         Call Int_Mem_3ov( LInt, LScr)
      EndIf
C
C     We also need to transpose the resulting integrals from Int to Scr
C     so we have them in a convenient organization.  So we must make
C     sure that Scr is at least the size of Int (though the integrals
C     will probably already require this).
C
      LScr = Max( LScr, LInt )
C
C     And we reuse Int as the intermediate in the two transformation
C     steps, where we need it to be large enough to hold the
C     (aq|s) block for all a, a shell of q, and a single s.
C
      LInt = Max( LInt, NA*MxBFSh)
C
      Status = .TRUE.
      Status = Status .AND. MA_Push_Get( MT_Dbl, LInt, 'Integral eval',
     $   H_Int, Int_Idx)
      Status = Status .AND. MA_Push_Get( MT_Dbl, LScr, 'Integral scr',
     $   H_Scr, Scr_Idx)
C
C     We also need to localize bits of the transformation matrix, C
C
      LCP = NA * MxBFSh
      LCQ = NI * MxBFSh
C
      Status = Status .AND. MA_Push_Get( MT_Dbl, LCP, 'trans. scr',
     $   H_CP, CP_Idx)
      Status = Status .AND. MA_Push_Get( MT_Dbl, LCQ, 'trans. scr',
     $   H_CQ, CQ_Idx)
C
      If (.NOT. Status ) then
         Write (LuErr, *) 'XF3cI_Step12: Problem allocating local ',
     $      'memory.'
         XF3cI_Step12 = .FALSE.
         Return
      EndIf
C
C     ********************************************
C     * Grab our panel of G_I to access directly *
C     ********************************************
C     This avoids the need for a scratch matrix of size NA*NI
C     and a copy operation.
C
      Call GA_Access( G_I, IAlo, IAhi, MySLo, MySHi, Local_Idx, LDLocal)
C
c$$$      Write (LuOut, *) 'Accessing ', IALo, IAHi, MySLo, MySHi, LDLocal,
c$$$     $   NA, NI,
c$$$     $    IAHi-IAlo+1, NA*NI, LdLocal, NA * NI
      If ( IAHi-IAlo+1 .lt. NA*NI .OR. LdLocal .lt. NA * NI ) then
         Write (LuErr, *) 'XF3cI_Step12: Local dimensions of G_I ',
     $      ' inconsistent with arguments. Distribution may be wrong.'
         Call GA_Release( G_I, IAlo, IAhi, MySLo, MySHi)
         XF3cI_Step12 = .FALSE.
         Return
      EndIf
C
C     ******************************
C     * Cover the whole (pq| space *
C     ******************************
C     Having these loops on the outside (instead of SSh) minimizes the
C     communication required to get the patches of C we need.
C
      Status = .TRUE.
C
      Do Psh = 1, NAOSh
         Status = Status .AND. Bas_Cn2BFR( AOBas, PSh, Plo, Phi)
         NBFPSh = Phi - Plo +1
C
C        Get the matching block of the transformation matrix
C
         Call GA_Get(G_C, Plo, Phi, Alo, Ahi,
     $      Dbl_MB( CP_Idx), NBFPSh)
c$$$         Call DGeWr(Phi-Plo+1, AHi-ALo+1, Dbl_MB( CP_Idx),
c$$$     $      NBFPSh, Plo-1, Alo-1, LuOut, 'P tranformation', 80,
c$$$     $      'E15.6', Junk)
C
         Do Qsh = 1, NAOSh
            Status = Bas_Cn2BFR( AOBas, QSh, Qlo, Qhi)
            NBFQSh = Qhi - Qlo +1
C
C           Get the matching block of the transformation matrix
C
c$$$            Call GA_Print(G_C)
            Call GA_Get(G_C, Qlo, Qhi, Ilo, Ihi,
     $         Dbl_MB( CQ_Idx), NBFQSh)
c$$$            Call DGeWr(Qhi-Qlo+1, IHi-ILo+1, Dbl_MB( CQ_Idx),
c$$$     $         NBFQSh, Qlo-1, Ilo-1, LuOut, 'Q tranformation', 80,
c$$$     $         'E15.6', Junk)
C
C           ********************************************            
C           * Step through all |s) shells on this node *
C           ********************************************
C
            Do SSh = SShlo, SShhi
               Status = Status .AND. Bas_Cn2BFR( FitBas, SSh, Slo, Shi)
               NBFSSh = Shi - Slo +1
C
C              ********************************
C              * Evaluate this integral block *
C              ********************************
C
               If ( UseERI ) then
                  Call Int_2e3c( FitBas, SSh, AOBas, PSh, QSh,
     $               LScr, Dbl_MB( Scr_Idx), LInt, Dbl_MB( Int_Idx))
               Else
                  Call Int_1e3ov( FitBas, SSh, AOBas, PSh, AOBAs, QSh,
     $               LScr, Dbl_MB( Scr_Idx), LInt, Dbl_MB( Int_Idx))
               EndIf
C
C              The integrals come out (s|pq), but we want (pq|s).
C
c$$$               Write (LuOut, *) 'PSh, QSh, SSh = ', PSh, QSh, SSh
c$$$               Write (LuOut, *) 'NBFPSh, NBFQSh, NBFSSh = ',
c$$$     $            NBFPSh, NBFQSh, NBFSSh
c$$$               Call DGeWr(NBFPSh*NBFQSh, NBFSSh, Dbl_MB( Int_Idx),
c$$$     $            NBFPSh*NBFQSh, 0, 0, LuOut, 'AO integrals', 80,
c$$$     $            'E15.6', Junk)
c$$$               Call DGECOP('T', NBFPSh*NBFQSh, NBFSSh,
c$$$     $            Dbl_MB( Int_Idx), NBFSSh,
c$$$     $            Dbl_MB( Scr_Idx), NBFPSh*NBFQSh)
               Call DGECOP('N', NBFPSh*NBFQSh, NBFSSh,
     $            Dbl_MB( Int_Idx), NBFPSh*NBFQSh,
     $            Dbl_MB( Scr_Idx), NBFPSh*NBFQSh)
C
c$$$               Write (LuOut, *) 'PSh, QSh, SSh = ', PSh, QSh, SSh
c$$$               Call DGeWr(NBFPSh*NBFQSh, NBFSSh, Dbl_MB( Scr_Idx),
c$$$     $            NBFPSh*NBFQSh, 0, 0, LuOut, 'AO integrals', 80,
c$$$     $            'E15.6', Junk)
C
C              ******************************************
C              * Transform this block (pq|s) --> (ai|s) *
C              ******************************************
C              Use Int as scratch (we've made sure its large enough)
C              and put the results directly into the global array G_I
C
               Do S = 1, NBFSSh
C
C                 C_{pa}^T (pq|s) --> (aq|s)
C
                  Call dgemm('T', 'N', NA, NBFQSh, NBFPSh,
     $               FP_One, Dbl_MB( CP_Idx), NBFPSh, 
     $               Dbl_MB( Scr_Idx + (S-1)*NBFPSh*NBFQSh), NBFPSh,
     $               FP_Zero, Dbl_MB( Int_Idx), NA)
C
C                 (aq|s) --> (ai|s)
C
                  Call dgemm('N', 'N', NA, NI, NBFQSh,
     $               FP_One, Dbl_MB( Int_Idx), NA,
     $                       Dbl_MB( CQ_Idx ), NBFQSh, 
C
C                   Offet below derives from
C                   a) Offset for this node, -(MySLo-1)
C                   b) Offset for this shell, SLo-1
C                   c) Offset within the shell, S
C                   d) Offest of first column must be 0, -1
C
     $               FP_One, Dbl_MB( Local_Idx+(S+SLo-MySLo-1)*LDLocal),
     $               NA)
c$$$     $               FP_One, Dbl_MB( Local_Idx + (S+SLo-2)*NA*NI), NA)
C
c$$$               Write (LuOut, *) 'PSh, QSh, SSh, S = ', PSh, QSh, SSh, S
c$$$                  Call GA_Print(G_C)
               EndDo ! S
C
            EndDo               ! SSh
c$$$C
c$$$            Call DGeWr(NA*NI, MySHi-MySLo+1, Dbl_MB( Local_Idx),
c$$$     $         NA*NI, 0, 0, LuOut, 'MO integrals', 80,
c$$$     $         'E15.6', Junk)
         EndDo                  ! QSh
C
      EndDo                     ! PSh
C
C     ****************************
C     * Clean up after ourselves *
C     ****************************
C
      Call GA_Release_Update( G_I, IAlo, IAhi, MySLo, MySHi)
C
      If (.NOT. Status ) then
         Write (LuErr, *) 'XF3cI_Step12: Problem with basis info ',
     $      'queries in loops.'
         XF3cI_Step12 = .FALSE.
         Return
      EndIf
C
      Status = Status .AND. MA_Pop_Stack( H_CQ)
      Status = Status .AND. MA_Pop_Stack( H_CP)
      Status = Status .AND. MA_Pop_Stack( H_Scr)
      Status = Status .AND. MA_Pop_Stack( H_Int)
C
      If (.NOT. Status ) then
         Write (LuErr, *) 'XF3cI_Step12: Problem freeing local ',
     $      'memory.'
         XF3cI_Step12 = .FALSE.
         Return
      EndIf
C
      Call GA_Sync
C
      XF3cI_Step12 = .TRUE.
C
      Return
      End

