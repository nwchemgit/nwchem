C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     XF3cI_Step12 -- Do the first two steps of the three-center
C     integral transformation
C
C REVISION
C     $Id: xf3ci_step12.F,v 1.7 1994-10-14 19:44:56 gg502 Exp $
C
C SYNOPSIS
      Logical Function XF3cI_Step12(IntTyp, Alo, Ahi, Ilo, Ihi,
     $   AOBas, G_C, FitBas, G_I, Geom, Neglect)
      Implicit NONE
      Character*1 IntTyp
      Integer Alo, Ahi, Ilo, Ihi, FitBas, AOBas, G_I, G_C, Geom
      Double precision Neglect
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
C     Neglect  Neglect raw integrals smaller than this value [IN]
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
#include "geom.fh"
#include "schwarz.fh"
#include "rimp2_timers.fh"
C
C EXTERNAL ROUTINES
      Logical LSame      ! From BLAS-[23]
      External LSame
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer IAlo, IAhi, NI, NA
      Integer MxBFSh, LInt, LIntScr, LScr, MySLo, MySHi, MxFitSh
      Integer SShlo, SShhi, SSh, First, Last
      Logical Status, UseERI
      Integer Local_Idx, LDLocal, H_Int, H_Scr, Int_Idx, Scr_Idx
      Integer H_IntScr, IntScr_Idx, H_Scr2, Scr2_Idx, LScr2
      Integer H_CP, H_CQ, CP_Idx, CQ_Idx, LCP, LCQ
      Integer NCent, CenP, CenQ, MxBFCe, PBFLo, PBFHi, QBFLo,
     $   QBFHi
      Integer H_CP2, H_CQ2, CP2_Idx, CQ2_Idx
      Integer PSh, PShLo, PShHi, PLo, PHi, NBFPSh
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
C
C     If there is no data on this node, there is no work fo us to do!
C
      If ( MySLo .gt. 0 .AND. MySHi .gt. 0) then
         Status = .TRUE.
         Status = Status .AND. Bas_BF2Cn( FitBas, MySLo, SShlo)
         Status = Status .AND. Bas_BF2Cn( FitBas, MySHi, SShhi)
C
C        Find the largest shell we own
C
         MxFitSh = 0
         Do SSh = SShLo, SShHi
            Status = Status .AND. Bas_Cn2Bfr( FitBas, SSh, First, Last)
            MxFitSh = Max( MxFitSh, Last-First+1)
         EndDo
C
C        We also need some info about the AO basis
C
         Status = Status .AND. Bas_NBF_Cn_Max( AOBas, MxBFSh)
         Status = Status .AND. Bas_NBF_Ce_Max( AOBas, MxBFCe)
C
         If (.NOT. Status ) then
            Write (LuErr, *) 'XF3cI_Step12: Problems with basis set ',
     $         'queries.'
            XF3cI_Step12 = .FALSE.
            Return
         EndIf
C
C        **********************************
C        * Allocate local "scratch" space *
C        **********************************
C        These guys are both fed into the integral package, and reused
C        during the transformation, so we must be sure to size for both
C        uses.
C
C        This is what the integral package needs
C
         If ( UseERI ) then
            Call Int_Mem_2e3c( LInt, LIntScr)
         Else
            Call Int_Mem_3ov( LInt, LIntScr)
         EndIf
C
         Status = .TRUE.
         Status = Status .AND. MA_Push_Get( MT_Dbl, LInt,
     $      'Integral eval', H_Int, Int_Idx)
         Status = Status .AND. MA_Push_Get( MT_Dbl, LIntScr,
     $      'Integral scr', H_IntScr, IntScr_Idx)
C
C        We also need to localize bits of the transformation matrix, C
C
         LCP = NA * MxBFCe
         LCQ = NI * MxBFCe
C
         Status = Status .AND. MA_Push_Get( MT_Dbl, LCP, 'trans. scr'
     $      ,H_CP, CP_Idx)
         Status = Status .AND. MA_Push_Get( MT_Dbl, LCQ, 'trans. scr',
     $      H_CQ, CQ_Idx)
C
         Status = Status .AND. MA_Push_Get( MT_Dbl, LCQ, 'trans. scr',
     $      H_CP2, CP2_Idx)
         Status = Status .AND. MA_Push_Get( MT_Dbl, LCP, 'trans. scr',
     $      H_CQ2, CQ2_Idx)
C
C        And some scratch space
C
         LScr  = MxBFSh * NI * (MySHi-MySLo+1)
         LScr2 = NI * MxBFSh * MxFitSh
C
         Status = Status .AND. MA_Push_Get( MT_Dbl, LScr,
     $      'Transformation scr', H_Scr, Scr_Idx)
         Status = Status .AND. MA_Push_Get( MT_Dbl, LScr2,
     $      'Transformation scr', H_Scr2, Scr2_Idx)
C
         If (.NOT. Status ) then
            Write (LuErr, *) 'XF3cI_Step12: Problem allocating local ',
     $         'memory.'
            XF3cI_Step12 = .FALSE.
            Return
         EndIf
C
C        ********************************************
C        * Grab our panel of G_I to access directly *
C        ********************************************
C        This avoids the need for a scratch matrix of size NA*NI
C        and a copy operation.  But the leading dimension
C        of the local patch must be a multiple of NA, or the second
C        transformation step will be wrong.
C
         Call GA_Access( G_I, IAlo, IAhi, MySLo, MySHi, Local_Idx,
     $      LDLocal)
C
         If ( IAHi-IAlo+1 .lt. NA*NI .OR. LdLocal .ne. NA * NI ) then
            Write (LuErr, *) 'XF3cI_Step12: Local dimensions of G_I ',
     $         ' inconsistent with arguments. Distribution may be ',
     $         'wrong.'
            Call GA_Release( G_I, IAlo, IAhi, MySLo, MySHi)
            XF3cI_Step12 = .FALSE.
            Return
         EndIf
C
         If ( .NOT. Geom_NCent( Geom, NCent) ) then
            Write (LuErr, *) 'XF3cI_Step12: Problem getting number of ',
     $         'centers.'
            XF3cI_Step12 = .FALSE.
            Return
         EndIf
c$$$C
c$$$         Write (0, *) 'Starting actual transformation'
c$$$         Call Util_Flush(0)
C
         Do CenP = 1, NCent
            Status = Bas_Ce2BFR( AOBas, CenP, PBFLo, PBFHi)
            Call GA_Get(G_C, PBFlo, PBFhi, Alo, Ahi,
     $         Dbl_MB(CP_Idx), MxBFCe)
            Call GA_Get(G_C, PBFlo, PBFhi, Ilo, Ihi,
     $         Dbl_MB(CP2_Idx), MxBFCe)
C
            Status = Bas_Ce2CnR( AOBas, CenP, PShLo, PShHi)
            Do PSh = PShLo, PShHi
               Status = Bas_Cn2BFR( AOBas, PSh, Plo, Phi)
               NBFPSh = Phi - Plo +1
C
               Call DFill(LScr, FP_Zero, Dbl_MB(Scr_Idx), 1)
C
               Do CenQ = 1, NCent
                  Status = Bas_Ce2BFR( AOBas, CenQ, QBFLo, QBFHi)
C
                  If ( Schwarz_Atom(CenP, CenQ) * Schwarz_Max() .gt.
     $               Neglect ) then
                     Call GA_Get(G_C, QBFlo, QBFhi, Ilo, Ihi,
     $                  Dbl_MB(CQ_Idx), MxBFCe)
                     Call GA_Get(G_C, QBFlo, QBFhi, Alo, Ahi,
     $                  Dbl_MB(CQ2_Idx), MxBFCe)
C
                     Call XF3cI_Step12a(UseERI, AOBas, CenP, PSh, CenQ,
     $                  Dbl_MB(CP_Idx), MxBFCe, NA,
     $                  Dbl_MB(CP2_Idx), MxBFCe, NI,
     $                  Dbl_MB(CQ_Idx), MxBFCe, NI,
     $                  Dbl_MB(CQ2_Idx), MxBFCe, NA,
     $                  FitBas, SShLo, SShhi,
     $                  Dbl_MB(Local_Idx), LDLocal, 
     $                  Dbl_MB(Int_Idx), LInt,
     $                  Dbl_MB(IntScr_Idx), LIntScr,
     $                  Dbl_MB(Scr_Idx), MxBFSh, MySLo, MySHi,
     $                  Dbl_MB(Scr2_Idx), LScr2,
     $                  Neglect)
                  EndIf ! Schwarz screening
C
c$$$                  Write(0,*) 'Done with Q center ', CenP, CenQ
c$$$                  Call Util_Flush(0)
               EndDo    ! CenQ
C
C              C_{pa}^T (pi|s) --> (ai|s)
C
C              Notice that we are thinking of the result we generate
C              as an NA by NI*NS matrix, while in reality, it is
C              an NA*NI by NS matrix.  Therefore, of the leading dimension
C              is not exactly NA*NI, we're in trouble!  This is checked
C              when we get access to the local patch, above.
C
               Call PStat_On( PS_Step12)
C
               Call sgemm('T', 'N', NA, NI*(MySHi-MySLo+1), NBFPSh,
     $            FP_One, Dbl_MB(CP_Idx + PLo-PBFLo), MxBfCe,
     $            Dbl_MB(Scr_Idx), NBFPSh,
     $            FP_One, Dbl_MB(Local_Idx), NA)
               Call MA_Verify_Allocator_Stuff
C
               Call PStat_Off( PS_Step12)
C
            EndDo               ! PSh
c$$$            Write(0,*) 'Done with P center ', CenP
c$$$            Call Util_Flush(0)
         EndDo                  ! CenP
C
C        ****************************
C        * Clean up after ourselves *
C        ****************************
C
         Call GA_Release_Update( G_I, IAlo, IAhi, MySLo, MySHi)
C
         Status = Status .AND. MA_Pop_Stack( H_Scr2)
         Status = Status .AND. MA_Pop_Stack( H_Scr)
         Status = Status .AND. MA_Pop_Stack( H_CQ2)
         Status = Status .AND. MA_Pop_Stack( H_CP2)
         Status = Status .AND. MA_Pop_Stack( H_CQ)
         Status = Status .AND. MA_Pop_Stack( H_CP)
         Status = Status .AND. MA_Pop_Stack( H_IntScr)
         Status = Status .AND. MA_Pop_Stack( H_Int)
C
         If (.NOT. Status ) then
            Write (LuErr, *) 'XF3cI_Step12: Problem freeing local ',
     $         'memory.'
            XF3cI_Step12 = .FALSE.
            Return
         EndIf
C
      EndIf
C
      Call GA_Sync
C
      XF3cI_Step12 = .TRUE.
c$$$C
c$$$      Write( 0,*) 'Leaving first 2/3'
c$$$      Call Util_Flush(0)
C
      Return
      End

