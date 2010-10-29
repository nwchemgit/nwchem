C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     XF3cI_Step12 -- Do the first two steps of the three-center
C     integral transformation
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function XF3cI_Step12(RTDB, Geom, IntTyp,
     $   Alo, Ahi, Ilo, Ihi, AOBas, G_C, FitBas, G_I, Neglect)
      Implicit NONE
#include "errquit.fh"
      Character*1 IntTyp
      Integer RTDB, Alo, Ahi, Ilo, Ihi, FitBas, AOBas, G_I, G_C, Geom
      Double precision Neglect
C
C ARGUMENTS
C     RTDB     Handle to the run-time database [IN]
C     Geom     Handle to the geometry [IN]
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
C     Primary memory requirements are for G_I, the global array that 
C     will hold the 2/3 transformed integrals.  Local memory is used
C     to increase matrix multiplication performance for the actual
C     transformation.  Routine XF3cI_Mem12 will do its best to use
C     all available local memory subject to external limitations
C     which may be specified in the RTDB.
C     
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "numerical_constants.fh"
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "rimp2_timers.fh"
#include "rtdb.fh"
C
C EXTERNAL ROUTINES
      Logical LSame      ! From BLAS-[23]
      External LSame
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer IAlo, IAhi, NI, NA, NF
      Integer MxBFSh, LInt, LScr, LOut, NodeSLo, NodeSHi, MxFitSh
      Integer SShlo, SShhi, SSh, First, Last
      Logical Status, UseERI
      Integer Local_Idx, LDLocal, H_Int, H_Out, I_Int, I_Out
      Integer H_Scr, I_Scr, H_Interm, I_Interm, LInterm
      Integer H_CP, H_CQ, I_CP, I_CQ, LCP, LCQ
      Integer MxBFCe
      Integer NAOBF
      Integer LIBuf, H_IBuf, I_IBuf
      Integer PMin, QMin, SMin, PMax, QMax, SMax, PDim, QDim, SDim
      Integer MemAvail, MemUsed
#ifdef  BAD_GACCESS 
      integer l_local
#endif
      Double Precision FLOPs
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
C     ********************************
C     * Some AO basis size inquiries *
C     ********************************
C     Needed before RTDB calls
C
      Status = .TRUE.
      Status = Status .AND. Bas_NBF_Cn_Max( AOBas, MxBFSh)
      Status = Status .AND. Bas_NBF_Ce_Max( AOBas, MxBFCe)
      Status = Status .AND. Bas_NumBF( AOBas, NAOBF)
C
      If (.NOT. Status ) Call ErrQuit(
     $   'XF3cI_Step12: Problems with basis set queries', 0, BASIS_ERR)
C
C     ********************************************
C     * Determine our range of the fitting basis *
C     ********************************************
C     Each node will fill the range that lives in its memory.  So we
C     get the basis functions that define the local panel and translate
C     them to shell blocks.
C
      Call GA_Distribution( G_I, GA_NodeID(), IAlo, IAhi,
     $   NodeSLo, NodeSHi)
      NF = NodeSHi - NodeSLo +1
C
C     ***************************************************************
C     * Make some database inquiries before we go our separate ways *
C     ***************************************************************
C     This allows the user to manpulate the way we do things, but it
C     will probaly happen only in pretty wierd circumstances.
C
C     Note that Q is the first index transformed, and most important to
c     performance, therefore is it labeled as "AO 1".
C
      If ( .NOT. RTDB_Get( RTDB, 'xf3ci:AO 1 batch size',
     $   MT_Int, 1, QMax) ) QMax = NAOBF
      If ( .NOT. RTDB_Get( RTDB, 'xf3ci:AO 2 batch size',
     $   MT_Int, 1, PMax) ) PMax = NAOBF
      If ( .NOT. RTDB_Get( RTDB, 'xf3ci:fit batch size',
     $   MT_Int, 1, SMax) ) SMax = NF
C
      PMax = Min( PMax, NAOBF)
      QMax = Min( QMax, NAOBF)
      SMax = Min( SMax, NF)
C
C     *******************
C     * Exit if no work *
C     *******************
C     Must have everyone to read database, but nodes with no fitting
C     basis functions don't want to do the basis set inquiries nor
C     try to allocate local memory (MA doesn't like zero lengths).
C
      If ( NF .le. 0 ) then
         XF3cI_Step12 = .TRUE.
         Call GA_Sync
         Return
      EndIf
C
C     ***************************************
C     * Get other fitting basis information *
C     ***************************************
C     Find the largest shell in our range of the fitting basis
C      
      Status = .TRUE.
      Status = Status .AND. Bas_BF2Cn( FitBas, NodeSLo, SShlo)
      Status = Status .AND. Bas_BF2Cn( FitBas, NodeSHi, SShhi)
C
      MxFitSh = 0
      Do SSh = SShLo, SShHi
         Status = Status .AND. Bas_Cn2Bfr( FitBas, SSh, First, Last)
         MxFitSh = Max( MxFitSh, Last-First+1)
      EndDo
C
      If (.NOT. Status ) Call ErrQuit(
     $   'XF3cI_Step12: Problems with basis set queries', 1, BASIS_ERR)
C
C     ********************************
C     * Setup for 3-center integrals *
C     ********************************
C
      If ( UseERI ) then
         Call Int_Mem_2e3c( LInt, LScr)
      Else
         Call Int_Mem_3ov( LInt, LScr)
      EndIf
C
      Status = .TRUE.
      Status = Status .AND. MA_Push_Get( MT_Dbl, LInt,
     $   'Integral eval', H_Int, I_Int)
      Status = Status .AND. MA_Push_Get( MT_Dbl, LScr,
     $   'Integral scr', H_Scr, I_Scr)
C
C     ***************************************
C     * Determine scratch space size to use *
C     ***************************************
C     This transformation can use varying amounts of local memory.  More
c     memory allows largers blocks of integrals to be evaluated and
c     transformed at a time, which increases the efficiency of matrix
c     multiplication.  This can be a major performance win, so we try
C     to use as much as we can -- intelligently.
C
C     The minimum sizes are the largest shell block in each basis
C
      PMin = MxBFSh
      QMin = MxBFSh
      SMin = MxFitSh
C
C     Insure that the maximums are reasonable. Maybe we should report to
c     the user when we can't obey their minimums?
C
      PMax = Max( PMax, PMin)
      QMax = Max( QMax, QMin)
      SMax = Max( SMax, SMin)
C
C     Now inquire how much space is available to us (everything we need
c     to allocate is of type double)
C
      MemAvail = MA_Inquire_Avail( MT_Dbl)
C
C     Figure out how best to use this memory.  Encapsulated in a
c     separate routine for readibility and returns the individual array
c     sizes required in addition to the base dimensions for safety.
c     That way you are less likely to try to sneak in an extra array
c     without telling the routine that does the negotiation about it.
C     Not that I don't trust you, of course...
C
      Call XF3cI_Mem12(MemAvail, PMin, QMin, SMin, PMax, QMax, SMax, NA,
     $   NI, PDim, QDim, SDim, MemUsed, LCP, LCQ, LInterm, LOut,
     $   LIBuf)
c$$$C
c$$$      Write (LuOut, *) 'XF3cI_Step12: Node ', GA_NodeID(), ' using ',
c$$$     $   MemUsed, ' of ', MemAvail, ' with ', PDim, QDim, SDim
C
C     If there isn't enough memory the actual dimensions will be zero
c     and LUsed will contain the minimum memory required to proceed.
C
      If ( PDim * QDim * SDim .eq. 0) Call ErrQuit(
     $   'XF3cI_Step12: Need more local memory to continue (bytes)',
     $      MA_SizeOf(MT_Dbl, MemUsed-MemAvail, MT_Byte) , MA_ERR)
C
C     *************************
C     * Allocate local memory *
C     *************************
C
      Status = .TRUE.
C
      Status = Status .AND. MA_Push_Get( MT_Dbl, LCP,
     $   'local xform matrix 2', H_CP, I_CP)
      Status = Status .AND. MA_Push_Get( MT_Dbl, LCQ,
     $   'local xform matrix 1',  H_CQ, I_CQ)
C
      Status = Status .AND. MA_Push_Get( MT_Dbl, LOut,
     $   'xf3ci step1 ouput', H_Out, I_Out)
      Status = Status .AND. MA_Push_Get( MT_Dbl, LInterm,
     $   'xf3ci step1 intermediate', H_Interm, I_Interm)
      Status = Status .AND. MA_Push_Get( MT_Dbl, LIBuf,
     $   'xf3ci int buffer', H_IBuf, I_IBuf)
C
      If (.NOT. Status ) Call ErrQuit(
     $   'XF3cI_Step12: Problem allocating local memory', 0, MA_ERR)
C
C     ********************************************
C     * Grab our panel of G_I to access directly *
C     ********************************************
C     This avoids the need for a scratch matrix of size NA*NI
C     and a copy operation.  But the leading dimension
C     of the local patch must be a multiple of NA, or the second
C     transformation step will be wrong.
C
#ifdef  BAD_GACCESS 
      if(.not.ma_push_get(MT_DBL,(iahi-ialo+1)*(nodeshi-nodeslo+1),
     $  'scratch buff', l_local, local_idx)) call
     $  errquit('step12: pushget failed',0,0)
      ldlocal=na*ni
      call ga_get(g_i,ialo,iahi,nodeslo,nodeshi,
     $     dbl_mb(local_idx),ldlocal)
#else
      Call GA_Access( G_I, IAlo, IAhi, NodeSLo, NodeSHi, Local_Idx,
     $   LDLocal)
C
      If ( IAHi-IAlo+1 .lt. NA*NI .OR. LdLocal .ne. NA * NI ) then
         Write (LuErr, *) 'XF3cI_Step12: Local dimensions of G_I ',
     $      ' inconsistent with arguments. Distribution may be ',
     $      'wrong.'
         Call GA_Release( G_I, IAlo, IAhi, NodeSLo, NodeSHi)
         XF3cI_Step12 = .FALSE.
         Return
      EndIf
#endif
C
C     ***********************************************************
C     * After all this preparation, we finally get to the work! *
C     ***********************************************************
C
      FLOPs = FP_Zero
      Call XF3cI_Step12a(UseERI, AOBas, PDim,
     $   QDim, FitBas, SDim, Neglect,
     $   Dbl_MB(I_Int), LInt,
     $   Dbl_MB(I_Scr), LScr,
     $   Dbl_MB(I_IBuf),
     $   Dbl_MB(I_CP), ALo, AHi,
     $   Dbl_MB(I_CQ), ILo, IHi, G_C,
     $   Dbl_MB(I_Interm), Dbl_MB(I_Out), PDim,
     $   Dbl_MB(Local_Idx), LDLocal, NodeSLo, NodeSHi, FLOPs)
c$$$C
c$$$C     Eventually, this should be accumulated into PStat-managed memory
c$$$C
c$$$      Write (LuOut, *) 'XF3cI_Step12: Node ', GA_NodeID(),
c$$$     $   ' did approximately ', FLOPs, ' FLOPs'
C
C     ****************************
C     * Clean up after ourselves *
C     ****************************
C
#ifdef  BAD_GACCESS 
      call ga_put(g_i,ialo,iahi,nodeslo,nodeshi,
     $     dbl_mb(local_idx),ldlocal)
      if(.not.ma_pop_stack(l_local)) call 
     $  errquit('step12: popstack failed',0,0)
#else
      Call GA_Release_Update( G_I, IAlo, IAhi, NodeSLo, NodeSHi)
#endif
C
      Status = .TRUE.
      Status = Status .AND. MA_Pop_Stack( H_IBuf)
      Status = Status .AND. MA_Pop_Stack( H_Interm)
      Status = Status .AND. MA_Pop_Stack( H_Out)
      Status = Status .AND. MA_Pop_Stack( H_CQ)
      Status = Status .AND. MA_Pop_Stack( H_CP)
      Status = Status .AND. MA_Pop_Stack( H_Scr)
      Status = Status .AND. MA_Pop_Stack( H_Int)
C
      If (.NOT. Status ) Call ErrQuit(
     $   'XF3cI_Step12: Problem freeing local memory', 0, MA_ERR)
C
      Call GA_Sync
C
      XF3cI_Step12 = .TRUE.
C
      Return
      End

