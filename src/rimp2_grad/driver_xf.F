C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_XF -- Three-center integral transformation driver for RI-MP2
C
C REVISION
C     $Id: driver_xf.F,v 1.1 1998-07-06 21:09:46 bernhold Exp $
C
C SYNOPSIS
      Logical Function RIMP2_Driver_XF(MaxSpin, RTDB, Geom, AOBas,
     $     FitBas,  TopSpin, D_Int, ALo, AHi, ILo, IHi,
     $   NILimit, Neglect, MinEVal, Approx, G_C, FN_Int)
      Implicit NONE
      Integer MaxSpin, RTDB, Geom, AOBas, FitBas, TopSpin, NILimit
      Integer ALo(TopSpin), AHi(TopSpin), ILo(TopSpin), IHi(TopSpin)
      Double precision Neglect, MinEVal
      Character*(*) Approx
      Integer G_C(TopSpin), D_Int(MaxSpin, 2, 2)
      Character*(*) FN_Int(MaxSpin, 2, 2)
C
C ARGUMENTS
C     AOBas     Handle for AO basis [IN]
C     FitBas    Handle for fitting basis [IN]
C     EVal      SCF eigenvalues [OUT]
C     NAOBF     Number of functions in AOBas [IN]
C     TopSpin   Number of spin cases (1 for RHF, 2 for UHF...) [IN]
C     D_Int     Handles for integral files [INOUT]
C     FNInt     File name for transformed integrals [IN]
C     ALo       Lower limit of virtual orbital range, by spin case [OUT]
C     AHi       Upper limit of virtual orbital range, by spin case [OUT]
C     ILo       Lower limit of occupied orbital range, by spin case [OUT]
C     IHi       Upper limit of occupied orbital range, by spin case [OUT]
C     NILiimit  Limit on the number of I blocks which be treated
C               simultaneously.  Determines maximum size of G_Int
C               and G_Out [IN]
C     Neglect   Threshold for neglect of integrals [IN]
C     MinEVal   Smallest eigenvalue to be retained in matrix inversions for
C               fiting procedure [IN]
C     Approx    RI approximation choice [IN]
C     G_C       Handles for transformation matrices [IN]
C
C MEMORY USE
C     Brief use of NAOBF local memory.
C     Creates several large global arrays for use in the transformation
C     Transformation routines themselves use additional memory.
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "numerical_constants.fh"
#include "stdio.fh"
#include "context.fh"
#include "bas.fh"
#include "geom.fh"
#include "rimp2g_timers.fh"
#include "msgids.fh"
#include "util.fh"
#include "inp.fh"
#include "dra.fh"
C
C EXTERNAL ROUTINES
      Logical GA_Create_One_By_Shells, XForm_3cInts
      External GA_Create_One_By_Shells, XForm_3cInts
C
      Logical Mk_Fit_XF
      External Mk_Fit_XF
C
c$$$      Logical Mk_V, Bas_NumBF
c$$$      External Mk_V, Bas_NumBF
c$$$      Double precision Max_V_Diagonal
c$$$      External Max_V_Diagonal
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
C     GA_DGEMM is currently using 128x128 blocking internally, so we
C     want a minimum block size that is going to make sure we are
C     not forcing too many small messages.
C
      Integer MinChunk
      Parameter ( MinChunk = 64 )
C
C LOCAL VARIABLES
C
      Logical Status
      Integer ISpin, NA, NI, FitLen
      Integer G_V, G_Out
      Double precision MyMaxV, Thresh
      Integer Initial_Avail, Final_Avail, Block_Size, NBlocks
      Integer LInt, LScr, H_Int, H_Scr, I_Int, I_Scr
      Character*1 Sp
      Character*256 TmpName, FilNam
C
C     ********************************************
C     * Compute the fitting basis transformation *
C     ********************************************
C
      Call GA_Sync
      Call PStat_On( PS_Fit_xf)
      Status = .TRUE.
      Status = Status .AND. Mk_Fit_XF( Approx, 'Y',
     $   FitBas, MinEVal, G_V)
C
      If (.NOT. Status )
     $   Call ErrQuit('Failed to compute V^{-1/2}', 0)
      Call PStat_Off( PS_Fit_xf )
C
      If ( Util_Print('fit xf', Print_Debug) ) Call GA_Print(G_V)
C
C     *************************************
C     * Initialize for integral screening *
C     *************************************
C     Schwarz' inequality says that |(pq|s)|^2 <= (pq|pq)(s|s), so we
c     can use the same mechanism as a 4c integral code would for the
c     first part...
C
      Call PStat_On( PS_Schwarz)
      Call Schwarz_Init( Geom, AOBas)
C
C     Along with an idea of the largest (s|s).  Need an appropriately
c     distributed GA to to act as a template so we can determine which s
c     shells each node will have to worry about.
C
      Call Max_V_Diagonal(FitBas, MyMaxV)
C
C     Adjust the integral neglect threshold by the maximum value
C     of the (s|s) integral to simplify Schwarz screening within the
C     children.  Note that if this node has no work to do, Max_v_diagonal
C     will return 0 so we should protect against it.
C
      If ( MyMaxV .gt. FP_Zero) then
         Thresh = Neglect / Sqrt( MyMaxV )
      Else
         Thresh = Neglect
      EndIf
      Call PStat_Off( PS_Schwarz)
C
C     *********************************************
C     * Perform transformation for each spin case *
C     *********************************************
C
      Do ISpin = 1, TopSpin
C
C        Compute sizes for convenience
C
         NA = AHi(ISpin) - ALo(ISpin) + 1
         NI = IHi(ISpin) - Ilo(ISpin) + 1
C
C        ********************************************
C        * Determine space how large out GAs can be *
C        ********************************************
C        Temporarily allocate space for the integral evaluation so we
C        are certain there will be enough later, when we really need it.
C
         If ( Approx .eq. 'V' ) then
            Call Int_Mem_2e3c( LInt, LScr)
         ElseIf (Approx .eq. 'SVS') then
            Call Int_Mem_3ov( LInt, LScr)
         Else
            Write (LuErr, *) 'Unknown approximation ', Approx
            Call ErrQuit('Rimp2_driver_xf: Unrecognized approximation',
     $         0)
         EndIf
C         Write (0, *) 'Node ', GA_NodeID(), LInt, LScr
C
         If ( .NOT. MA_Push_Get( MT_Dbl, LInt, 'Integral eval',
     $      H_Int, I_Int)) Call ErrQuit(
     $      'Rimp2_driver_xf: Unable to allocate local mem.', LInt)
         If ( .NOT. MA_Push_Get( MT_Dbl, LScr, 'Integral scr',
     $      H_Scr, I_Scr)) Call ErrQuit(
     $      'Rimp2_driver_xf: Unable to allocate local mem.', LScr)
C
C        Find out how much space a single occupied orbital
C        takes up.
C
         Initial_Avail = MA_SizeOf(MT_Byte, GA_Memory_Avail(), MT_Dbl)
C
         If ( .NOT. GA_Create_One_By_Shells( MT_Dbl, NA * 1,
     $      FitBas, 'Output integral matrix', G_Out) ) Call ErrQuit(
     $      'Rimp2_driver_xf: Unable to allocate minimal GA', ISpin)
C
         Final_Avail = MA_SizeOf(MT_Byte, GA_Memory_Avail(), MT_Dbl)
         Block_Size = Initial_Avail - Final_Avail
C
         If (.NOT. GA_Destroy( G_Out ) ) Call ErrQuit(
     $      'Rimp2_driver_xf: Unable to destroy a GA!', ISpin)
C
         If ( .NOT. MA_Pop_Stack( H_Scr) ) Call ErrQuit(
     $      'Rimp2_driver_xf: Unable to free local memory!', ISpin)
         If ( .NOT. MA_Pop_Stack( H_Int) ) Call ErrQuit(
     $      'Rimp2_driver_xf: Unable to free local memory!', ISpin)
C
C        Figure out how many blocks we can hold.  Remember we need
C        room for two matching arrays.
C
C        Also, block_size will be zero if there is no data on a node.
C        Must handle that possibility.
C
         If ( Block_Size .gt. 0) then
            NBlocks = Initial_Avail / (2 * Block_Size)
         Else
            NBlocks = NI
         EndIf
C
C        Find the smallest number of blocks across the machine
C
         Call GA_IGOp(Msg_RIMP2G_Mem, NBlocks, 1, 'min')
C
         If ( NILimit .le. 0) then
            NILimit = NBlocks
         Else
            NILimit = Min( NILimit, NBlocks)
         EndIf
C
C        **************************************************
C        * Create global arrays for use in transformation *
C        **************************************************
C
         Status = .TRUE.
         Status = Status .AND.
     $      GA_Create_One_By_Shells( MT_Dbl, NA * Min(NI, NILimit),
     $      FitBas, 'Output integral matrix', G_Out)
C
         If ( .NOT. Status )
     $      Call ErrQuit(
     $      'Rimp2_driver_xf:  Can''t allocate global arrays',
     $      ISpin)
C
C        *****************************************************
C        * Create the file to hold the transformed integrals *
C        *****************************************************
C        Setup filename. This depends on ASCII character set, but is that
C        a real limitation these days?
C
         If (.NOT. Bas_NumBF(FitBas, FitLen) ) Call ErrQuit(
     $      'Rimp2_driver_xf: can''t get basis info', FitBas)
         If ( DRA_Create(MT_Dbl, NI*NA, FitLen, 'MO 3-center integrals',
     $      FN_Int(ISpin, 1, 1)(:Inp_StrLen(FN_Int(ISpin, 1, 1))),
     $        DRA_W, NA, FitLen, D_Int(ISpin, 1, 1) ) .ne. 0)
     $      Call Errquit('Rimp2_driver_xf: can''t open integral file',
     $      ISpin)
C
C        *************************
C        * Do the transformation *
C        *************************
C
         If (Approx .eq. 'SVS') then
            Status = XForm_3cInts(RTDB, Geom, 'OVERLAP', AOBas, FitBas,
     $         ALo(ISpin), AHi(ISpin), ILo(ISpin), IHi(ISpin),
     $         G_C(ISpin), 'YES', G_V, 'DISK', G_Out,
     $           D_Int(ISpin, 1, 1), Thresh )
         ElseIf (Approx .eq. 'V') then
            Status = XForm_3cInts(RTDB, Geom, 'ERI', AOBas, FitBas,
     $         ALo(ISpin), AHi(ISpin), ILo(ISpin), IHi(ISpin),
     $         G_C(ISpin), 'Yes', G_V, 'DISK', G_Out,
     $           D_Int(ISpin, 1, 1), Thresh )
         Else
            Write (LuErr, *) 'Unknown approximation ', Approx
            Call ErrQuit('Rimp2_driver_xf: Unrecognized approximation',
     $         0)
         EndIf
C
C        Do ket integrals if necessary
C
         If ( FN_Int(ISpin, 1, 1) .ne. FN_Int(ISpin, 2, 1) ) then
            If ( DRA_Create(MT_Dbl, NI*NA, FitLen,
     $           'MO 3-center integrals',
     $           FN_Int(ISpin, 2, 1)(:Inp_StrLen(FN_Int(ISpin, 2, 1))),
     $           DRA_W, NA, FitLen, D_Int(ISpin, 2, 1) ) .ne. 0)
     $           Call Errquit(
     $           'Rimp2_driver_xf: can''t open integral file', ISpin)
C
            If (Approx .eq. 'SVS') then
               Status = XForm_3cInts(RTDB, Geom, 'OVERLAP', AOBas,
     $              FitBas, ALo(ISpin), AHi(ISpin), ILo(ISpin),
     $              IHi(ISpin), G_C(ISpin), 'YES', G_V, 'DISK', G_Out,
     $              D_Int(ISpin, 2, 1), Thresh )
            ElseIf (Approx .eq. 'V') then
               Status = XForm_3cInts(RTDB, Geom, 'ERI', AOBas,
     $              FitBas, ALo(ISpin), AHi(ISpin), ILo(ISpin),
     $              IHi(ISpin), G_C(ISpin), 'YES', G_V, 'DISK', G_Out,
     $              D_Int(ISpin, 2, 1), Thresh )
            Else
               Write (LuErr, *) 'Unknown approximation ', Approx
               Call ErrQuit(
     $              'Rimp2_driver_xf: Unrecognized approximation', 0)
            EndIf
         EndIf ! Bra and ket different
C
         If ( .NOT. Status ) Call ErrQuit(
     $      'Rimp2_driver_xf: Transformation trouble for spin', ISpin)
C
C        Free global arrays we are done with
C
         Status = .TRUE.
         Status = Status .AND. GA_Destroy( G_Out )
C
         If ( .NOT. Status ) Call ErrQuit(
     $      'Rimp2_driver_xf:  Unable to free global arrays',
     $      ISpin)
C
         If ( DRA_Close(D_Int(ISpin, 1, 1) ) .ne. 0) Call ErrQuit(
     $      'Rimp2_driver_xf: can''t close integral file', ISpin)
C
         If ( FN_Int(ISpin, 1, 1) .ne. FN_Int(ISpin, 2, 1) ) then
            If ( DRA_Close(D_Int(ISpin, 2, 1) ) .ne. 0) Call ErrQuit(
     $           'Rimp2_driver_xf: can''t close integral file', ISpin)
        EndIf
      EndDo ! ISpin
C
C     Cleanup after the whole transformation is complete
C
      Status = .TRUE.
      Status = Status .AND. GA_Destroy( G_V )
C
      If ( .NOT. Status ) Call ErrQuit(
     $   'Rimp2_driver_xf:  Unable to free GAs ', 0)
C
      Call Schwarz_Tidy
C
C     ***********
C     * Go home *
C     ***********
C
      Rimp2_driver_xf = .TRUE.
C
      Return
      End


