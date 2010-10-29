C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_XF -- Three-center integral transformation driver for RI-MP2
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function RIMP2_Driver_XF(RTDB, Geom, AOBas,
     $     FitBas,  TopSpin, MaxSpin, D_Int,
     $     NFrzO, NAct, NVir, NFrzV, KeepFrz,
     $     ALo, AHi, ILo, IHi, 
     $   NILimit, Neglect, MinEVal, Approx, G_C,
     $     BraKetSame, DoGrad, IntRange)
      Implicit NONE
#include "errquit.fh"
      Integer RTDB, Geom, AOBas, FitBas, TopSpin, NILimit
      Integer ALo(TopSpin), AHi(TopSpin), ILo(TopSpin), IHi(TopSpin)
      Double precision Neglect, MinEVal
      Character*(*) Approx
      Integer G_C(TopSpin), MaxSpin, D_Int(MaxSpin, 2, 2)
      Logical BraKetSame
      Logical DoGrad
      Integer NFrzO ![in]
      Integer NAct(TopSpin) ![in]
      Integer NVir(TopSpin) ![in]
      Integer NFrzV(TopSpin) ![in]
      Logical KeepFrz ![in]
      Integer IntRange(2, 2, MaxSpin, 2, 2) ![in]
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
C DESCRIPTION
C     Integrals generated based on inputs computed in RIMP2_Prep_DRA
C
c$$$C     BraKetSame   DoGrad   D_Int()   Integrals
c$$$C
c$$$C     T            F        *, 1, 1   (ai|D')
c$$$C                           *, 2, 1      -
c$$$C                           *, 1, 2      -
c$$$C                           *, 2, 2      -
c$$$C
c$$$C     F            F        *, 1, 1   (ai|D)
c$$$C                           *, 2, 1   (ai|D")
c$$$C                           *, 1, 2      -
c$$$C                           *, 2, 2      -
c$$$C
c$$$C     T            T        *, 1, 1   (ap|D')
c$$$C                           *, 2, 1      -
c$$$C                           *, 1, 2   (ij|D')
c$$$C                           *, 2, 2      -
c$$$C
c$$$C     F            T        *, 1, 1   (ap|D)
c$$$C                           *, 2, 1   (ap|D")
c$$$C                           *, 1, 2   (ij|D)
c$$$C                           *, 2, 2   (ij|D")
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
      Logical GA_Create_One_By_Shells, XForm_3cIntsG
      External GA_Create_One_By_Shells, XForm_3cIntsG
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
      Integer MyMaxSpin
      Parameter (MyMaxSpin = 2)
C
C LOCAL VARIABLES
C
      Logical Status
      Integer ISpin, FitLen
      Integer G_V, G_Out
      Double precision MyMaxV, Thresh
      Character*1 Sp, SpinLabel
      Character*256 String1
      Character*3 XFFit
      Integer Class, BraKet, N1, N2, NPass
      Logical PrInfo, PrPrgRpt
c$$$C
c$$$      Integer Lim1(MyMaxSpin, 2, 2, 2)
c$$$      Integer Lim2(MyMaxSpin, 2, 2, 2)
C
C STATEMENT FUNCTIONS
      Character*1 SpinItoA
      SpinItoA(ISpin) = Char( ISpin + IChar('a') - 1 )
C
C     Sanity check
C
      If ( MyMaxSpin .ne. MaxSpin ) Call ErrQuit(
     $     'RIMP2_Driver_XF: MaxSpin does not match parent routine',
     $     MyMaxSpin, INPUT_ERR)
C
      NPass = 0
      PrInfo =   Util_Print('information', Print_Low)
      PrPrgRpt = Util_Print('progress reports', Print_Default)
C
      If ( PrInfo .or. PrPrgRpt .AND. GA_NodeID() .eq. 0) then
         If ( DoGrad ) then
            Write (LuOut, 9005) 
         Else
            Write (LuOut, 9000) 
         EndIf
      EndIf
C
 9000 Format(1X, 'Integral transformation for energy only')
 9005 Format(1X, 'Integral transformation for energy + gradient')
C
C     ********************************************
C     * Compute the fitting basis transformation *
C     ********************************************
C
C
      If (PrPrgRpt) Call Progress_Report(LuOut,
     $        '   Computing fitting transformation matrix')
      Call GA_Sync
      Call PStat_On( PS_Fit_xf)
      Sp = 'N'
      If ( BraKetSame ) Sp = 'Y'
      Status = .TRUE.
      Status = Status .AND. Mk_Fit_XF( Approx, Sp, FitBas, MinEVal, G_V)
C
      If (.NOT. Status )
     $   Call ErrQuit('Failed to compute fitting transformation', 0,
     &       CALC_ERR)
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
      If (.NOT. Bas_NumBF(FitBas, FitLen) ) Call ErrQuit(
     $     'RIMP2_Driver_XF: can''t get basis info', FitBas, BASIS_ERR)
C
      Do Class = 1, 2
         Do BraKet = 1, 2
            Do ISpin = 1, TopSpin
               SpinLabel = SpinItoA(ISpin)
c$$$C
c$$$C        Compute sizes for convenience
c$$$C
c$$$         NA = AHi(ISpin) - ALo(ISpin) + 1
c$$$         NI = IHi(ISpin) - Ilo(ISpin) + 1
c$$$C
c$$$C        ********************************************
c$$$C        * Determine space how large out GAs can be *
c$$$C        ********************************************
c$$$C        Temporarily allocate space for the integral evaluation so we
c$$$C        are certain there will be enough later, when we really need it.
c$$$C
c$$$         If ( Approx .eq. 'V' ) then
c$$$            Call Int_Mem_2e3c( LInt, LScr)
c$$$         ElseIf (Approx .eq. 'SVS') then
c$$$            Call Int_Mem_3ov( LInt, LScr)
c$$$         Else
c$$$            Write (LuErr, *) 'Unknown approximation ', Approx
c$$$            Call ErrQuit('Rimp2_driver_xf: Unrecognized approximation',
c$$$     $         0)
c$$$         EndIf
c$$$C         Write (0, *) 'Node ', GA_NodeID(), LInt, LScr
c$$$C
c$$$         If ( .NOT. MA_Push_Get( MT_Dbl, LInt, 'Integral eval',
c$$$     $      H_Int, I_Int)) Call ErrQuit(
c$$$     $      'Rimp2_driver_xf: Unable to allocate local mem.', LInt)
c$$$         If ( .NOT. MA_Push_Get( MT_Dbl, LScr, 'Integral scr',
c$$$     $      H_Scr, I_Scr)) Call ErrQuit(
c$$$     $      'Rimp2_driver_xf: Unable to allocate local mem.', LScr)
c$$$C
c$$$C        Find out how much space a single occupied orbital
c$$$C        takes up.
c$$$C
c$$$         Initial_Avail = MA_SizeOf(MT_Byte, GA_Memory_Avail(), MT_Dbl)
c$$$C
c$$$         If ( .NOT. GA_Create_One_By_Shells( MT_Dbl, NA * 1,
c$$$     $      FitBas, 'Output integral matrix', G_Out) ) Call ErrQuit(
c$$$     $      'Rimp2_driver_xf: Unable to allocate minimal GA', ISpin)
c$$$C
c$$$         Final_Avail = MA_SizeOf(MT_Byte, GA_Memory_Avail(), MT_Dbl)
c$$$         Block_Size = Initial_Avail - Final_Avail
c$$$C
c$$$         If (.NOT. GA_Destroy( G_Out ) ) Call ErrQuit(
c$$$     $      'Rimp2_driver_xf: Unable to destroy a GA!', ISpin)
c$$$C
c$$$         If ( .NOT. MA_Pop_Stack( H_Scr) ) Call ErrQuit(
c$$$     $      'Rimp2_driver_xf: Unable to free local memory!', ISpin)
c$$$         If ( .NOT. MA_Pop_Stack( H_Int) ) Call ErrQuit(
c$$$     $      'Rimp2_driver_xf: Unable to free local memory!', ISpin)
c$$$C
c$$$C        Figure out how many blocks we can hold.  Remember we need
c$$$C        room for two matching arrays.
c$$$C
c$$$C        Also, block_size will be zero if there is no data on a node.
c$$$C        Must handle that possibility.
c$$$C
c$$$         If ( Block_Size .gt. 0) then
c$$$            NBlocks = Initial_Avail / (2 * Block_Size)
c$$$         Else
c$$$            NBlocks = NI
c$$$         EndIf
c$$$C
c$$$C        Find the smallest number of blocks across the machine
c$$$C
c$$$         Call GA_IGOp(Msg_RIMP2G_Mem, NBlocks, 1, 'min')
c$$$C
c$$$         If ( NILimit .le. 0) then
c$$$            NILimit = NBlocks
c$$$         Else
c$$$            NILimit = Min( NILimit, NBlocks)
c$$$         EndIf
C
C              Sizes for this batch
C
               N1 = IntRange(2, 1, ISpin, BraKet, Class)
     $              - IntRange(1, 1, ISpin, BraKet, Class) + 1
               N2 = IntRange(2, 2, ISpin, BraKet, Class)
     $              - IntRange(1, 2, ISpin, BraKet, Class) + 1
C
C              If there is real work to do, do it!
C
               If ( N1 * N2 .gt. 0) then
C
C                 Informational stuff
C
                  NPass = NPass + 1
                  If ( PrPrgRpt ) then
                     If ( BraKet .eq. 1) then
                        Write ( String1, 9100) 
     $                       IntRange(1, 1, ISpin, BraKet, Class),
     $                       IntRange(2, 1, ISpin, BraKet, Class),
     $                       IntRange(1, 2, ISpin, BraKet, Class),
     $                       IntRange(2, 2, ISpin, BraKet, Class),
     $                       SpinLabel
                     Else
                        Write ( String1, 9105) 
     $                       IntRange(1, 1, ISpin, BraKet, Class),
     $                       IntRange(2, 1, ISpin, BraKet, Class),
     $                       IntRange(1, 2, ISpin, BraKet, Class),
     $                       IntRange(2, 2, ISpin, BraKet, Class),
     $                       SpinLabel
                     EndIf
                     Call Util_Sqz_Spc(String1)
                     String1 = '   '//String1
                     Call Progress_Report(LuOut,
     $                    String1(:Inp_StrLen(String1)) )
                  EndIF
 9100             Format('Beginning (', I5, ':', I5, ', ', I5, ':',
     $                 I5, '| fit ) spin ', A )
 9105             Format('Beginning ( fit |', I5, ':', I5, ', ',
     $                 I5, ':', I5, ') spin ', A )
C
C                 Create output GA
C
                  If ( .NOT. GA_Create_One_By_Shells( MT_Dbl, N1*N2,
     $                 FitBas, 'Output integral matrix', G_Out) )
     $                 Call ErrQuit(
     $                 'RIMP2_Driver_XF: Failed to allocate g_out',
     $                 ISpin, GA_ERR)
C
C                 Do the transformation.  Note that for the S
C                 approximation, there is an arbitrary choice to make
C                 the bra side the ERI and the ket side the overlap.
C
C                 For assymetric approximations, we transform the
C                 fitting basis index of the ket side only.
C
                  XFFit = 'YES'
                  If ( (.NOT. BraKetSame) .AND. BraKet .eq. 1)
     $                 XFFit = 'NO'
C
                  If ( (Approx .eq. 'SVS')
     $            .OR. (Approx .eq. 'S' .and. BraKet .eq. 2) ) then
                     If ( .NOT. XForm_3cIntsG(RTDB, Geom, 'OVERLAP',
     $                    AOBas, FitBas,
     $                    IntRange(1, 1, ISpin, BraKet, Class),
     $                    IntRange(2, 1, ISpin, BraKet, Class), 
     $                    IntRange(1, 2, ISpin, BraKet, Class),
     $                    IntRange(2, 2, ISpin, BraKet, Class), 
     $                    G_C(ISpin), XFFit, G_V, 'DISK', G_Out,
     $                    D_Int(ISpin, BraKet, Class), Thresh ) )
     $                    Call ErrQuit(
     $                    'RIMP2_Driver_XF: transformation failed', 0,
     &       CALC_ERR)
                  ElseIf ( (Approx .eq. 'V')
     $                .OR. (Approx .eq. 'S' .and. BraKet .eq. 1) ) then
                     If ( .NOT. XForm_3cIntsG(RTDB, Geom, 'ERI',
     $                    AOBas, FitBas,
     $                    IntRange(1, 1, ISpin, BraKet, Class),
     $                    IntRange(2, 1, ISpin, BraKet, Class), 
     $                    IntRange(1, 2, ISpin, BraKet, Class),
     $                    IntRange(2, 2, ISpin, BraKet, Class), 
     $                    G_C(ISpin), XFFit, G_V, 'DISK', G_Out,
     $                    D_Int(ISpin, BraKet, Class), Thresh ) )
     $                    Call ErrQuit(
     $                    'RIMP2_Driver_XF: transformation failed', 0,
     &       CALC_ERR)
                  Else
                     Write (LuErr, *) 'Unknown approximation ', Approx
                     Call ErrQuit(
     $                    'RIMP2_Driver_XF: Unrecognized approximation',
     $                    0, INPUT_ERR)
                  EndIf
C
C                 Clean up memory and disk
C
                  If ( .NOT. GA_Destroy( G_Out ) ) Call ErrQuit(
     $                 'RIMP2_Driver_XF: Can''t free global array',
     $                 ISpin, GA_ERR)
C
               EndIf            ! Transformation must be done
            EndDo               ! ISpin
         EndDo                  ! BraKet
      EndDo                     ! Class
C
C     Cleanup after the whole transformation is complete
C
      Status = .TRUE.
      Status = Status .AND. GA_Destroy( G_V )
C
      If ( .NOT. Status ) Call ErrQuit(
     $   'Rimp2_driver_xf:  Unable to free GAs ', 0, GA_ERR)
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


