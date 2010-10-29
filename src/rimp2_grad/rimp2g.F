      Logical Function RIMP2G( RTDB )
C$Id$
      Implicit NONE
#include "errquit.fh"
      Integer RTDB
C
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "pstat.fh"
#include "rimp2g_timers.fh"
#include "numerical_constants.fh"
#include "util.fh"
#include "dra.fh"
#include "inp.fh"
C
C EXTERNAL ROUTINES
      Logical  Int_Normalize, RIMP2_Driver_XF, RIMP2_Driver_E,
     $     Corr_Mk_Ref
      External Int_Normalize, RIMP2_Driver_XF, RIMP2_Driver_E, 
     $     Corr_Mk_Ref
C
C CONSTANTS
      Integer MaxSpin
      Parameter (MaxSpin = 2)
C
      Integer MinChunk
      Parameter (MinChunk = 64)
C
      Integer Kilobyte, Megabyte
      Parameter (Kilobyte = 1024, Megabyte = 1024*1024)
C
      Integer MxTBatch, MxEBatch(2)
C
      Logical DoDens, DoGrad, BraKetSame
c$$$      Parameter (DoGrad = .TRUE., BraKetSame = .FALSE.)
C
      Character*(*) Def_Int, Def_MOVec, Def_Gam3
      Parameter (Def_Int = '3ci', Def_MOVec = 'movecs',
     $     Def_Gam3 = 'gam')
      Character*(*) Geom_Name, AOBasName, FitBasName
      Parameter (Geom_Name = 'geometry', AOBasName = 'ao basis',
     $   FitBasName = 'ri-mp2 basis')
C
      Character*(*) My_Ctx
      Parameter (My_Ctx = 'mp2')
C
C LOCAL VARIABLES
      Integer AOBas, FitBas, Bases(2), Geom
      Integer NAOBF, TopSpin, NFit, RefNAO, NFrzOcc, NFrzVir
      Integer ALo(MaxSpin), AHi(MaxSpin), ILo(MaxSpin), IHi(MaxSpin)
      Integer NFrzO(MaxSpin), NAct(MaxSpin), NVir(MaxSpin),
     $     NFrzV(MaxSpin)
      Integer I, P, Me
      Integer H_Contrib, Contrib_Idx
      Integer NCent, H_Grad, I_Grad
      Double precision Thresh, CorEner, RefEner, TotEner, MinEVal
      Character*(80) FNInt, Approx, FNGam3
      Character*(255) AOName, AONameT, FitName, FitNameT, JobTitle
      Character*255 String1
      Logical Status, KeepFrz
      Character*(1) Mem_Opt
      Integer G_P2(MaxSpin), g_W2(MaxSpin)
C
      Double precision Tol2e_Fock
C
C     Handles and indices to for the "correlated" basis
C
      Integer H_Occ, H_EVal
      Integer I_Occ, I_EVal
      Integer G_AOXF(MaxSpin)
      Integer NMO(MaxSpin), MxNCorBF
C
C     "Handles" for the integral file(s).  Indices are:
C         1 = spin (alpha, beta)
C         2 = bra, ket (may be distinct)
C         3 = integral group (pa, ij)
C
      Integer D_Int(MaxSpin, 2, 2), D_Gam(MaxSpin)
      Character*80 FN_Int(MaxSpin, 2, 2), FN_Gam(MaxSpin)
C
C     Range specifications for transformations
C        1   = lower/upper limit
C        2   = index 1/2
C        3+ = as in D_Int, D_Gam 1+
C
      Integer IntRange(2, 2, MaxSpin, 2, 2), GamRange(2, 2, MaxSpin)
C
C     Print controls
C
      Logical PrInfo, PrPrgRpt
C
      Double precision MaxDASiz, TotDASiz
      Integer NrDA
C
C     ******************
C     * Initialization *
C     ******************
C
      Tol2e_Fock = 1.0d-9
C
      Status = .TRUE.
C
      Call RIMP2G_PStat_Init
      Call PStat_On( PS_Overall)
C
      Me = GA_NodeID()
C
C     Get necessary info from database
C
      Call Util_Print_Push
      Call Util_Print_RTDB_Load(RTDB, My_Ctx)
C
      If ( .NOT. RTDB_Get(RTDB, 'mp2:rimp2 density', MT_Log, 1,
     $     DoDens) ) DoDens = .FALSE.
      If ( .NOT. RTDB_Get(RTDB, 'mp2:rimp2 gradient', MT_Log, 1,
     $     DoGrad) ) DoGrad = .FALSE.
C
C     Current gradients require breaking bra-ket symmetry, but
C     anything short of that will still work.
C
      If ( .NOT. RTDB_Get(RTDB, 'mp2:rimp2 braketsame', MT_Log, 1,
     $     BraKetSame) ) BraKetSame = .NOT. DoGrad
      Write (LuOut, *) 'DoDens, DoGrad, BraKetSame = ',
     $     DoDens, DoGrad, BraKetSame
C
      If (.NOT. RTDB_CGet(RTDB, 'title', 1, JobTitle)) JobTitle = ' '
C
      PrInfo = Util_Print('information', Print_Low) .AND. (Me .eq. 0)
      PrPrgRpt = Util_Print('progress reports', Print_Default)
C
C     *************************************************
C     * Start creating & initializing data structures *
C     *************************************************
C
      Status = Status .AND. Geom_Create( Geom, Geom_Name)
      Status = Status .AND. Geom_RTDB_Load( RTDB, Geom, Geom_Name)
C
      Status = Status .AND. Bas_Create( AOBas, AOBasName)
      Status = Status .AND. Bas_RTDB_Load( RTDB, Geom, AOBas, AOBasName)
      Status = Status .AND. Int_Normalize( RTDB,AOBas)
C
      Status = Status .AND. Bas_Create( FitBas, FitBasName)
      Status = Status .AND. Bas_RTDB_Load( RTDB, Geom, FitBas,
     $   FitBasName)
      Status = Status .AND. Int_Normalize( RTDB, FitBas)
C
C     Get other controls -- threshold for neglect of integrals, etc.
C
      Call RIMP2G_Get_Param( RTDB, Thresh, MxTBatch, Mem_Opt, MxEBatch,
     $   Approx, MinEVal)
C
C     Figure out file names.  Munging to get the full name is done later!
C
      if (.not. RTDB_CGet(rtdb, 'mp2:mo 3-center integral file',
     $     1, FNInt)) FNInt = Def_Int
      if (.not. RTDB_CGet(rtdb, 'mp2:mo 3-center gamma file',
     $     1, FNGam3)) FNGam3 = Def_Gam3
C
C     Get some info about the system & basis sets for later use
C
      If (.NOT. Geom_NCent( Geom, NCent) ) Call ErrQuit(
     $     'RIMP2g: geom_ncent failed', Geom, GEOM_ERR)
C
      Status = Status .AND. Bas_NumBF( AOBas, NAOBF)
      Status = Status .AND. Bas_NumBF( FitBas, NFit)
      Status = Status .AND. Bas_Name(AOBas, AOName, AONameT)
      Status = Status .AND. Bas_Name(FitBas, FitName, FitNameT)
C
C     *********************
C     * Get SCF reference *
C     *********************
C     For the time being, we must keep frozen orbitals if we're doing
C     gradients, but otherwise we don't _need_ them, so we drop them.
C
      KeepFrz = DoDens .OR. DoGrad
      If ( .NOT. Corr_Mk_Ref(RTDB, My_Ctx, Geom, AOBas, .TRUE., 
     $     RefNAO, TopSpin, NMO, G_AOXF, MinChunk, MinChunk, MaxSpin,
     $     H_Occ, H_EVal, NFrzOcc, NFrzVir, KeepFrz) ) Call ErrQuit(
     $     'RIMP2G: Failed to get reference', 0, INPUT_ERR)
C
      String1 = 'RIMP2G: density & gradients support frozen core only, '
     $     // 'not frozen virtuals.'
      If ( (DoDens .OR. DoGrad) .AND. (NFrzVir .ne. 0) ) Call ErrQuit(
     $     String1(:Inp_StrLen(String1)), NFrzVir, UNKNOWN_ERR)
C
      If (RefNAO .ne. NAOBF) Call ErrQuit(
     $     'RIMP2G: Incompatible reference AO basis size', RefNAO-NAOBF,
     &       BASIS_ERR)
      MxNCorBF = NAOBF
      If ( .NOT. MA_Get_Index(H_Occ, I_Occ) ) Call ErrQuit(
     $     'RIMP2G: Failed to get index for occ vector', H_Occ, MA_ERR)
      If ( .NOT. MA_Get_Index(H_EVal, I_EVal) ) Call ErrQuit(
     $     'RIMP2G: Failed to get index for eval vector', H_EVal,
     &       MA_ERR)
C
C     Count number of frozen/active occupied/virtual.
C
C     Corr_Mk_Ref does not support separate freeze directives by
C     spin, but in the rest of the code, we can easily handle that
C     to make the transition here.
C
C     IMPORTANT NOTE: If .NOT. KeepFrz, then the frozen orbitals have been
C     removed from the occupation vector (and other quantities), but
C     NFrzOcc and NFrzVir can be non-zero because they are meant for
C     reporting to the user.  For the time being, we can setup the new
C     NFrzO and NFrzV arrays to the correct values for the code, while
C     NFrzOcc and NFrzVir continue to carry the original information
C     which can be reported to the user.
C
      Do I = 1, TopSpin
         If ( KeepFrz) then
            NFrzO(i) = NFrzOcc
            NFrzV(i) = NFrzVir
         Else
            NFrzO(i) = 0
            NFrzV(i) = 0
         EndIf
         NAct(I) = 0
         NVir(I) = 0
C
         Do P = 1, NMO(I)
            If ( Dbl_MB( I_Occ+MxNCorBF*(I-1)+P-1 ) .ne. 0.0d0 ) then
               NAct(I) = NAct(I) + 1
            Else
               NVir(I) = NVir(I) + 1
            EndIf
         EndDo
C
C        If the Occ array includes the frozen orbitals, we've overcounted
C        for NAct and NVir.
C
         If ( KeepFrz ) then
            NAct(I) = NAct(I) - NFrzO(i)
            NVir(I) = NVir(I) - NFrzV(i)
         EndIf
C
         Call Occ_Ranges( NMO(I), Dbl_MB(I_Occ + MxNCorBF * (I-1) ),
     $      ILo(I), IHi(I), ALo(I), AHi(I) )
      EndDo
C
      If ( .NOT. MA_Free_Heap( H_Occ) ) Call ErrQuit(
     $   'RIMP2G: Unable to free local array H_Occ', H_Occ, MA_ERR)
C
C     ***********************
C     * Print module header *
C     ***********************
C
      If ( PrInfo) then
         Write (LuOut, '(//)')
         Call Util_Print_Centered(LuOut,
     $      'NWChem RI-MP2 Energy/Gradient Module', 40, .TRUE.)
         If ( JobTitle .ne. ' ') Call Util_Print_Centered( LuOut,
     $      JobTitle, 40, .FALSE.)
         Write (LuOut, '(//)')
      EndIf
C
C     ******************************
C     * Prepare file names & sizes *
C     ******************************
C     We have enough info to do this now, though most of it doesn't
C     get used for quite some time.  However PrInfo will print
C     largest and total, so the user has some idea what to expect.
C
      Call RIMP2_Prep_DRA(MaxSpin, TopSpin, NFrzO, NAct, NVir,
     $     NFrzV, NFit, DoDens .OR. DoGrad, BraKetSame, NrDA, MaxDASiz,
     $     TotDASiz, FNInt(:Inp_StrLen(FNInt)),
     $     FNGam3(:Inp_StrLen(FNGam3)), FN_Int, FN_Gam, IntRange,
     $     GamRange)
C
C     **********************************************
C     * Print out calculation information for user *
C     **********************************************
C
      If ( PrInfo) then
C
c$$$         Write (LuOut, '(A)') 
c$$$         Call Banner(LuOut, 'RI-MP2 Calculation Information', '-',
c$$$     $      .TRUE., .TRUE., .FALSE.)
         If ( DoGrad ) then
            Write (LuOut, 8505) 'Calculation Type:', 'Energy + Gradient'
         ElseIf ( DoDens .AND. .NOT. DoGrad ) then
            Write (LuOut, 8505) 'Calculation Type:', 'Energy + Density'
         Else
            Write (LuOut, 8505) 'Calculation Type:', 'Energy only'
         EndIf
         Write (LuOut, 8520) 'Nr. of AO basis functions:', NAOBF
         Write (LuOut, 8520) 'Nr. of fitting basis functions:', NFit
         Write (LuOut, 8520) 'Nr. of MO spin cases:', TopSpin
         Write (LuOut, 8520) 'Nr. of MOs to be correlated:',
     $      (NMO(i), i = 1, TopSpin)
         Write (LuOut, 8520) 'Nr. of frozen occupied MOs',
     $        (NFrzOcc, i = 1, TopSpin)
         Write (LuOut, 8520) 'Nr. of occupied MOs correlated:',
     $      (NAct(i), i = 1, TopSpin)
         Write (LuOut, 8520) 'Nr. of virtual MOs correlated:',
     $      (NVir(i), i = 1, TopSpin)
         Write (LuOut, 8520) 'Nr. of frozen virtual MOs',
     $        (NFrzVir, i = 1, TopSpin)
         Write (LuOut, *)
C
         If ( BraKetSame ) then
            Write (LuOut, 8505) 'RI integral approximation:',
     $           Approx(:Inp_StrLen(Approx)), '(symmetric)'
         Else
            Write (LuOut, 8505) 'RI integral approximation:',
     $           Approx(:Inp_StrLen(Approx)), '(asymmetric)'
         EndIf
         Write (LuOut, 8510) 'Threshold for Schwarz screening:', Thresh
         Write (LuOut, 8510)
     $      'Smallest eigenvalue retained in fit:', MinEVal
         If ( MxTBatch .gt. 0) then
            Write (LuOut, 8520) 'Max. batch size for transformation:',
     $         MxTBatch
         EndIf
         If ( MxEBatch(1) .gt. 0 .OR. MxEBatch(2) .gt. 0) then
            Write (LuOut, 8520) 'Max. batch size for energy:',
     $         MxEBatch
         EndIf
         Write (LuOut, *)
C
C        Quick sanity check
C
         Do I = 1, TopSpin
            If ( NAct(i) .eq. 0 ) Call ErrQuit(
     $           'rimp2g: no active orbitals to correlate', i,
     &       BASIS_ERR)
            If ( NVir(i) .eq. 0 ) Call ErrQuit(
     $           'rimp2g: no virtual orbitals to correlate into', i,
     &       BASIS_ERR)
         EndDo
C
C        File information
C
         If ( DoDens .OR. DoGrad ) then
            Write (LuOut, 8505) 'Base file suffixes:', 
     $           FNInt(:Inp_StrLen(FNInt)), FNGam3(:Inp_StrLen(FNGam3))
         Else
            Write (LuOut, 8505) 'Base file suffixes:', 
     $           FNInt(:Inp_StrLen(FNInt))
         EndIf
         Write (LuOut, 8530) 'Largest file (MB):', MaxDASiz/Megabyte
         Write (LuOut, 8530) 'Total disk space (MB):', TotDASiz/Megabyte
         Write (LuOut, *)
C
C        If user requests, this is the place to dump out a lot of
C        additional information
C
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom))
     $         call errquit('RIMP2G: geom_print', 0, GEOM_ERR)
            write(LuOut,*)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(AOBas)) 
     $           call errquit('RIMP2G: bas_print', 0, BASIS_ERR)
            write(LuOut,*)
            if (.not. bas_print(FitBas)) 
     $           call errquit('RIMP2G: bas_print', 0, BASIS_ERR)
            write(LuOut,*)
         endif
         if (util_print('geombas', print_debug)) then
            if (.not. gbs_map_print(AOBas)) 
     $           call errquit('RIMP2G:gmp', 0, UNKNOWN_ERR)
            write(LuOut,*)
            if (.not. gbs_map_print(FitBas)) 
     $           call errquit('RIMP2G:gmp', 0, UNKNOWN_ERR)
            write(LuOut,*)
         endif
         call util_flush(LuOut)
      EndIf
C
 8500 Format(1X, A35, 1X, '''', A, '''' )
 8505 Format(1X, A35, 4(2X, A), (/36X, 4(2X, A)) )
 8510 Format(1X, A35, 4(1X, 1PE10.3), (/36X, 4(1X, 1PE10.3)) )
 8520 Format(1X, A35, 4(1X, I10), (/36X, 4(1X,I10)) )
 8530 Format(1X, A35, 4(1X, F10.3), (/36X, 4(1X, F10.3)) )
C
 9050 Format(1X, 'The following mapping relates the ''correlated ',
     $       'basis'' to the'/
     $   1X, 'original SCF basis.')
 9051 Format(1X, 'This correlated calculation is based on an SCF ',
     $       'reference in which'/
     $   1X, 'orbitals have been frozen or rearranged relative to ',
     $        'the result of the'/
     $   1X, 'SCF calculation.  From this point on, the orbitals ',
     $       'are labeled in the'/
     $   1X, '''correlated basis''.'/)
C
C     Initialize DAs. Size limits computed earlier; let system
C     decide how much memory to give it.
C
      If (DRA_Init(NrDA, MaxDASiz, TotDASiz, -1) .ne. 0)
     $   Call ErrQuit('RIMP2G: DRA_Init failed', 0, UNKNOWN_ERR)
C
C     Create all necessary DRAs (assumes no per-DRA memory cost)
C
      Call RIMP2_Create_DRA(MaxSpin, TopSpin, NFit, BraKetSame, 
     $     FN_Int, FN_Gam, IntRange, GamRange, D_Int, D_Gam)
C
C     Initialize integrals
C
      Bases(1) = AOBas
      Bases(2) = FitBas
C
      Call Int_Init(RTDB, 2, Bases)
C
C
      Status = Status .AND.
     $   RIMP2_Driver_XF( RTDB, Geom, AOBas, FitBas, TopSpin,
     $     MaxSpin, D_Int, NFrzO, NAct, NVir, NFrzV, KeepFrz,
     $     ALo, AHi, ILo, IHi, MxTBatch, Thresh, MinEVal, Approx,
     $     G_AOXF, BraKetSame, DoDens .OR. DoGrad,
     $     IntRange)
C
C     If we're not doing gradients, we can clean up some stuff now
C
      If ( .NOT. (DoDens .OR. DoGrad) ) then
         Do I = TopSpin, 1, -1
            If ( .NOT. GA_Destroy( G_AOXF(I)) ) Call ErrQuit(
     $           'RIMP2G: can''t destroy AO transformation matrix', I,
     &       GA_ERR)
         EndDo
      EndIf
C
C     If we're doing gradients, we'll have to re-initialize again
C     to do derivative integrals
C
      Call Int_Terminate
C
C     Shift ranges of virtual and occupieds into the _number_ of each
C
      Do I = 1, TopSpin
         AHi(I) = AHi(I) - ALo(I) + 1
         ALo(I) = 1
         IHi(I) = IHi(I) - ILo(I) + 1
         ILo(I) = 1
      EndDo
C
C     Create matrix to hold spin case breakdown of E(2)
C
      Status = Status .AND. MA_Push_Get( MT_Dbl, TopSpin*TopSpin,
     $   'E(2) spin case breakdown', H_Contrib, Contrib_Idx)
C
C     **********************
C     * Energy calculation *
C     **********************
C     Print & store results before proceeding to gradient to give
C     something in case of later failure.
C
      Status = .TRUE.
      Status = Status .AND.
     $   RIMP2_Driver_E( MaxSpin, BraKetSame, D_Int, 
     $     D_Gam, TopSpin, NMO, NFrzO, NAct, NVir, NFit,
     $     FitBas, Dbl_MB(I_Eval), MxNCorBF, Dbl_MB(Contrib_Idx),
     $     CorEner, Mem_Opt, MxEBatch(1), MxEBatch(2),
     $     DoDens .OR. DoGrad, G_P2)
C
C     Also get the SCF energy so we can report the total energy
C
      Status = RTDB_Get(RTDB, 'scf:energy', MT_Dbl, 1, RefEner)
      If ( .NOT. Status ) then
         RefEner = FP_Zero
         If ( Me .eq. 0) Write ( LuOut, *)
     $      'RIMP2G: Warning! Can''t get scf:energy from RTDB. ',
     $        'Continuing.'
      EndIf
C
      TotEner = RefEner + CorEner
C
      If ( Me .eq. 0) then
         Write (LuOut, '(A)') ' '
         Call Banner(LuOut, 'RI-MP2 energy', '-',
     $      .TRUE., .TRUE., .FALSE.)
         Write (LuOut, '(1X, A, F30.15)') 'Reference energy:    ',
     $      RefEner
         Write (LuOut, '(1X, A, F30.15)') 'Correlation energy:  ',
     $      CorEner
         Write (LuOut, '(1X, A, F30.15)') 'Total RI-MP2 energy: ',
     $      TotEner
      EndIf
C
C     Write our energy out to the database
C
      Status = .TRUE.
      Status = Status .AND. RTDB_Put( RTDB,
     $   'rimp2_grad:correlation energy', MT_Dbl, 1, CorEner)
      Status = Status .AND. RTDB_Put( RTDB, 'rimp2_grad:energy',
     $   MT_Dbl, 1, TotEner)
      If (.NOT. Status .AND. Me .eq. 0) Write (LuOut, *)
     $     'RIMP2G: WARNING! Can''t write energy to RTDB. Continuing.'
C
C     Clean up some memory
C
      Status = .TRUE.
      Status = Status .AND. MA_Pop_Stack( H_Contrib)
c$$$      Status = Status .AND.  MA_Free_Heap(H_Eval)
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2G: Unable to destroy local arrays', 0, MA_ERR)
C
C     ***********
C     * Density *
C     ***********
C
      If ( DoDens .OR. DoGrad ) then
         Call RIMP2_Driver_P(RTDB, Geom, AOBas, NAOBF,
     $        MaxSpin, TopSpin, NMO, 
     $        NFrzO, NAct, NVir,
     $        FitBas, g_P2, g_W2,
     $        Dbl_MB(I_Eval), MxNCorBF, g_AOXF, Tol2e_Fock,
     $        JobTitle(:Inp_StrLen(JobTitle)), D_Int, D_Gam )
      EndIf
C
C     ************
C     * Gradient *
C     ************
C
      If ( DoGrad ) then
C
         Call RIMP2_Driver_G(RTDB, Geom, AOBas, NAOBF,
     $        MaxSpin, TopSpin, NMO, 
     $        NFrzO, NAct, NVir,
     $        FitBas, BraKetSame, g_P2, g_W2,
     $        Dbl_MB(I_Eval), MxNCorBF, g_AOXF, Tol2e_Fock,
     $        JobTitle(:Inp_StrLen(JobTitle)), D_Int, D_Gam, NCent)
C
      EndIf
C
C     ***********
C     * Cleanup *
C     ***********
C
      If ( DoDens .OR. DoGrad ) then
         Do I = TopSpin, 1, -1
            If ( .NOT. GA_Destroy( G_AOXF(I)) ) Call ErrQuit(
     $           'RIMP2G: can''t free AO transformation matrix', I,
     &       GA_ERR)
            If ( .NOT. GA_Destroy( G_P2(I)) ) Call ErrQuit(
     $           'RIMP2G: can''t free P(2)', I, GA_ERR)
            If ( .NOT. GA_Destroy( G_W2(I)) ) Call ErrQuit(
     $           'RIMP2G: can''t free W(2)', I, GA_ERR)
         EndDo
      EndIf
C
      If ( DoGrad ) then
         Do I = TopSpin, 1, -1
         EndDo
      EndIf
C
C     Close, delete, and terminate DRAs
C
      Call RIMP2_Delete_DRA(MaxSpin, TopSpin, BraKetSame, D_Int, D_Gam)
      If ( DRA_Terminate() .ne. 0) Call ErrQuit(
     $   'RIMP2G: DRA_Terminate failed', 0, CALC_ERR)
C
      Status = .TRUE.
      Status = Status .AND.  MA_Free_Heap(H_Eval)
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2G: Unable to destroy local arrays', 0, MA_ERR)
C
      Status = .TRUE.
      Status = Status .AND. Geom_Destroy( Geom)
      Status = Status .AND. Bas_Destroy( FitBas)
      Status = Status .AND. Bas_Destroy( AOBas)
      If (.NOT. Status) Call ErrQuit(
     $   'RIMP2G: Failed to destroy geometry and  basis sets', 0,
     &       GEOM_ERR)
C
      Call PStat_Off( PS_Overall)
      Call RIMP2G_PStat_Print
      Call RIMP2G_PStat_Term
C
      RIMP2G = .TRUE.
      Return
      End
