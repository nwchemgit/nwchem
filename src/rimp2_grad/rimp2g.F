      Logical Function RIMP2G( RTDB )
C$Id: rimp2g.F,v 1.6 1998-07-13 02:38:10 bernhold Exp $
      Implicit NONE
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
      Logical  Int_Normalize, RIMP2_Driver_XF, RIMP2_Driver_E
      External Int_Normalize, RIMP2_Driver_XF, RIMP2_Driver_E
C
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
      Logical DoGrad, BraKetSame
      Parameter (DoGrad = .TRUE., BraKetSame = .TRUE.)
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
      Integer AOBas, FitBas, Bases(2), Geom
      Integer NAOBF, TopSpin, NFit, RefNAO, NFrzOcc, NFrzVir
      Integer ALo(MaxSpin), AHi(MaxSpin), ILo(MaxSpin), IHi(MaxSpin)
      Integer NFrz(MaxSpin), NAct(MaxSpin), NVir(MaxSpin)
      Integer I, Me
      Integer H_Contrib, Contrib_Idx
      Double precision Thresh, CorEner, RefEner, TotEner, MinEVal
      Character*(80) FNInt, Approx, FNGam3
      Character*(255) AOName, AONameT, FitName, FitNameT, JobTitle
      Logical Status
      Character*(1) Mem_Opt
      Integer G_P2(MaxSpin)
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
C     Print controls
C
      Logical PrInfo
C
      Integer MaxDASiz, TotDASiz, Size, SizeB, NrDA
C
      Data NFrz/0, 0/
C
C     Initialize this module
C
      Status = .TRUE.
C
C     For debugging purposes
C
      Call IFill(MaxSpin*2*2, -99352, D_Int, 1)
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
      PrInfo = Util_Print('information', Print_Low)
C
      If (.NOT. RTDB_CGet(RTDB, 'title', 1, JobTitle)) JobTitle = ' '
C
C     *****************
C     * Module header *
C     *****************
C
      If ( Me .eq. 0 .AND. PrInfo) then
C
         Call Util_Print_Centered(LuOut,
     $      'NWChem RI-MP2 Energy/Gradient Module', 40, .TRUE.)
         If ( JobTitle .ne. ' ') Call Util_Print_Centered( LuOut,
     $      JobTitle, 40, .FALSE.)
         Write (LuOut, '(//)')
      EndIf
C
      Status = Status .AND. Geom_Create( Geom, Geom_Name)
      Status = Status .AND. Geom_RTDB_Load( RTDB, Geom, Geom_Name)
C
      Status = Status .AND. Bas_Create( AOBas, AOBasName)
      Status = Status .AND. Bas_RTDB_Load( RTDB, Geom, AOBas, AOBasName)
      Status = Status .AND. Int_Normalize( AOBas)
C
      Status = Status .AND. Bas_Create( FitBas, FitBasName)
      Status = Status .AND. Bas_RTDB_Load( RTDB, Geom, FitBas,
     $   FitBasName)
      Status = Status .AND. Int_Normalize( FitBas)
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
C     Get some info about the basis sets for later use
C
      Status = Status .AND. Bas_NumBF( AOBas, NAOBF)
      Status = Status .AND. Bas_NumBF( FitBas, NFit)
      Status = Status .AND. Bas_Name(AOBas, AOName, AONameT)
      Status = Status .AND. Bas_Name(FitBas, FitName, FitNameT)
C
      Call Corr_Mk_Ref(RTDB, My_Ctx, Geom, RefNAO, TopSpin, NMO,
     $     G_AOXF, MinChunk, MinChunk, MaxSpin, H_Occ, H_EVal,
     $     NFrzOcc, NFrzVir, .False.)
C
      If (RefNAO .ne. NAOBF) Call ErrQuit(
     $     'RIMP2G: Incompatible reference AO basis size', RefNAO-NAOBF)
      MxNCorBF = NAOBF
      If ( .NOT. MA_Get_Index(H_Occ, I_Occ) ) Call ErrQuit(
     $     'RIMP2G: Failed to get index for occ vector', H_Occ)
      If ( .NOT. MA_Get_Index(H_EVal, I_EVal) ) Call ErrQuit(
     $     'RIMP2G: Failed to get index for eval vector', H_EVal)
C
C     Determine ranges to transform.
C
      Do I = 1, TopSpin
         Call Occ_Ranges( NMO(I), Dbl_MB(I_Occ + MxNCorBF * (I-1) ),
     $      ILo(I), IHi(I), ALo(I), AHi(I) )
      EndDo
C
      If ( .NOT. MA_Free_Heap( H_Occ) ) Call ErrQuit(
     $   'RIMP2G: Unable to free local array H_Occ', H_Occ)
C
C     Prepare DRA sizes and names
C
      Call RIMP2_Prep_DRA(TopSpin, ILo, IHi, ALo, AHi, NFit,
     $     DoGrad, BraKetSame, NrDA, MaxDASiz, TotDASiz,
     $     FNInt(:Inp_StrLen(FNInt)), FNGam3(:Inp_StrLen(FNGam3)),
     $     FN_Int, FN_Gam)
C
C     **********************************************
C     * Print out calculation information for user *
C     **********************************************
C
      If ( Me .eq. 0 .AND. PrInfo) then
C
         Write (LuOut, '(A)') 
         Call Banner(LuOut, 'RI-MP2 Calculation Information', '-',
     $      .TRUE., .TRUE., .FALSE.)
         If ( DoGrad ) then
            Write (LuOut, 8505) 'Calculation Type:', 'Energy + Gradient'
         Else
            Write (LuOut, 8505) 'Calculation Type:', 'Energy only'
         EndIf
         Write (LuOut, 8520) 'Nr. of AO basis functions:', NAOBF
         Write (LuOut, 8520) 'Nr. of fitting basis functions:', NFit
         Write (LuOut, 8520) 'Nr. of MO spin cases:', TopSpin
         Write (LuOut, 8520) 'Nr. of MOs to be correlated:',
     $      (NMO(i), i = 1, TopSpin)
         Write (LuOut, 8520) 'Nr. of occupied MOs correlated:',
     $      (IHi(i)-ILo(i)+1, i = 1, TopSpin)
         Write (LuOut, 8520) 'Nr. of virtual MOs correlated:',
     $      (AHi(i)-ALo(i)+1, i = 1, TopSpin)
         Write (LuOut, *)
C
         If ( BraKetSame ) then
            Write (LuOut, 8505) 'RI integral approximation:',
     $           Approx(:Inp_StrLen(Approx)), '(symmetric)'
         Else
            Write (LuOut, 8505) 'RI integral approximation:',
     $           Approx(:Inp_StrLen(Approx)), '(assymetric)'
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
C        File information
C
         If ( DoGrad ) then
            Write (LuOut, 8505) 'Base file suffixes:', 
     $           FNInt(:Inp_StrLen(FNInt)), FNGam3(:Inp_StrLen(FNGam3))
         Else
            Write (LuOut, 8505) 'Base file suffixes:', 
     $           FNInt(:Inp_StrLen(FNInt))
         EndIf
         Write (LuOut, 8530) 'Largest file (MB):',
     $        Dble(MaxDASiz)/Megabyte
         Write (LuOut, 8530) 'Total disk space (MB):',
     $        Dble(TotDASiz)/Megabyte
         Write (LuOut, *)
C
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom))
     $         call errquit('RIMP2G: geom_print', 0)
            write(LuOut,*)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(AOBas)) 
     $           call errquit('RIMP2G: bas_print', 0)
            write(LuOut,*)
            if (.not. bas_print(FitBas)) 
     $           call errquit('RIMP2G: bas_print', 0)
            write(LuOut,*)
         endif
         if (util_print('geombas', print_debug)) then
            if (.not. gbs_map_print(AOBas)) 
     $           call errquit('RIMP2G:gmp', 0)
            write(LuOut,*)
            if (.not. gbs_map_print(FitBas)) 
     $           call errquit('RIMP2G:gmp', 0)
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
     $   Call ErrQuit('RIMP2G: DRA_Init failed', 0)
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
     $     MaxSpin, D_Int, ALo, AHi, ILo, IHi, MxTBatch, Thresh,
     $     MinEVal, Approx, G_AOXF, FN_Int, BraKetSame, DoGrad)
C
C     Finished with the AO transformation
C
      Status = .TRUE.
      Do I = TopSpin, 1, -1
         Status = Status .AND.  GA_Destroy( G_AOXF(I))
      EndDo
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2G: Unable to destroy AO transformation', 0)
C
C     Done computing integrals
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
C
         NAct(I) = IHi(I)
         NVir(I) = AHi(I)
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
     $   RIMP2_Driver_E( MaxSpin, BraKetSame, FN_Int, D_Int, FN_Gam,
     $     D_Gam, TopSpin, NMO, NFrzOcc, IHi, AHi, NFit,
     $     FitBas, Dbl_MB(I_Eval), MxNCorBF, Dbl_MB(Contrib_Idx),
     $     CorEner, Mem_Opt, MxEBatch(1), MxEBatch(2), DoGrad, G_P2)
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
      Status = Status .AND.  MA_Free_Heap(H_Eval)
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2G: Unable to destroy local arrays', 0)
C
C     ************
C     * Gradient *
C     ************
C
      If ( DoGrad ) then
         Call RIMP2_Driver_G(TopSpin, NFrzOcc, NAct, NVir, FitBas)
      EndIf
C
C     ***********
C     * Cleanup *
C     ***********
C
      If ( DRA_Terminate() .ne. 0) Call ErrQuit(
     $   'RIMP2G: DRA_Terminate failed', 0)
C
      Status = .TRUE.
      Status = Status .AND. Geom_Destroy( Geom)
      Status = Status .AND. Bas_Destroy( FitBas)
      Status = Status .AND. Bas_Destroy( AOBas)
      If (.NOT. Status) Call ErrQuit(
     $   'RIMP2G: Failed to destroy geometry and  basis sets', 0)
C
      Call PStat_Off( PS_Overall)
      Call RIMP2G_PStat_Print
      Call RIMP2G_PStat_Term
C
      RIMP2G = .TRUE.
      Return
      End
