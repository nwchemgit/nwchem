      Logical Function RIMP2( RTDB )
C$Id: rimp2.F,v 1.43 2000-02-04 01:45:12 mg201 Exp $
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
#include "rimp2_timers.fh"
#include "numerical_constants.fh"
#include "util.fh"
#ifdef RIMP2_USE_DA
#include "dra.fh"
#endif
C
      Logical MOVecs_Read_Header, MOVecs_Read, MOVecs_Gather,
     $   Int_Normalize
      External MOVecs_Read_Header, MOVecs_Read, MOVecs_Gather,
     $   Int_Normalize
C
      Integer LNBlnk
      External LNBlnk
C
      Logical RIMP2_XF_Driver, RIMP2_V_E2, SCF
      External RIMP2_XF_Driver, RIMP2_V_E2, SCF
C
      Integer MaxSpin
      Parameter (MaxSpin = 2)
C
      Integer MinChunk
      Parameter (MinChunk = 64)
C
      Integer MxTBatch, MxEBatch(2)
C
      Character*(*) Def_Int, Def_MOVec
      Parameter (Def_Int = 'mo3cint', Def_MOVec = 'movecs')
      Character*(*) Geom_Name, AOBasName, FitBasName
      Parameter (Geom_Name = 'geometry', AOBasName = 'ao basis',
     $   FitBasName = 'ri-mp2 basis')
C
      Character*(*) My_Ctx
      Parameter (My_Ctx = 'mp2')
C
      Integer AOBas, FitBas, Bases(2), Geom
      Integer NAOBF, TopSpin, RefSpin, NFit, RefNAO, NFrzOcc, NFrzVir
      Integer ALo(MaxSpin), AHi(MaxSpin), ILo(MaxSpin), IHi(MaxSpin)
      Integer RefALo(MaxSpin), RefAHi(MaxSpin), RefILo(MaxSpin),
     $   RefIHi(MaxSpin)
      Integer NMO(MaxSpin), RefNMO(MaxSpin),I, Me
      Integer H_Contrib, Contrib_Idx
      Double precision Thresh, CorEner, RefEner, TotEner, MinEVal
      Character*(80) Title, SCFBasNam, FNVec, FNInt, Approx
      character*20 scftype
      Character*(255) AOName, AONameT, FitName, FitNameT, JobTitle
      Logical Status
      Character*(1) Mem_Opt
C
C     Handles and indices to for the "correlated" basis
C
      Integer H_Occ, H_EVal, H_CorMap
      Integer I_Occ, I_EVal, I_CorMap
      Integer G_AOXF(MaxSpin)
      Integer NCorBF(MaxSpin), MxNCorBF
C
C     "Handles" for the integral file(s) (either logical unit numbers
C     or disk array handles)
C
      Integer D_Int(MaxSpin)
C
C     Print controls
C
      Logical PrInfo, PrWarn, PrRef
C
      Logical FullRef, FunnyRef
#ifdef RIMP2_USE_DA
C
      Integer MaxDASiz, TotDASiz, Size, SizeB
#endif
C
C     Stuff for SCF reference
C
      Double precision ThreshSCF
C
C     Initialize this module
C
      Status = .TRUE.
C
      Call RIMP2_PStat_Init
      Call PStat_On( PS_Overall)
C
      Me = GA_NodeID()
C
#ifndef RIMP2_USE_DA
C
C     Initialize integral file LUNs if not using DAs
C
      Do I = 1, MaxSpin
         D_Int(i) = 9 + i
      EndDo
#endif
C
C     Get necessary info from database
C
      Call Util_Print_Push
      Call Util_Print_RTDB_Load(RTDB, My_Ctx)
      PrInfo = Util_Print('information', Print_Low)
C
      call ecce_print_module_entry('mp2')
C
      If (.NOT. RTDB_CGet(RTDB, 'title', 1, JobTitle)) JobTitle = ' '
C
C     ************************************
C     * Check for adequate SCF reference *
C     ************************************
c
      if (.not. rtdb_get(rtdb,'scf:thresh',mt_dbl,1,threshscf))
     $     threshscf = 1.0d-4
      if (threshscf .gt. 1.0d-6) then
         threshscf = 1.0d-6
         if (.not. rtdb_put(rtdb,'scf:thresh',mt_dbl,1,threshscf))
     $        call errquit('mp2_grad: failed setting scf:thresh',0)
         status = rtdb_delete(rtdb,'scf:converged') ! Status irrelevant
      endif
c
      if (.not. scf(rtdb)) then
         rimp2 = .false.
         call ecce_print_module_exit('mp2', 'failed')
         Call Util_Print_pop
         return
      endif
C
C     *****************
C     * Module header *
C     *****************
C
      If ( Me .eq. 0 .AND. PrInfo) then
C
         Call Util_Print_Centered(LuOut, 'NWChem RI-MP2 Module', 40,
     $   .TRUE.)
         If ( JobTitle .ne. ' ') Call Util_Print_Centered( LuOut,
     $      JobTitle, 40, .FALSE.)
         Write (LuOut, '(//)')
      EndIf
C
      If ( .NOT. Geom_Create( Geom, Geom_Name) )
     $   Call ErrQuit('rimp2: geom_create failed', Geom)
      If ( .NOT. Geom_RTDB_Load( RTDB, Geom, Geom_Name) )
     $   Call ErrQuit('rimp2: geom_rtdb_load failed', Geom)
C
      If ( .NOT. Bas_Create( AOBas, AOBasName) )
     $   Call ErrQuit('rimp2: aobas bas_create failed', AOBas)
      If ( .NOT. Bas_RTDB_Load( RTDB, Geom, AOBas, AOBasName) )
     $   Call ErrQuit('rimp2: aobas bas_rtdb_load failed', AOBas)
      If ( .NOT. Int_Normalize( AOBas, RTDB) )
     $   Call ErrQuit('rimp2: aobas int_normalize failed', AOBas)
C
      If ( .NOT. Bas_Create( FitBas, FitBasName) )
     $   Call ErrQuit('rimp2: fitbas bas_create failed', FitBas)
      If ( .NOT. Bas_RTDB_Load( RTDB, Geom, FitBas, FitBasName) )
     $   Call ErrQuit('rimp2: fitbas bas_rtdb_load failed', FitBas)
      If ( .NOT. Int_Normalize( FitBas, RTDB) )
     $   Call ErrQuit('rimp2: fitbas int_normalize failed', FitBas)
C
C     Get other controls -- threshold for neglect of integrals, etc.
C
      Call RIMP2_Get_Param( RTDB, Thresh, MxTBatch, Mem_Opt, MxEBatch,
     $   Approx, MinEVal)
C
C     Figure out file names.  Munging to get the full name is done in
C     open_3ci, not here!
C
      if (.not. RTDB_CGet(rtdb, 'mp2:mo 3-center integral file',
     $     1, FNInt)) FNInt = Def_Int
C
C     Get some info about the basis sets for later use
C
      If ( .NOT. Bas_NumBF( AOBas, NAOBF) )
     $   Call ErrQuit('rimp2: aobas bas_numbf failed', AOBas)
      If ( .NOT. Bas_NumBF( FitBas, NFit) )
     $   Call ErrQuit('rimp2: fitbas bas_numbf failed', FitBas)
      If ( .NOT. Bas_Name( AOBas, AOName, AONameT) )
     $   Call ErrQuit('rimp2: aobas bas_name failed', AOBas)
      If ( .NOT. Bas_Name( FitBas, FitName, FitNameT) )
     $   Call ErrQuit('rimp2: fitbas bas_name failed', FitBas)
C
      Call Corr_Mk_Ref(RTDB, My_Ctx, Geom, AOBas, .false.,
     . RefNAO, TopSpin, NCorBF,G_AOXF, MinChunk, MinChunk, MaxSpin, 
     .     H_Occ, H_EVal, NFrzOcc, NFrzVir, .False.)
C
      If (RefNAO .ne. NAOBF) Call ErrQuit(
     $     'RIMP2: Incompatible reference AO basis size', RefNAO-NAOBF)
      MxNCorBF = NAOBF
      If ( .NOT. MA_Get_Index(H_Occ, I_Occ) ) Call ErrQuit(
     $     'RIMP2: Failed to get index for occ vector', H_Occ)
      If ( .NOT. MA_Get_Index(H_EVal, I_EVal) ) Call ErrQuit(
     $     'RIMP2: Failed to get index for eval vector', H_EVal)
C
C     Determine ranges to transform.
C
      Do I = 1, TopSpin
         Call Occ_Ranges( NCorBF(I), Dbl_MB(I_Occ + MxNCorBF * (I-1) ),
     $      ILo(I), IHi(I), ALo(I), AHi(I) )
      EndDo
C
      If ( .NOT. MA_Free_Heap( H_Occ) ) Call ErrQuit(
     $   'RIMP2: Unable to free local array H_Occ', H_Occ)
C
C     **********************************************
C     * Print out calculation information for user *
C     **********************************************
C
      If ( Me .eq. 0 .AND. PrInfo) then
C
         Write (LuOut,'(1X, A)' ) 'RI-MP2 Calculation Information:'
         Write (LuOut, 8520) 'Nr. of AO basis functions:', NAOBF
         Write (LuOut, 8520) 'Nr. of fitting basis functions:', NFit
         Write (LuOut, 8520) 'Nr. of MO spin cases:', TopSpin
         Write (LuOut, 8520) 'Nr. of MOs to be correlated:',
     $      (NCorBF(i), i = 1, TopSpin)
         Write (LuOut, 8520) 'Nr. of occupied MOs correlated:',
     $      (IHi(i)-ILo(i)+1, i = 1, TopSpin)
         Write (LuOut, 8520) 'Nr. of virtual MOs correlated:',
     $      (AHi(i)-ALo(i)+1, i = 1, TopSpin)
         Write (LuOut, *)
C
         Write (LuOut, 8500) 'RI integral approximation:',
     $      Approx(:LNBlnk(Approx))
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
         Write (LuOut, 8500) 'Base file name for MO integrals',
     $      FNInt(:LNBlnk(FNInt))
#ifdef RIMP2_USE_DA
         Write (LuOut, 8500) 'Type of storage for MO integrals',
     $      'disk arrays'
#else
         Write (LuOut, 8500) 'Type of storage for MO integrals',
     $      'Fortran direct access files'
#endif
         Write (LuOut, *)
C
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom))
     $         call errquit('RIMP2: geom_print', 0)
            write(LuOut,*)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(AOBas)) 
     $           call errquit('RIMP2: bas_print', 0)
            write(LuOut,*)
            if (.not. bas_print(FitBas)) 
     $           call errquit('RIMP2: bas_print', 0)
            write(LuOut,*)
         endif
         if (util_print('geombas', print_debug)) then
            if (.not. gbs_map_print(AOBas)) 
     $           call errquit('RIMP2:gmp', 0)
            write(LuOut,*)
            if (.not. gbs_map_print(FitBas)) 
     $           call errquit('RIMP2:gmp', 0)
            write(LuOut,*)
         endif
         call util_flush(LuOut)
      EndIf
C
 8500 Format(1X, A35, 1X, '''', A, '''')
 8505 Format(1X, A35, 2X, A)
 8510 Format(1X, A35, 4(1X, 1PE10.3), (/36X, 4(1X, I10)) )
 8520 Format(1X, A35, 4(1X, I10), (/36X, 4(1X,I10)) )
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
C     *************************************
C     * Integral transformation (to disk) *
C     *************************************
C
#ifdef RIMP2_USE_DA
C
C     Compute largest and total DA sizes (in bytes)
C
      MaxDASiz = 0
      TotDASiz = 0
      Do I = 1, TopSpin
         Size =  (IHi(i)-ILo(i)+1) * (AHi(i)-ALo(i)+1) * NFit
         SizeB = MA_SizeOf(MT_Dbl, Size, MT_Byte)
         MaxDASiz = Max( MaxDASiz, SizeB)
         TotDASiz = TotDASiz + SizeB
      EndDo
C
C     Initialize DAs.  We know everything but how much memory to give it.
C
      If (DRA_Init(TopSpin, MaxDASiz, TotDASiz, -1) .ne. 0)
     $   Call ErrQuit('RIMP2: DRA_Init failed', 0)
#endif
C
C
C     Initialize integrals
C
      Bases(1) = AOBas
      Bases(2) = FitBas
C
      Call Int_Init(RTDB, 2, Bases)
C
C
      If ( .NOT. RIMP2_XF_Driver( RTDB, Geom, AOBas, FitBas,
     $   TopSpin, D_Int, FNInt, ALo, AHi, ILo, IHi, MxTBatch, Thresh,
     $   MinEVal, Approx, G_AOXF) ) Call ErrQuit(
     $   'rimp2: rimp2_xf_driver failed', 911)
C
C     Finished with the AO transformation
C
      Status = .TRUE.
      Do I = TopSpin, 1, -1
         Status = Status .AND.  GA_Destroy( G_AOXF(I))
      EndDo
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2: Unable to destroy AO transformation', 0)
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
      EndDo
C
C     Create matrix to hold spin case breakdown of E(2)
C
      If ( .NOT. MA_Push_Get( MT_Dbl, TopSpin*TopSpin,
     $   'E(2) spin case breakdown', H_Contrib, Contrib_Idx) ) Call
     $   ErrQuit('rimp2: can''t allocate spin contribution matrix', 911)
C
C     **********************
C     * Energy calculation *
C     **********************
C
      If ( .NOT. RIMP2_V_E2( D_Int, FNInt, TopSpin, NCorBF, IHi, AHi,
     $   NFit, FitBas, Dbl_MB(I_Eval), MxNCorBF, Dbl_MB(Contrib_Idx),
     $   CorEner, Mem_Opt, MxEBatch(1), MxEBatch(2)) ) Call ErrQuit(
     $   'rimp2: rimp2_v_e2 failed', 911)
C
C     Also get the SCF energy so we can report the total energy
C
      Status = RTDB_Get(RTDB, 'scf:energy', MT_Dbl, 1, RefEner)
      If ( .NOT. Status ) then
         RefEner = FP_Zero
         If ( Me .eq. 0) Write ( LuOut, *)
     $      'RIMP2: Warning! Unable to get scf:energy from RTDB.'
      EndIf
C
      TotEner = RefEner + CorEner
C
      If ( Me .eq. 0) then
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
      Status = RTDB_Put( RTDB, 'rimp2:correlation energy', MT_Dbl, 1,
     $   CorEner)
      Status = Status .AND. RTDB_Put( RTDB, 'rimp2:energy', MT_Dbl, 1,
     $     TotEner)
      If ( .NOT. Status ) Call ErrQuit(
     $   'rimp2: unable to store energies in rtdb', 911)
C
#ifdef RIMP2_USE_DA
C
      If ( DRA_Terminate() .ne. 0) Call ErrQuit(
     $   'RIMP2: DRA_Terminate failed', 0)
#endif
C
C     ***********
C     * Cleanup *
C     ***********
C
      Status = .TRUE.
      Status = Status .AND. MA_Pop_Stack( H_Contrib)
      Status = Status .AND.  MA_Free_Heap(H_Eval)
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2: Unable to destroy local arrays', 0)
C
      Status = Status .AND. Geom_Destroy( Geom)
      Status = Status .AND. Bas_Destroy( FitBas)
      Status = Status .AND. Bas_Destroy( AOBas)
      If (.NOT. Status) Call ErrQuit(
     $   'RIMP2: Failed to destroy geometry and  basis sets', 0)
C
      Call Util_Print_Pop
      call ecce_print_module_exit('mp2', 'ok')
C
      Call PStat_Off( PS_Overall)
      Call RIMP2_PStat_Print
      Call RIMP2_PStat_Term
C
      RIMP2 = .TRUE.
      Return
      End
