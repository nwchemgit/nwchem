      Logical Function RIMP2( RTDB )
C$Id: rimp2.F,v 1.38 1998-01-31 20:14:20 bernhold Exp $
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
      Logical RIMP2_XF_Driver, RIMP2_V_E2
      External RIMP2_XF_Driver, RIMP2_V_E2
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
      Integer NAOBF, TopSpin, RefSpin, NFit
      Integer ALo(MaxSpin), AHi(MaxSpin), ILo(MaxSpin), IHi(MaxSpin)
      Integer RefALo(MaxSpin), RefAHi(MaxSpin), RefILo(MaxSpin),
     $   RefIHi(MaxSpin)
      Integer NMO(MaxSpin), RefNMO(MaxSpin),I, Me, SpinMap(MaxSpin)
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
C     Handles and indices for the full SCF reference basis
C     These are used only briefly to hold the full info
C     while the map to the correlated basis is constructed and
C     everything is gathered into the correlated basis stuff
C     above.
C
      Integer H_Full_Occ, H_Full_EVal
      Integer I_Full_Occ, I_Full_EVal
      Integer G_Full_AOXF(MaxSpin)
      Integer MxNMO
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
      If (.NOT. RTDB_CGet(RTDB, 'title', 1, JobTitle)) JobTitle = ' '
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
      Call RIMP2_Get_Param( RTDB, Thresh, MxTBatch, Mem_Opt, MxEBatch,
     $   Approx, MinEVal, TopSpin, SpinMap, MaxSpin)
C
C     Figure out file names.  Munging to get the full name is done in
C     open_3ci, not here!
C
      if (.not. RTDB_CGet(rtdb, 'mp2:mo 3-center integral file',
     $     1, FNInt)) FNInt = Def_Int
C
C     Get some info about the basis sets for later use
C
      Status = Status .AND. Bas_NumBF( AOBas, NAOBF)
      Status = Status .AND. Bas_NumBF( FitBas, NFit)
      Status = Status .AND. Bas_Name(AOBas, AOName, AONameT)
      Status = Status .AND. Bas_Name(FitBas, FitName, FitNameT)
C
C     Figure out where to get the reference from
C
      if (.not. RTDB_CGet(rtdb, 'mp2:input vectors',
     $     1, FNVec)) then
         if (.not. RTDB_CGet(rtdb, 'scf:output vectors',
     $      1, FNVec)) call util_file_name(Def_movec,
     $        .false.,.false., FNVec)
      EndIf
C
C     Check MOVecs file for number of spin cases, and the number of MOs
C     contained.
C
      If ( .NOT.  MOVecs_Read_Header(FNVec, Title, SCFBasNam, 
     $     scftype, NAOBF,
     $   RefSpin, RefNMO, MaxSpin) )
     $   Call ErrQuit('RIMP2: Unable to open SCF eigenvectors', 0)
C
C     MOVecs_Read/MOVecs_Gather (for historical reasons) require the
C     occupation & eigenvalue vectors to be NAOBF rather than the actual
C     number of MOs.  Since MxNMO sizes these quantities in this part of
C     the code, we must insure it is large enough.
C
      MxNMO = NAOBF
c$$$      MxNMO = 0
      Do I = 1, RefSpin
         MxNMO = Max( MxNMO, RefNMO(I) )
      EndDo
C
C     Compare reference information with requested spin mapping
C
      If (TopSpin .eq. -1) then
         TopSpin = RefSpin
         Do I = 1, TopSpin
            SpinMap(i) = i
         EndDo
      Else
         Do I = 1, TopSpin
            If ( SpinMap(i) .lt. 1 .OR. SpinMap(i) .gt. RefSpin)
     $         Call ErrQuit('RIMP2: Bad spin map entry', SpinMap(i))
         EndDo
      EndIf
C
C     Produce a correct NMO array and decide if this really is an
C     abnormal reference so we can remind the user below.
C
      FunnyRef = TopSpin .ne. RefSpin
      Do I = 1, TopSpin
         NMO(i) = RefNMO( SpinMap(i) )
         FunnyRef = FunnyRef .OR. (SpinMap(i) .ne. i)
      EndDo
C
C     *************************************
C     * Get info about full SCF reference *
C     *************************************
C     All this effort is to get the occupation vevctors and reduce
C     them to ranges so we can tell the user as a purely informational
C     point.  I think this is useful to helping verify that the
C     calculation is setup correctly, so I'm willing to do the work.
C
C     It might be convenient if we had a way to read only the
C     occupation vector or just the 1-d data from the movecs file.
C
      Status = .TRUE.
      Status = Status .AND. MA_Alloc_Get(MT_Dbl, MxNMO * TopSpin,
     $   'SCF occupation', H_Full_Occ, I_Full_Occ)
      Status = Status .AND. MA_Alloc_Get(MT_Dbl, MxNMO * TopSpin,
     $   'SCF eigenvalues', H_Full_Eval, I_Full_Eval)
C
      Do I = 1, TopSpin
         Status = Status .AND.  GA_Create(MT_Dbl, NAOBF, NMO(I),
     $      'AO transformation', MinChunk, MinChunk, G_Full_AOXF(I))
      EndDo
C
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2: Insufficient space for SCF reference info', 0)
C
C     Read in the MOVECS data
C
      Status = .TRUE.
      Do I = 1, TopSpin
         Status = Status .AND. MOVecs_Read( FNVec, SpinMap(I),
     $      Dbl_MB( I_Full_Occ + MxNMO * (I-1) ),
     $      Dbl_MB(I_Full_EVal + MxNMO * (I-1) ), G_Full_AOXF(I) )
      EndDo
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2: Unable to read SCF reference info', 0)
C
C     Determine size of occupied and virtual spaces in full SCF ref.
C
      Do I = 1, TopSpin
         Call Occ_Ranges( NMO(I), Dbl_MB(I_Full_Occ + MxNMO * (I-1) ),
     $      RefILo(I), RefIHi(I), RefALo(I), RefAHi(I) )
      EndDo
C
      Status = .TRUE.
      Status = Status .AND. MA_Free_Heap( H_Full_Occ  )
      Status = Status .AND. MA_Free_Heap( H_Full_Eval )
      Do I = 1, TopSpin
         Status = Status .AND. GA_Destroy(G_Full_AOXF(I) )
      EndDo
      If (.NOT. Status) Call ErrQuit(
     $   'RIMP2: Unable destroy local arrays', 1)
C
C     **********************************
C     * Create wave reference function *
C     **********************************
C     Create the gather vector which maps the full SCF reference into
C     the "correlated basis" in which some orbitals may be frozen.
C
      If ( .NOT. MA_Push_Get(MT_Int, MxNMO * TopSpin,
     $      'Correlation->SCF basis map', H_CorMap, I_CorMap) )
     $   Call ErrQuit('RIMP2: ma failed', MxNMO*TopSpin)
C
      Call Mk_Correl_Map(RTDB, Geom, TopSpin, NMO, NCorBF, MxNMO,
     $   Int_MB(I_CorMap), FullRef )
C
C     MOVecs_Read/MOVecs_Gather (for historical reasons) require the
C     occupation & eigenvalue vectors to be NAOBF rather than the actual
C     number of MOs.  Since MxNCorBF sizes these quantities in this part of
C     the code, we must insure it is large enough.
C
      MxNCorBF = NAOBF
c$$$      MxNCorBF = 0
      Do I = 1, TopSpin
         MxNCorBF = Max( MxNCorBF, NCorBF(I) )
      EndDo
C
      Status = .TRUE.
      Status = Status .AND. MA_Alloc_Get(MT_Dbl, MxNCorBF * TopSpin,
     $   'SCF occupation', H_Occ, I_Occ)
      Status = Status .AND. MA_Push_Get(MT_Dbl, MxNCorBF * TopSpin,
     $   'SCF eigenvalues', H_Eval, I_Eval)
C
      Do I = 1, TopSpin
         Status = Status .AND.  GA_Create(MT_Dbl, NAOBF, NCorBF(I),
     $      'AO transformation', MinChunk, MinChunk, G_AOXF(I))
      EndDo
C
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2: Insufficient space for SCF reference info', 1)
C
C     Read in the correlated basis occupation and eigenstuff.
C
      Do I = 1, TopSpin
         If ( .NOT. MOVecs_Gather( FNVec, SpinMap(I), NCorBF(I),
     $      Int_MB(I_CorMap + MxNMO * (I-1) ),
     $      Dbl_MB(I_Occ + MxNCorBF * (I-1) ),
     $      Dbl_MB(I_Eval + MxNCorBF * (I-1) ), G_AOXF(I) ) )
     $      Call ErrQuit('RIMP2: Unable to read reference', I)
      EndDo
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
         Call Util_Print_Centered(LuOut, 'NWChem RI-MP2 Module', 40,
     $   .TRUE.)
         If ( JobTitle .ne. ' ') Call Util_Print_Centered( LuOut,
     $      JobTitle, 40, .FALSE.)
         Write (LuOut, '(//)')
C
         Write (LuOut, 8520) 'Nr. of AO basis functions:', NAOBF
         Write (LuOut, 8520) 'Nr. of fitting basis functions:', NFit
         If ( FunnyRef ) then
            Write (LuOut, 8520) 'Nr. of spin cases in ref.:', RefSpin
            Write (LuOut, 8520) 'Nr. of spin cases in calc.:', TopSpin
            Write (LuOut, 8520) 'Mapping of ref. spin cases to calc.:',
     $         (SpinMap(i), i = 1, TopSpin)
         Else
            Write (LuOut, 8520) 'Nr. of MO spin cases:', TopSpin
         EndIf
         If ( .NOT. FullRef ) then
            Write (LuOut, 8520) 'Nr. of MOs in SCF reference:',
     $         (NMO(i), i = 1, TopSpin)
            Write (LuOut, 8520) 'Nr. of occupied MOs in SCF ref.:',
     $         (RefIHi(i)-RefILo(i)+1, i = 1, TopSpin)
            Write (LuOut, 8520) 'Nr. of virtual MOs in SCF ref.:',
     $         (RefAHi(i)-RefALo(i)+1, i = 1, TopSpin)
         EndIf
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
         Write (LuOut, 8500) 'SCF reference file:',
     $      FNVec(:LNBlnk(FNVec))
         Write (LuOut, 8500) 'SCF reference title:',
     $      Title(:LNBlnk(Title))
         Write(LuOut,*)
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
C     Warn users that everything from this point is labeled
C     in a different basis.
C
      PrWarn = Util_Print('warnings', Print_Low)
      PrRef  = Util_Print('reference', Print_High)
C
      If ( .NOT. FullRef .AND. Me .eq. 0) then
         If ( PrRef) then
            Write (LuOut, '(A)')
            Call Banner(LuOut, 'IMPORTANT NOTE', '*',
     $         .TRUE., .TRUE., .TRUE.)
C
            Write (LuOut, 9050)
C
            Call Print_Correl_Map(LuOut, TopSpin, NCorBF,
     $         Int_MB(I_CorMap), MxNCorBF, 'I5',
     $         Dbl_MB(I_Occ), MxNCorBF, 'F3.1',
     $         Dbl_MB(I_Eval), MxNCorBF, 'F12.3', 80)
         ElseIf( PrWarn ) then
            Write (LuOut, '(A)')
            Call Banner(LuOut, 'IMPORTANT NOTE', '*',
     $         .TRUE., .TRUE., .TRUE.)
C
            Write (LuOut, 9051)
         EndIf
         Call Util_Flush(LuOut)
      EndIF
C
 8500 Format(1X, A35, 1X, '''', A, '''')
 8505 Format(1X, A35, 2X, A)
 8510 Format(1X, A35, 4(1X, 1PE10.3), (/36X, 4(1X, I10)) )
 8520 Format(1X, A35, 4(1X, I10), (/36X, 4(1X,I10)) )
C
 9050 Format(1X, 'This correlated calculation is based on an SCF ',
     $       'reference in which'/
     $   1X, 'orbitals have been frozen or rearranged relative to ',
     $        'the result of the'/
     $   1X, 'SCF calculation.  The following mapping relates the ',
     $       '''correlated'/
     $   1X, 'basis'' used from here on out to the original SCF ',
     $       'basis.  From this'/
     $   1X, 'point on, orbitals are labeled in the ''correlated ',
     $       'basis''.'/)
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
      Status = Status .AND.
     $   RIMP2_XF_Driver( RTDB, Geom, AOBas, FitBas,
     $   TopSpin, D_Int, FNInt, ALo, AHi, ILo, IHi, MxTBatch, Thresh,
     $   MinEVal, Approx, G_AOXF)
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
      Status = Status .AND. MA_Push_Get( MT_Dbl, TopSpin*TopSpin,
     $   'E(2) spin case breakdown', H_Contrib, Contrib_Idx)
C
C     **********************
C     * Energy calculation *
C     **********************
C
      Status = Status .AND.
     $   RIMP2_V_E2( D_Int, FNInt, TopSpin, NCorBF, IHi, AHi, NFit,
     $   FitBas, Dbl_MB(I_Eval), MxNCorBF, Dbl_MB(Contrib_Idx),
     $   CorEner, Mem_Opt, MxEBatch(1), MxEBatch(2))
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
      Status = Status .AND. RTDB_Put( RTDB, 'rimp2:correlation energy',
     $     MT_Dbl, 1, CorEner)
      Status = Status .AND. RTDB_Put( RTDB, 'rimp2:energy', MT_Dbl, 1,
     $     TotEner)
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
      Status = Status .AND.  MA_Pop_Stack(H_Eval)
      Status = Status .AND.  MA_Pop_Stack(H_CorMap)
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2: Unable to destroy local arrays', 0)
C
      Status = Status .AND. Geom_Destroy( Geom)
      Status = Status .AND. Bas_Destroy( FitBas)
      Status = Status .AND. Bas_Destroy( AOBas)
      If (.NOT. Status) Call ErrQuit(
     $   'RIMP2: Failed to destroy geometry and  basis sets', 0)
C
      Call PStat_Off( PS_Overall)
      Call RIMP2_PStat_Print
      Call RIMP2_PStat_Term
C
c$$$C
c$$$C     Make sure the user gets everything we've written before going on
c$$$C
c$$$      Call Util_Flush(LuOut)
c$$$      Call Util_Flush(LuErr)
C
      RIMP2 = .TRUE.
      Return
      End
