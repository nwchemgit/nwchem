      Logical Function RIMP2( RTDB )
C$Id: rimp2.F,v 1.15 1995-02-02 23:21:17 d3g681 Exp $
      Implicit NONE
      Integer RTDB
C
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.h"
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "pstat.fh"
#include "rimp2_timers.fh"
#include "numerical_constants.fh"
C
      Logical MOVecs_Read_Header, MOVecs_Read, Int_Normalize
      External MOVecs_Read_Header, MOVecs_Read, Int_Normalize
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
      Integer MxTBatch, MxEBatch
C
      Character*(*) Def_Int, Def_MOVec
      Parameter (Def_Int = 'mo3cint', Def_MOVec = 'movecs')
      Character*(*) Geom_Name, AOBasName, FitBasName
      Parameter (Geom_Name = 'geometry', AOBasName = 'mo basis',
     $   FitBasName = 'ri-mp2 basis')
C
      Character*(*) My_Ctx
      Parameter (My_Ctx = 'mp2')
C
      Integer AOBas, FitBas, Bases(2), Geom
      Integer NAOBF, TopSpin, NFit
      Integer ALo(MaxSpin), AHi(MaxSpin), ILo(MaxSpin), IHi(MaxSpin)
      Integer NMO(MaxSpin), I, Me
      Integer H_Contrib, Contrib_Idx
      Double precision Thresh, CorEner, RefEner, TotEner
      Character*(80) Title, SCFBasNam, FNVec, FNInt, String, Approx
      Character*(255) AOName, AONameT, FitName, FitNameT
      Logical Status
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
      Logical FullRef
C
      Status = .TRUE.
      Status = Status .AND. Context_Push(My_Ctx)
C
      Call RIMP2_PStat_Init
      Call PStat_On( PS_Overall)
C
      Me = GA_NodeID()
C
C     Print a headline
C
      If ( Me .eq. 0)
     $   Call Util_Print_Centered(LuOut, 'NWCHEM RI-MP2 Energy', 40,
     $   .TRUE.)
C
C     Load the geometry and basis sets
C
      Status = Status .AND. Geom_Create( Geom, Geom_Name)
      Status = Status .AND. Geom_RTDB_Load( RTDB, Geom, Geom_Name)
c$$$      If ( Me .eq. 0 )
c$$$     $   Status = Status .AND. Geom_Print(Geom)
C
      Status = Status .AND. Bas_Create( AOBas, AOBasName)
      Status = Status .AND. Bas_RTDB_Load( RTDB, Geom, AOBas, AOBasName)
      Status = Status .AND. Int_Normalize( AOBas)
C
      Status = Status .AND. Bas_Create( FitBas, FitBasName)
      Status = Status .AND. Bas_RTDB_Load( RTDB, Geom, FitBas,
     $   FitBasName)
      Status = Status .AND. Int_Normalize( FitBas)
c$$$      If ( Me .eq. 0) then
c$$$         Call Banner( LuOut, 'Fitting basis after normalization',
c$$$     $      '-', .TRUE., .TRUE., .FALSE.)
c$$$         Status = Status .AND. Bas_Print(FitBas)
c$$$c$$$         Status = Status .AND. GBS_Map_Print( FitBas)
c$$$      EndIf
C
C     Get other controls -- threshold for neglect of integrals, etc.
C
      Call RIMP2_Get_Param( RTDB, Thresh, MxTBatch, MxEBatch, Approx)
C
C     Figure out file names
C
      if (.not. RTDB_Get(rtdb, 'mp2:mo 3-center integral file',
     $     MT_Int, 1, FNInt)) FNInt = Def_Int
C
C     Get some info about the basis sets for later use
C
      Status = Status .AND. Bas_NumBF( AOBas, NAOBF)
      Status = Status .AND. Bas_NumBF( FitBas, NFit)
      Status = Status .AND. Bas_Name(AOBas, AOName, AONameT)
      Status = Status .AND. Bas_Name(FitBas, FitName, FitNameT)
C
C     ****************************
C     * Get/create SCF reference *
C     ****************************
C     If the user hasn't specified a set of vectors, then default
C     to the SCF output vectors.  
C
      if (.not. RTDB_CGet(rtdb, 'mp2:input vectors',
     $     1, FNVec)) then
         if (.not. RTDB_CGet(rtdb, 'scf:output vectors',
     $      1, FNVec)) FNVec = Def_movec
      EndIf
C
C     Check MOVecs file for number of spin cases
C
      If ( .NOT.  MOVecs_Read_Header(FNVec, Title, SCFBasNam, NAOBF,
     $   TopSpin, NMO, MaxSpin) )
     $   Call ErrQuit('RIMP2: Unable to open SCF eigenvectors', 0)
C
C
C     Create space for all of the data from the MOVECS file
C     and for the mapping functions
C
      MxNMO = 0
      Do I = 1, TopSpin
         MxNMO = Max( MxNMO, NMO(I) )
      EndDo
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
      Status = Status .AND. MA_Push_Get(MT_Int, MxNMO * TopSpin,
     $      'Correlation->SCF basis map', H_CorMap, I_CorMap)
C
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2: Insufficient space for SCF reference info', 0)
C
C     Read in the MOVECS data
C
      Status = .TRUE.
      Do I = 1, TopSpin
         Status = Status .AND. MOVecs_Read( FNVec, I,
     $      Dbl_MB( I_Full_Occ + MxNMO * (I-1) ),
     $      Dbl_MB(I_Full_EVal + MxNMO * (I-1) ), G_Full_AOXF(I) )
      EndDo
      If ( .NOT. Status) Call ErrQuit(
     $   'RIMP2: Unable to read SCF reference info', 0)
C
C     Create the mapping functions
C
      Call Mk_Correl_Map(RTDB, MxNMO, TopSpin, NMO, Dbl_MB(I_Full_Occ),
     $   Dbl_MB(I_Full_Eval), NCorBF, Int_MB(I_CorMap), FullRef )
C
C     Create space for all of the data we'll use in the calculation
C
      MxNCorBF = 0
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
C     Form the new occupation, eigenvalues, and orbitals for
C     the correlated calculation
C
      Do I = 1, TopSpin
         Call Drop_Orbitals(
     $      Dbl_MB(I_Full_Occ + MxNMO * (I-1) ),
     $      Dbl_MB(I_Full_Eval + MxNMO * (I-1) ), G_Full_AOXF(I),
     $      NCorBF(I), Int_MB(I_CorMap + MxNMO * (I-1) ),
     $      Dbl_MB(I_Occ + MxNCorBF * (I-1) ),
     $      Dbl_MB(I_Eval + MxNCorBF * (I-1) ), G_AOXF(I) )
      EndDo
C
C     Now we can get rid of the full reference
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
C     **********************************************
C     * Print out calculation information for user *
C     **********************************************
C
      If ( Me .eq. 0) then
         Write (LuOut, 8500) 'RI integral approximation:',
     $      Approx(:LNBlnk(Approx))
         Write (LuOut, 8510) 'Threshold for Schwarz screening:', Thresh
         Write (LuOut, 8520) 'Max. batch size for transformation:',
     $      MxTBatch
         Write (LuOut, 8520) 'Max. batch size for energy:',
     $      MxEBatch
         Write (LuOut, 8500) 'Base file name for MO integrals',
     $      FNInt(:LNBlnk(FNInt))
         Write (LuOut, 8500) 'AO basis set name:',
     $      AONameT(:LNBlnk(AONameT))
         Write (LuOut, 8520) 'Nr. of AO functions:', NAOBF
         Write (LuOut, 8500) 'Fitting basis set name:',
     $      FitName(:LNBlnk(FitName))
         Write (LuOut, 8520) 'Nr. of fitting functions:', NFit
C
         Write (LuOut, 8500) 'SCF reference from file:',
     $      FNVec(:LNBlnk(FNVec))
         Write (LuOut, 8500) 'SCF reference title:',
     $      Title(:LNBlnk(Title))
         Write (LuOut, 8500) 'SCF reference basis:',
     $      SCFBasNam(:LNBlnk(SCFBasNam))
         Write (LuOut, 8520) 'Nr. of AOs in SCF reference:', NAOBF
         Write (LuOut, 8520) 'Nr. of sets of orbitals:', TopSpin
         Write (LuOut, 8520) 'Nr. of MOs in full reference:',
     $      (NMO(i), i = 1, TopSpin)
         Write (LuOut, 8520) 'Nr. of MOs to be correlated:',
     $      (NCorBF(i), i = 1, TopSpin)
C
C        Warn users that everything from this point is labeled
C        in a different basis.
C
         If ( .NOT. FullRef ) then
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
         EndIf
      EndIF
      Call Util_Flush(LuOut)
C
 8500 Format(1X, A35, 1X, '''', A, '''')
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
C
C     *************************************
C     * Integral transformation (to disk) *
C     *************************************
C
C     Determine ranges to transform
C
      Do I = 1, TopSpin
         Call Occ_Ranges( NCorBF(I), Dbl_MB(I_Occ + MxNCorBF * (I-1) ),
     $      ILo(I), IHi(I), ALo(I), AHi(I) )
      EndDo
C
C     Done with occupation vector
C
      If ( .NOT. MA_Free_Heap( H_Occ) ) Call ErrQuit(
     $   'RIMP2: Unable to free local array H_Occ', H_Occ)
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
     $   TopSpin, FNInt, ALo, AHi, ILo, IHi, MxTBatch, Thresh,
     $   Approx, G_AOXF)
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
     $   RIMP2_V_E2( FNInt, TopSpin, NCorBF, IHi, AHi, NFit, FitBas,
     $   Dbl_MB(I_Eval), MxNCorBF, Dbl_MB(Contrib_Idx), CorEner,
     $   MxEBatch, MxEBatch)
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
      Status = Status .AND. Context_Prefix('correlation energy', String)
      Status = Status .AND. RTDB_Put( RTDB, String, MT_Dbl, 1, CorEner)
      Status = Status .AND. Context_Prefix('total energy', String)
      Status = Status .AND. RTDB_Put( RTDB, String, MT_Dbl, 1, TotEner)
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
      Status = Status .AND. Context_Pop(My_Ctx)
      If (.NOT. Status) Call ErrQuit(
     $   'RIMP2: Failed to pop context', 0)
C
      RIMP2 = .TRUE.
      Return
      End
