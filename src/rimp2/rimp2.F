      Logical Function RIMP2( RTDB )
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
      Logical MOVecs_Read_Header, Int_Normalize
      External MOVecs_Read_Header, Int_Normalize
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
      Integer NILimit
c$$$      Parameter (NILimit = 10)
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
      Integer H_Eval, Eval_Idx, H_Contrib, Contrib_Idx
      Double precision Thresh, CorEner, RefEner, TotEner
      Character*(80) Title, SCFBasNam, FNVec, FNInt, String, Approx
      Logical Status
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
      If ( Me .eq. 0 )
     $   Status = Status .AND. Geom_Print(Geom)
C
      Status = Status .AND. Bas_Create( AOBas, AOBasName)
      Status = Status .AND. Bas_RTDB_Load( RTDB, Geom, AOBas, AOBasName)
c$$$      Call Banner( LuOut, 'AO basis prior to normalization', '-',
c$$$     $   .TRUE., .TRUE., .FALSE.)
c$$$      Status = Status .AND. Bas_Print( AOBas)
      Status = Status .AND. Int_Normalize( AOBas)
c$$$      Call Banner( LuOut, 'AO basis after normalization', '-',
c$$$     $   .TRUE., .TRUE., .FALSE.)
c$$$      Status = Status .AND. Bas_Print( AOBas)
c$$$      Status = Status .AND. GBS_Map_Print( AOBas)
C
      Status = Status .AND. Bas_Create( FitBas, FitBasName)
      Status = Status .AND. Bas_RTDB_Load( RTDB, Geom, FitBas,
     $   FitBasName)
c$$$      If ( Me .eq. 0)
c$$$     $   Call Banner( LuOut, 'Fitting basis prior to normalization',
c$$$     $   '-', .TRUE., .TRUE., .FALSE.)
c$$$      Status = Status .AND. Bas_Print( FitBas)
      Status = Status .AND. Int_Normalize( FitBas)
      If ( Me .eq. 0) then
         Call Banner( LuOut, 'Fitting basis after normalization',
     $      '-', .TRUE., .TRUE., .FALSE.)
         Status = Status .AND. Bas_Print(FitBas)
c$$$         Status = Status .AND. GBS_Map_Print( FitBas)
      EndIf
C
C     Initialize integrals
C
      Bases(1) = AOBas
      Bases(2) = FitBas
C
      Call Int_Init( 2, Bases)
C
C     Figure out file names
C
*      if (.not. context_rtdb_match(rtdb, 'mo vector file', movecs))
*     $     movecs = file_prefix_name('movecs')
      if (.not. context_rtdb_match(rtdb, 'input mo vector file',
     $     FNVec)) FNVec = Def_movec
C
      if (.not. context_rtdb_match(rtdb, 'mo 3-center integral file',
     $     FNInt)) FNInt = Def_Int
C
C     Check MOVecs file for number of spin cases
C
      If ( .NOT.  MOVecs_Read_Header(FNVec, Title, SCFBasNam, NAOBF,
     $   TopSpin, NMO, MaxSpin) )
     $   Call ErrQuit('RIMP2: Unable to open SCF eigenvectors', 0)
C
      If ( Me .eq. 0) then
         Call Banner(LuOut, 'SCF Eigenvector information', '-',
     $      .TRUE., .TRUE., .FALSE.)
C
         Write (LuOut, 9000) Title(:LNBlnk(Title)),
     $      SCFBasNam(:LNBlnk(SCFBasNam)), NAOBF, TopSpin,
     $      (NMO(i), i = 1, TopSpin)
      EndIF
C
 9000 Format(1X, 'Source job: ''', A, ''''
     $   /1X,    'Basis:      ''', A, '''',
     $   /1X,    'Nr. of AOs:  ', I5,
     $   /1X,    'Nr. of sets: ', I5,
     $   /1X,    'Nr. of MOs:  ', 10I5)

C
C     Get the basis set sizes we need
C
      Status = Status .AND. Bas_NumBF( AOBas, NAOBF)
      Status = Status .AND. Bas_NumBF( FitBas, NFit)
C
C     Create space for eigenvalues
C
      Status = Status .AND. MA_Push_Get( MT_Dbl, NAOBF*TopSpin,
     $   'SCF eigenvalues', H_EVal, EVal_Idx)
C
C     Get other controls -- threshold for neglect of integrals, etc.
C
      Call RIMP2_Get_Param( RTDB, Thresh, NILimit, Approx)
C
C     *************************************
C     * Integral transformation (to disk) *
C     *************************************
C
c$$$      Write (LuOut, *) 'About to call transformation'
      Status = Status .AND.
     $   RIMP2_XF_Driver( Geom, AOBas, FitBas, Dbl_MB(EVal_Idx), NAOBF,
     $   TopSpin, FNVec, FNInt, ALo, AHi, ILo, IHi, NILimit, Thresh,
     $   Approx)
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
c$$$      Write (LuOut, *) 'About to compute energy'
      Status = Status .AND.
     $   RIMP2_V_E2( FNInt, TopSpin, NAOBF, IHi, AHi, NFit, FitBas,
     $   Dbl_MB(EVal_Idx), Dbl_MB(Contrib_Idx), CorEner, NILimit,
     $   NILimit)
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
c$$$      Write (LuOut, *) 'Cleaning up'
      Status = Status .AND. MA_Pop_Stack( H_Contrib)
      Status = Status .AND. MA_Pop_Stack( H_EVal)
C
      Call Int_Terminate
C
      Status = Status .AND. Geom_Destroy( Geom)
      Status = Status .AND. Bas_Destroy( FitBas)
      Status = Status .AND. Bas_Destroy( AOBas)
C
      Call PStat_Off( PS_Overall)
      Call RIMP2_PStat_Print
      Call RIMP2_PStat_Term
C
      Status = Status .AND. Context_Pop(My_Ctx)
C
      RIMP2 = .TRUE.
      Return
      End
