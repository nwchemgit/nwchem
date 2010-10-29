C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_E_IJ -- Evaluate RI-MP2 energy contributions for IJ batch
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_E_IJ(TopSpin, SameSpin, ILo, IHi, NVirI,
     $   JLo, JHi, NVirJ, OccEigI, OccEigJ, VirEigI, VirEigJ,
     $   G_I, KeepI, G_J, KeepJ, G_Int, KeepInt, G_T2, KeepT2, G_PairE,
     $   SpinLabel, TopN, TopN_Ind, TopN_Val, ReqI, ReqJ)
      Implicit NONE
#include "errquit.fh"
C
      Integer TopSpin  ![in] RHF=1 UHF=2
      Logical SameSpin ![in] I and J have same spin
      Integer ILo      ![in] Lower limit of I in current batch
      Integer IHi      ![in] Upper limit of I in current batch
      Integer NVirI    ![in] Number of virtuals in I spin case
      Integer JLo      ![in] Lower limit of J in current batch
      Integer JHi      ![in] Upper limit of J in current batch
      Integer NVirJ    ![in] Number of virtuals in J spin case
      Double Precision OccEigI(*)      ![in] Occ eigenvalues for I's spin
      Double Precision OccEigJ(*)      ![in] Occ eigenvalues for J's spin
      Double Precision VirEigI(NVirI)  ![in] Vir eigenvalues for I's spin
      Double Precision VirEigJ(NVirJ)  ![in] Vir eigenvalues for I's spin
      Integer G_I      ![in]  3c2e integrals for I
      Logical KeepI    ![in]  Keep G_I around for caller
      Integer G_J      ![in]  3c2e integrals for J
      Logical KeepJ    ![in]  Keep G_J around for caller
      Integer G_Int    ![out] (approx) antisym 4c2e integral array (out)
      Logical KeepInt  ![in]  Keep G_Int around for caller
      Integer G_T2     ![out] T2 amplitude array (out)
      Logical KeepT2   ![in]  Keep G_T2 around for caller
      Integer G_PairE  ![in]  Pair energies for this spin case (inout)
      Character*(*) SpinLabel ![in] text label for current spin case
      Integer TopN
      Integer TopN_Ind(TopN, 4)
      Double precision TopN_Val(TopN)
      Integer ReqI ![in] DRA request id for I integral block
      Integer ReqJ ![in] DRA request id for J integral block
C
C DESCRIPTION
C     Given batches of 3c2e integrals, (ILo:IHi, JLo:JHi), produce the
c     RI-MP2 pair energies and (optionally) approximate antisymmetrized
c     4c2e integrals and T2 amplitudes.
C
C     Input 3c2e integrals are organized as (ia) x t' matrix with the
C     combined index ia having all a for some range ilo:ihi (a running
C     fastest).
C
C     4-index integrals & amplitudes are organized as blocks of all AB
c     for a particular IJ. A batch will contain blocks for multiple Is
c     (row dimension) and Js (column dimension) laid out in rectangular
c     fashion.  If SameSpin is .TRUE., then permutational symmetry of
c     I <-> J should be used, so that only the I > J triangle of the
c     batch will actually be calculated.
C
C EXTERNAL ROUTINES
      Logical  RIMP2_Asym_I, RIMP2_Mk_T
      External RIMP2_Asym_I, RIMP2_Mk_T
C
C INCLUDED FILES
#include "global.fh"
#include "mafdecls.fh"
#include "dra.fh"
#include "util.fh"
#include "rimp2g_timers.fh"
#include "stdio.fh"
#include "inp.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer Chunk               ! chunk size for GA distr. of G_Int, G_T2
      Parameter (Chunk = -1)      ! Chosen to give even distr.
      Logical IJTri               ! G_Int, G_T2 contain only unique IJ
      Parameter (IJTri = .FALSE.) ! Keep them rectangular
C
C LOCAL VARIABLES
      Integer NI, NJ, NFit, Type, Junk
      Logical PrPrgRpt
      Character*255 String1
C
      NI = IHi - ILo + 1
      NJ = JHi - JLo + 1
C
c$$$      PrPrgRpt = Util_Print('progress reports', Print_Default)
C     
C     Create GA for (approximate) 4c2e integrals.
C     
      String1 = 'Approx 4e2c ints ' // SpinLabel
      If ( .NOT. GA_Create( MT_Dbl, NI*NVirI, NJ*NVirJ,
     $   String1(:Inp_StrLen(String1)), Chunk, Chunk, G_Int) )
     $   Call ErrQuit( 'RIMP2_E_IJ: can''t allocate integral GA',
     $   NI*NVirI*NJ*NVirJ , GA_ERR)
c$$$C
c$$$C     At this point, we must insure that the integrals are actually
c$$$C     finished loading.
c$$$C
c$$$      Call PStat_On( PS_Rd_3ci)
c$$$      If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
c$$$     $     'RIMP2_IJ_Int: Error from DRA_Wait on I integrals', 1)
c$$$      If ( DRA_Wait(ReqJ) .ne. 0) Call ErrQuit(
c$$$     $     'RIMP2_IJ_Int: Error from DRA_Wait on J integrals', 2)
c$$$      Call PStat_Off( PS_Rd_3ci)
C     
C     Create approximate 4c2e integrals from these three-center things:
C     (ia|jb) = (ia|s')(s'|jb)
C     where the |s') = V^{-1/2} |s) has been done as part of the 
C     transformation.
C     
      Call GA_Sync
      Call PStat_On( PS_form_4ci)
      Call GA_Inquire(G_I, Type, Junk, NFit)
      Call GA_DGEMM('N', 'T',NI*NVirI, NJ* NVirJ, NFit,
     $   1.0d0, G_I, G_j, 0.0d0, G_Int)
      Call PStat_Off( PS_form_4ci)
c$$$      If ( PrPrgRpt ) Call Progress_Report(LuOut,
c$$$     $   '       Approximate 4c2e integrals done')
      If ( Util_Print('4c ints', Print_Debug) )
     $   Call GA_Print( G_Int )
C
C     Don't need I and J any more
C
      If ( .NOT. KeepJ) then
         If ( .NOT. GA_Destroy( G_j) ) Call ErrQuit(
     $      'RIMP2_E_IJ: can''t destroy integral GA J', G_j, GA_ERR)
      EndIf
      If ( .NOT. KeepI) then
         If ( .NOT. GA_Destroy( G_I) ) Call ErrQuit(
     $      'RIMP2_E_IJ: can''t destroy integral GA I', G_I, GA_ERR)
      EndIf
C     
C     Antisymmetrize integrals
C     
      Call PStat_On( PS_form_ab)
      If ( .NOT. RIMP2_Asym_I(TopSpin, SameSpin, ILo, IHi, NVirI, 
     $   JLo, JHi, NVirJ, .FALSE., IJTri, G_Int, SpinLabel) )
     $   Call ErrQuit('RIMP2G_V_E2: can''t antisymmetrize integrals', 0,
     &       INT_ERR)
      Call PStat_Off( PS_form_ab)
c$$$      If ( PrPrgRpt ) Call Progress_Report(LuOut,
c$$$     $   '       Antisymmetrized integrals done')
C     
C     Make amplitudes
C     
      Call PStat_On( PS_denom)
      If (.NOT. RIMP2_Mk_T(ILo, IHi, NVirI, JLo, JHI, NVirJ, .FALSE.,
     $   IJTri, G_Int, G_T2, OccEigI, OccEigJ, VirEigI, VirEigJ,
     $   SpinLabel) )
     $   Call ErrQuit('RIMP2G_V_E2: can''t create amplitude GA', 0,
     &       GA_ERR)
      Call PStat_Off( PS_denom)
c$$$      If ( PrPrgRpt ) Call Progress_Report(LuOut,
c$$$     $   '       T2 amplitudes done')
C
C     Collect statistics on t
C
      Call RIMP2_t_Stats(.FALSE., ILo, IHi, NVirI, JLo, JHi, NVirJ,
     $     G_T2, TopN, TopN_Ind, TopN_Val)
C     
C     Form IJ pair energy
C     
      Call PStat_On(PS_PrEn)
      Call RIMP2_PairE(ILo, IHi, NVirI, JLo, JHi, NVirJ,
     $   .FALSE., IJTri, G_Int, G_T2, G_PairE )
      Call PStat_Off(PS_PrEn)
c$$$      If ( PrPrgRpt ) Call Progress_Report(LuOut,
c$$$     $   '       Pair energies done')
C
C     Don't need integrals or amplitudes anymore
C
      If ( .NOT. KeepInt ) then
         If (.NOT. GA_Destroy(G_Int) ) Call ErrQuit(
     $      'RIMP2_E_IJ: can''t destroy 4c2e integral GA', G_Int,
     &       GA_ERR)
      EndIf
      If ( .NOT. KeepT2 ) then
         If (.NOT. GA_Destroy(G_T2) ) Call ErrQuit(
     $      'RIMP2_E_IJ: can''t destroy T2 GA', G_T2, GA_ERR)
      EndIf
C
      Return
      End
