C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_G_IJ -- Evaluate RI-MP2 gradient contributions for IJ batch
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_G_IJ(TopSpin, SameSpin, NFrzI, ILo, IHi, NVirI,
     $     NFrzJ, JLo, JHi, NVirJ, G_I, KeepI, G_J, KeepJ, G_T2,
     $     KeepT2, G_PairE, NoccI, ISpin, NOccJ, JSpin, G_P2,
     $     SpinLabel, D_Gam3, BraKetSame, D_Int)
      Implicit NONE
#include "errquit.fh"
C
      Integer TopSpin  ![in] RHF=1 UHF=2
      Logical SameSpin ![in] I and J have same spin
      Integer NFrzI    ![in] Nr of frozen occupied I orbitals
      Integer ILo      ![in] Lower limit of I in current batch
      Integer IHi      ![in] Upper limit of I in current batch
      Integer NVirI    ![in] Number of virtuals in I spin case
      Integer NFrzJ    ![in] Nr of frozen occupied J orbitals
      Integer JLo      ![in] Lower limit of J in current batch
      Integer JHi      ![in] Upper limit of J in current batch
      Integer NVirJ    ![in] Number of virtuals in J spin case
      Integer G_I      ![in]  3c2e integrals for I
      Logical KeepI    ![in]  Keep G_I around for caller
      Integer G_J      ![in]  3c2e integrals for J
      Logical KeepJ    ![in]  Keep G_J around for caller
c$$$      Integer G_Int    ![out] (approx) antisym 4c2e integral array (out)
c$$$      Logical KeepInt  ![in]  Keep G_Int around for caller
      Integer G_T2     ![out] T2 amplitude array (out)
      Logical KeepT2   ![in]  Keep G_T2 around for caller
      Integer G_PairE  ![in]  Pair energies for this spin case (inout)
C
      Integer NOccI
      Integer ISpin
      Integer NOccJ
      Integer JSpin
      Integer G_P2(TopSpin)
      Character*(*) SpinLabel
      Integer D_Gam3(TopSpin)
      Logical BraKetSame
      Integer D_Int
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
C
C INCLUDED FILES
#include "global.fh"
#include "mafdecls.fh"
#include "dra.fh"
#include "util.fh"
#include "rimp2g_timers.fh"
#include "stdio.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Integer Chunk               ! chunk size for GA distr. of G_Int, G_T2
      Parameter (Chunk = -1)      ! Chosen to give even distr.
      Logical IJTri               ! G_Int, G_T2 contain only unique IJ
      Parameter (IJTri = .FALSE.) ! Keep them rectangular
C
C LOCAL VARIABLES
      Integer NI, NJ, NFit, Type, Junk
      Integer Vir1, Occ1
      Logical PrPrgRpt
C
      NI = IHi - ILo + 1
      NJ = JHi - JLo + 1
C
c$$$      PrPrgRpt = Util_Print('progress reports', Print_Default)
C
C     Compute P^{(2)}_{ab} contributions
C
      Vir1 = NFrzI+NOccI+1
      Call RIMP2_P2AB(SameSpin, .FALSE., IJTri, ILo, IHi, NVirI,
     $     JLo, JHi, NVirJ, .FALSE., G_T2, G_P2(ISpin), Vir1)
C
C     If this is the mixed spin case, we must also get the contribution
C     to the beta density.
C
      If ( TopSpin .eq. 2 .AND. .NOT. SameSpin) then
         Vir1 = NFrzJ+NOccJ+1
         Call RIMP2_P2AB(Samespin, .FALSE., IJTri, ILo, IHi, NVirI,
     $        JLo, JHi, NVirJ, .TRUE., G_T2, G_P2(JSpin), Vir1)
      EndIf
c$$$      If ( PrPrgRpt ) Call Progress_Report(LuOut,
c$$$     $   '       P(2)_{ab} done')
C
C     Compute P^{(2)}_{ij} contributions
C
      Occ1 = NFrzI+1
      Call RIMP2_P2ij(SameSpin, .FALSE., IJTri, ILo, IHi, NVirI,
     $     JLo, JHi, NVirJ, .FALSE., G_T2, G_P2(ISpin), Occ1)
C
C     If this is the mixed spin case, we must also get the contribution
C     to the beta density.
C
      If ( TopSpin .eq. 2 .AND. .NOT. SameSpin) then
         Occ1 = NFrzJ+1
         Call RIMP2_P2ij(SameSpin, .FALSE., IJTri, ILo, IHi, NVirI,
     $        JLo, JHi, NVirJ, .TRUE., G_T2, G_P2(JSpin), Occ1)
      EndIf
c$$$      If ( PrPrgRpt ) Call Progress_Report(LuOut,
c$$$     $   '       P(2)_{ij} done')
C
C     Compute G(D,ia) = (D|jb) t(ij,ab) aka gamma3
C
C     This is an intermediate in L1(pi), L2(ap), and the contraction
C     of the 2-PDM (gamma) with derivative integrals.
C
      Call RIMP2_Gamma3(.FALSE., IJTri, ILo, IHi, NVirI, JLo, JHi,
     $   NVirJ, G_I, G_J, G_T2, D_Gam3(ISpin), .FALSE., D_Int )
C
C     If this is the mixed spin case, we must also get the contribution
C     to the beta version.  In case we have to read integrals, shift
C     J range to account for frozen orbitals
C
      If ( TopSpin .eq. 2 .AND. .NOT. SameSpin) then
C
         Call RIMP2_Gamma3(.FALSE., IJTri, JLo, JHi, NVirJ,
     $        NFrzI+ILo, NFrzI+IHi, NVirI, G_J, G_I, G_T2,
     $        D_Gam3(JSpin), .NOT. BraKetSame, D_Int)
      EndIf
c$$$      If ( PrPrgRpt ) Call Progress_Report(LuOut,
c$$$     $   '       G(D,ia) intermediate done')
C
C
C     Don't need integrals or amplitudes anymore
C
      If ( .NOT. KeepI ) then
         If (.NOT. GA_Destroy(G_I) ) Call ErrQuit(
     $      'RIMP2_G_IJ: can''t destroy 3c2e integral GA', G_I, GA_ERR)
      EndIf
      If ( .NOT. KeepJ ) then
         If (.NOT. GA_Destroy(G_J) ) Call ErrQuit(
     $      'RIMP2_G_IJ: can''t destroy 3c2e integral GA', G_J, GA_ERR)
      EndIf
c$$$      If ( .NOT. KeepInt ) then
c$$$         If (.NOT. GA_Destroy(G_Int) ) Call ErrQuit(
c$$$     $      'RIMP2_G_IJ: can''t destroy 4c2e integral GA', G_Int)
c$$$      EndIf
      If ( .NOT. KeepT2 ) then
         If (.NOT. GA_Destroy(G_T2) ) Call ErrQuit(
     $      'RIMP2_G_IJ: can''t destroy T2 GA', G_T2, GA_ERR)
      EndIf
C
      Return
      End
