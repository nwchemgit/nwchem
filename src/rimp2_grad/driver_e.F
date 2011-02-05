C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_E -- Master routine for RI-MP2 energy evaluation
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function RIMP2_DRIVER_E(MaxSpin, BraKetSame,
     $     D_Int, D_Gam,
     $     TopSpin, NMO, NFrzO, NOcc, NVir, NFit, FitBas, Eig, LDEig,
     $     Contrib, E2, Mem_Opt, IBlock, JBlock, DoGrad, G_P2)
      Implicit NONE
#include "errquit.fh"
C
      Integer MaxSpin
      Logical BraKetSame
      Integer TopSpin, NMO(TopSpin), NFit, FitBas,
     $   IBlock, JBlock, LDEig, G_P2(TopSpin)
      Integer NFrzO(TopSpin), NOcc(TopSpin), NVir(TopSpin)
      Double precision Eig(LDEig, TopSpin), Contrib(TopSpin, TopSpin)
      Double precision E2
      Character*(*) Mem_Opt
      Logical DoGrad
      Integer D_Int(MaxSpin, 2, 2), D_Gam(TopSpin)
C
C ARGUMENTS
C     D_Int    Handles for integral files [IN]
C     FNInt    File name for MO basis 3-center integrals [IN]
C     FNGam3   File name for MO basis 3-center gamma [IN]
C     TopSpin  Number of spins in reference (RHF=1, UHF=2) [IN]
C     NMO      Number of MOs [IN]
C     NOcc     Number of occupied orbitals for each spin [IN]
C     NVir     Number of virtual orbitals for each spin [IN]
C     NFit     Number of functions in fitting basis [IN]
C     FitBas   Handle for fitting basis [IN]
C     Eig      Eigenvalues of reference Fock matrix [IN]
C     Contrib  E(2) broken down by spin case. (The ba, or 2,1 element
C              is always zero) [OUT]
C     E2       The second order MBPT energy [OUT]
C     Mem_Opt  Memory allocation optimization mode [IN]
C     IBlock   User-requested blocksize [IN]
C     JBlock   User-requested blocksize [IN]
C     DoGrad   Gradient calculation required [IN]
C     G_P2     1-density [OUT]
C
C RETURN VALUE
C     .TRUE.   If calculation completed successfully.
C     .FALSE.  If an error occurred.
C 
C DESCRIPTION
C     Computes the second order perturbation theory (MBPT(2), MP2)
C     contribution to the correlation energy using 3-center integrals
C     to approximate the usual four-center two-electron integrals
C     in the "V-approximation" of Vahtras, Almlof, and Feyereisen.
C     The 3-center integrals are obtained through calls to the function
C     Get_3ci_block.  Both RHF and UHF canonical references are handled.
C
C     If DoGrad is .TRUE., terms required for the gradient which 
C     explicitly require t are also evaluated.
C
C REFERENCES
C     Feyereisen, Fitzgerald, and Komornicki, Chem. Phys. Lett. 208,
C     359 (1993).
C     Vahtras, Almlof, and Feyereisen, ?
C
C NOTES
C     There are a couple of problems with the present algorithm.
C     Mainly it does things one block at a time.  This is going to be
C     inefficient if one wants to try to do this directly rather than
C     from disk, because we can probably store more than one block at
C     a time in core.  The other point is that the load balance is not
C     necessarily spectacular -- only those nodes which end up with 
C     data on them will be involved in the work.  In particular, since
C     the G_A and G_B matrices are only O(N^2), and we want to keep
C     the block sizes reasonable, they may not be spread out across the
C     whole machine.
C
C MEMORY USE
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "numerical_constants.fh"
#include "stdio.fh"
#include "rimp2g_timers.fh"
#include "util.fh"
#include "inp.fh"
#include "dra.fh"
C
C EXTERNAL ROUTINES
      Logical  GA_Create_One_By_Shells
      External GA_Create_One_By_Shells
C
      external Sum_Pairs
      Double Precision Sum_Pairs
C
      Logical  RIMP2_Mk_T, RIMP2_Asym_I
      External RIMP2_Mk_T, RIMP2_Asym_I
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C PARAMETERS
C
      Integer MinChunk             ! Minimum panel size that is reasonable
      Parameter ( MinChunk = -64)  ! for global distribution of data
      Logical IJTri                ! Triangular arrangement for I and T arrays
      Parameter ( IJTri = .FALSE.) !
      Integer TopN                 ! Report largest TopN amplitudes
      Parameter (TopN = 10)        !
      Integer MyMaxSpin            ! Allows definition of local "adj arrays"
      Parameter (MyMaxSpin = 2)    !
C
C LOCAL VARIBLES
      Integer ISpin, JSpin, JTop, NPass 
      Integer Proc
      Logical Status
      Integer G_bigI, ILo, IHi, NI, G_bigJ, JLo, JHi, NJ, G_Int, G_T2
      Integer G_PairE(2*MyMaxSpin-1)
      Character*255 String1
      Integer NewIBlk, NewJBlk
      Logical PrInfo, PrPrgRpt
      Character*12 SpinLabel
      Integer H_TopN_Ind(2*MyMaxSpin-1), H_TopN_Val(2*MyMaxSpin-1)
      Integer I_TopN_Ind(2*MyMaxSpin-1), I_TopN_Val(2*MyMaxSpin-1)
      Integer ReqI, ReqJ
C
C STATEMENT FUNCTIONS
      Character*1 SpinItoA
      SpinItoA(ISpin) = Char( ISpin + IChar('a') - 1 )
C
      Call GA_Sync
      Call PStat_On(PS_Energy)
C
      PrInfo =   Util_Print('information', Print_Low)
      PrPrgRpt = Util_Print('progress reports', Print_Default)
C
      If ( MaxSpin .gt. MyMaxSpin ) Call ErrQuit(
     $     'RIMP2_Driver_E: fatal program error: MaxSpin > MyMaxSpin',
     $     MaxSpin, CALC_ERR)
C
C     ************************
C     * Check inputs & setup *
C     ************************
C
      If (TopSpin .lt. 1 .OR. TopSpin .gt. MaxSpin) then
         Write (LuErr, *) ' Bad value of TopSpin ', TopSpin
         RIMP2_DRIVER_E = .FALSE.
         Return
      EndIf
      If (LDEig .lt. 1) then
         Write (LuErr, *) ' Bad value of LDEig ', LDEig
         RIMP2_DRIVER_E = .FALSE.
         Return
      EndIf
      If (NOcc(1) .lt. 0 .OR. NOcc(TopSpin) .lt. 0) then
         Write (LuErr, *) ' Bad value of NOcc ', NOcc
         RIMP2_DRIVER_E = .FALSE.
         Return
      EndIf
      If (NVir(1) .lt. 0 .OR. NVir(TopSpin) .lt. 0) then
         Write (LuErr, *) ' Bad value of NVir ', NVir
         RIMP2_DRIVER_E = .FALSE.
         Return
      EndIf
      If (NOcc(1) + NVir(1) .gt. NMO(1)
     $   .OR. NOcc(TopSpin) + NVir(TopSpin) .gt. NMO(TopSpin) ) then
         Write (LuErr, *) ' Bad value of NOcc+NVir ', NOcc, '+', NVir
         RIMP2_DRIVER_E = .FALSE.
         Return
      EndIf
      If (NFit .lt. 1) then
         Write (LuErr, *) ' Bad value of NFit ', NFit
         RIMP2_DRIVER_E = .FALSE.
         Return
      EndIf

C
C     Get my processor number
C
      Proc = GA_NodeID()
C
C     Clear out the energy contribution matrix
C
      Call dlaset('Full', TopSpin, TopSpin, FP_Zero, FP_Zero,
     $     Contrib, TopSpin)
C
      Do ISpin = 1, TopSpin
         Do JSpin = ISpin, TopSpin
            String1 = 'Pair Energies, spin ' // SpinItoA(ISpin)
     $           // SpinItoA(JSpin)
            Status = Status .AND.
     $         GA_Create( MT_Dbl, NOcc(ISpin), NOcc(JSpin),
     $         String1(:Inp_StrLen(String1)), MinChunk, MinChunk,
     $         G_PairE(ISpin + JSpin - 1) )
            Call GA_Zero( G_PairE(ISpin + JSpin - 1 ) )
C
            Call RIMP2_t_Stats_Init(TopN, H_TopN_Ind(ISpin+JSpin-1),
     $           H_TopN_Val(ISpin+JSpin-1) )
            If ( .NOT. MA_Get_Index(H_TopN_Ind(ISpin+JSpin-1),
     $           I_TopN_Ind(ISpin+JSpin-1) ) ) Call ErrQuit(
     $           'RIMP2_Driver_E: problems with t stats arrays',
     $           ISpin+JSpin-1, MA_ERR)
            If ( .NOT. MA_Get_Index(H_TopN_Val(ISpin+JSpin-1),
     $           I_TopN_Val(ISpin+JSpin-1) ) ) Call ErrQuit(
     $           'RIMP2_Driver_E: problems with t stats arrays',
     $           ISpin+JSpin-1, MA_ERR)
         EndDo
      EndDo
C
C     ******************************************
C     * Create some GAs required for gradients *
C     ******************************************
C     P^{(2)}_{pq}
C
      If ( DoGrad ) then
         Do ISpin = 1, TopSpin
            String1 = 'P(2) spin ' // SpinItoA(ISpin)
            If (.NOT. GA_Create(MT_Dbl, NMO(ISpin), NMO(ISpin),
     $         String1(:Inp_StrLen(String1)), MinChunk, MinChunk,
     $         G_P2(ISpin) ) ) Call ErrQuit(
     $         'RIMP2_Driver_E: can''t allocate P(2)',
     $         NMO(ISpin)*NMO(ISpin) , GA_ERR)
            Call GA_Zero( G_P2(ISpin) )
         EndDo
      EndIf ! Gradient-only
C
C     ***************************
C     * Determine loop blocking *
C     ***************************
C
      Call RIMP2G_E_Mem(TopSpin, NOcc, NVir, FitBas, Mem_Opt,
     $   NewIBlk, NewJBlk)
      If ( IBlock .le. 0) then
         IBlock = NewIBlk
      Else
         IBlock = Min(IBlock, NewIBlk)
      EndIf
      If ( JBlock .le. 0) then
         JBlock = NewJBlk
      Else
         JBlock = Min(JBlock, NewJBlk)
      EndIf
C
      If ( (PrInfo .AND. GA_NodeID() .eq. 0) .or. PrPrgRpt ) then
         If ( DoGrad ) then
            Write (String1, 9105) IBlock, JBlock
         Else
            Write (String1, 9100) IBlock, JBlock
         EndIf
         Call Util_Sqz_Spc(String1)
         Write (LuOut, 9110) String1(:Inp_StrLen(String1))
      EndIf
 9100 Format('Energy evaluation using block size ', I5, ' x ', I5)
 9105 Format('Energy evaluation using block size ', I5, ' x ', I5,
     $     ' (gradient terms also evaluated)')
 9110 Format(1X, A)
C
C     *****************************************************************
C     * Loop over occupied orbitals of both spins for the first index *
C     *****************************************************************
C
      NPass = 0
      Do ISpin = 1, TopSpin
C
         Do ILo = 1, NOcc(ISpin), IBlock
            IHi = Min( NOcc(ISpin), ILo + IBlock - 1)
            NI = IHi - ILo + 1
C
C           *************************
C           * Second occupied index *
C           *************************
C           Should cover both spins if I is alpha, giving the pure 
C           alpha and mixed spin components, but only beta Js are 
C           required when I is beta -- only the pure beta 
C           contribution is left.
C
            Do JSpin = ISpin, TopSpin
               If ( TopSpin .gt. 1) then
                  SpinLabel = 'spin ' // SpinItoA(ISpin)
     $               // SpinItoA(JSpin)
               Else
                  Spinlabel = 'spin adapted'
               EndIf
C
C              Use permutational symmetry if possible
C
c$$$               If ( ISpin .eq. JSpin) then
c$$$                  JTop = IHi
c$$$               Else
c$$$                  JTop = NOcc(JSpin)
c$$$               EndIf
               JTop = NOcc(JSpin)
C              
               Do JLo = 1, JTop, JBlock
                  JHi = Min( JTop, JLo + JBlock - 1)
                  NJ = JHi - JLo + 1
C
                  NPass = NPass + 1
                  If ( PrPrgRpt ) then
                     String1 = ' ' ! Not all machines initialize mem.
                     If ( TopSpin .eq. 1) then
                        Write (String1, 9201) NPass, ILo, IHi, JLo, JHi
                     Else
                        Write (String1, 9202) NPass,
     $                     SpinLabel(:Inp_StrLen(SpinLabel)),
     $                     ILo, IHi, JLo, JHi
                     EndIf
                     Call Util_Sqz_Spc(String1)
                     String1 = '    '//String1
                     Call Progress_Report(LuOut,
     $                  String1(:Inp_StrLen(String1)))
 9201                Format('Beginning pass ', I5,
     $                  '(', I5, ':', I5, ',', I5, ':', I5, ')')
 9202                Format('Beginning pass ', I5, 1X, A, 1X,
     $                  '(', I5, ':', I5, ',', I5, ':', I5, ')')
                  EndIf
C
C                 Get the I and J integrals.  Use the NFrzO offsets
C                 to make sure we get the right ones
C
                  Call RIMP2_IJ_Int(ISpin .eq. JSpin .AND.
     $                 D_Int(ISpin, 1, 1) .eq. D_Int(JSpin, 2, 1),
     $                 D_Int(ISpin, 1, 1), ILo+NFrzO(ISpin),
     $                 IHi+NFrzO(ISpin), NVir(ISpin),
     $                 D_Int(JSpin, 2, 1), JLo+NFrzO(JSpin),
     $                 JHi+NFrzO(JSpin), NVir(JSpin),
     $                 FitBas, NFit, G_bigI, ReqI, G_bigJ, ReqJ)
C
C                 Do energy stuff for this IJ batch.  
C                 Keep 3c integrals & t if we're doing gradients.
C                 Don't keep 4c integrals ever.
C
                  If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $                 '       Beginning energy terms')
                  Call RIMP2_E_IJ(TopSpin, ISpin .eq. JSpin,
     $                 ILo, IHi, NVir(ISpin), JLo, JHi, NVir(JSpin),
     $                 Eig(1+NFrzO(ISpin), ISpin),
     $                 Eig(1+NFrzO(JSpin), JSpin),
     $                 Eig(1+NFrzO(ISpin)+NOcc(ISpin), ISpin),
     $                 Eig(1+NFrzO(JSpin)+NOcc(JSpin), JSpin), 
     $               G_bigI, DoGrad, G_bigJ, DoGrad,
     $               G_Int, .FALSE., G_T2, DoGrad,
     $               G_PairE(ISpin + JSpin - 1),
     $                 SpinLabel(:Inp_StrLen(SpinLabel)),
     $                 TopN, Int_MB(I_TopN_Ind(ISpin+JSpin-1)),
     $                 Dbl_MB(I_TopN_Val(ISpin+JSpin-1)), ReqI, ReqJ  )
C
C                 Do gradient stuff for this IJ batch.  3c integrals,
C                 t no longer needed.
C
                  If ( DoGrad ) then
                     If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $                    '       Beginning density/gradient terms')
                     Call RIMP2_G_IJ(TopSpin, ISpin .eq. JSpin,
     $                    NFrzO(ISpin), ILo, IHi, NVir(ISpin),
     $                    NFrzO(JSpin), JLo, JHi, NVir(JSpin),
     $                    G_bigI, .FALSE., G_bigJ, .FALSE.,
     $                    G_T2, .FALSE., G_PairE(ISpin + JSpin - 1),
     $                    NOcc(ISpin), ISpin, NOcc(JSpin), JSpin, G_P2,
     $                    SpinLabel(:Inp_StrLen(SpinLabel)), D_Gam,
     $                    BraKetSame, D_Int(JSpin, 1, 1) )
C     
                  EndIf         ! Gradient-only section
C                 
               EndDo            ! J (blocks)
            EndDo               ! JSpin
         EndDo                  ! I (blocks)
      EndDo                     ! ISpin
C
C     Print amplitude statistics if requested & clean up
C
      If ( TopN .gt. 0) then
         Write (LuOut, '(A)') ' '
         Call Banner(LuOut, 'Largest t Amplitudes', '-', .TRUE.,
     $        .TRUE., .FALSE.)
      EndIf
      Do ISpin = 1, TopSpin
         Do JSpin = ISpin, TopSpin
C
            If ( TopSpin .gt. 1) then
               SpinLabel = SpinItoA(ISpin) // SpinItoA(JSpin)
            Else
               Spinlabel = 'adapted'
            EndIf
C
            Call RIMP2_t_Stats_Merge(TopN,
     $           Int_MB(I_TopN_Ind(ISpin+JSpin-1) ),
     $           Dbl_MB(I_TopN_Val(ISpin+JSpin-1) ) )
            Call RIMP2_t_Stats_Print(TopN,
     $           Int_MB(I_TopN_Ind(ISpin+JSpin-1) ),
     $           Dbl_MB(I_TopN_Val(ISpin+JSpin-1) ),
     $           0, 0, NOcc(ISpin), NOcc(JSpin),
     $           SpinLabel(:Inp_StrLen(SpinLabel)) )
C
            Call RIMP2_t_Stats_Term(TopN, H_TopN_Ind(ISpin+JSpin-1),
     $           H_TopN_Val(ISpin+JSpin-1) )
         EndDo
      EndDo
C
C     ***********************************
C     * Compute total energy from pairs *
C     ***********************************
C
      Call RIMP2_E_Final(TopSpin, G_PairE, Contrib, E2)
C
      Do ISpin = 1, TopSpin
         Do JSpin = ISpin, TopSpin
            Status = GA_Destroy( G_PairE(ISpin+JSpin-1) )
         EndDo
      EndDo
C
      Call GA_Sync
      Call PStat_Off(PS_Energy)
C
      RIMP2_DRIVER_E = .TRUE.
      Return
      End
