C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Driver_E -- Master routine for RI-MP2 energy evaluation
C
C REVISION
C     $Id: driver_e.F,v 1.2 1998-07-13 02:38:07 bernhold Exp $
C
C SYNOPSIS
      Logical Function RIMP2_DRIVER_E(MaxSpin, BraKetSame, FN_Int,
     $     D_Int, FN_Gam, D_Gam,
     $     TopSpin, NMO, NFrzO, NOcc, NVir, NFit, FitBas, Eig, LDEig,
     $     Contrib, E2, Mem_Opt, IBlock, JBlock, DoGrad, G_P2)
      Implicit NONE
C
      Integer MaxSpin
      Logical BraKetSame
      Integer TopSpin, NMO(TopSpin), NFit, FitBas,
     $   IBlock, JBlock, LDEig, G_P2(TopSpin)
      Integer NFrzO, NOcc(TopSpin), NVir(TopSpin)
      Double precision Eig(LDEig, TopSpin), Contrib(TopSpin, TopSpin)
      Double precision E2
      Character*(*) Mem_Opt
      Logical DoGrad
      Integer D_Int(MaxSpin, 2, 2), D_Gam(TopSpin)
      Character*(*) FN_Int(MaxSpin, 2, 2), FN_Gam(TopSpin)
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
C     FN_Int   Integral filenames [IN]
C     FN_Gam   Gamma3 filenames [IN]
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
      Logical Get_3cInts, GA_Create_One_By_Shells
      External Get_3cInts, GA_Create_One_By_Shells
C
      external Sum_Pairs
      Double Precision Sum_Pairs
      External GA_DDot_Patch_DP
      Double precision GA_DDot_Patch_DP
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
      Parameter (TopN = 10)
C
C LOCAL VARIBLES
      Integer I, J, ISpin, JSpin, JTop, NPass 
      Integer Proc
      Logical Status
      Integer G_bigI, ILo, IHi, NI, G_bigJ, JLo, JHi, NJ, G_Int, G_T2
      Integer G_PairE(2*MaxSpin-1)
      Character*255 String1, TmpName, FilNam
      Integer NewIBlk, NewJBlk
      Logical PrInfo, PrPrgRpt
      Character*1 Sp
      Character*12 SpinLabel
      Integer H_TopN_Ind(2*MaxSpin-1), H_TopN_Val(2*MaxSpin-1)
      Integer I_TopN_Ind(2*MaxSpin-1), I_TopN_Val(2*MaxSpin-1)
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
     $           ISpin+JSpin-1)
            If ( .NOT. MA_Get_Index(H_TopN_Val(ISpin+JSpin-1),
     $           I_TopN_Val(ISpin+JSpin-1) ) ) Call ErrQuit(
     $           'RIMP2_Driver_E: problems with t stats arrays',
     $           ISpin+JSpin-1)
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
     $         NMO(ISpin)*NMO(ISpin) )
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
C     **************
C     * Open files *
C     **************
C     Open both bra and ket files if they are distinct, otherwise
C     simply duplicate the bra handle.
C
      Do ISpin = 1, TopSpin
C
C        Integrals
C
         If ( DRA_Open(
     $        FN_Int(ISpin, 1, 1)(:Inp_StrLen(FN_Int(ISpin, 1, 1))),
     $        DRA_R, D_Int(ISpin, 1, 1) ) .ne. 0) then 
            Write (LuErr, *) 'RIMP2_Driver_E:  Problems opening ',
     $         'integral file for spin ', ISpin, '!'
            RIMP2_Driver_E = .FALSE.
            Return
         EndIf
         If ( BraKetSame ) then
            D_Int(ISpin, 2, 1) = D_Int(ISpin, 1, 1)
         Else
            If ( DRA_Open(
     $           FN_Int(ISpin, 2, 1)(:Inp_StrLen(FN_Int(ISpin, 2, 1))),
     $           DRA_R, D_Int(ISpin, 2, 1) ) .ne. 0) then 
               Write (LuErr, *) 'RIMP2_Driver_E:  Problems opening ',
     $              'integral file for spin ', ISpin, '!'
               RIMP2_Driver_E = .FALSE.
               Return
            EndIf
         EndIf
         If ( DoGrad ) then
            If ( DRA_Create(MT_Dbl, NOcc(ISpin)*NVir(ISpin), NFit,
     $           'MO 3-center gamma',
     $           FN_Gam(ISpin)(:Inp_StrLen(FN_Gam(ISpin))),
     $           DRA_RW, NVir(ISpin), NFit, D_Gam(ISpin) ) .ne. 0) then
               Write (LuErr, *) 'RIMP2_Driver_E:  Problems creating ',
     $              'gamma3 file for spin ', ISpin, '!'
               RIMP2_Driver_E = .FALSE.
               Return
            EndIf ! DRA_Create
         EndIf    ! DoGrad
      EndDo       ! ISpin
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
C                 Get the I and J integrals
C
                  Call RIMP2_IJ_Int(ISpin .eq. JSpin .AND.
     $                 D_Int(ISpin, 1, 1) .eq. D_Int(JSpin, 2, 1),
     $                 D_Int(ISpin, 1, 1), ILo, IHi, NVir(ISpin),
     $                 D_Int(JSpin, 2, 1), JLo, JHi, NVir(JSpin),
     $                 FitBas, NFit, G_bigI, ReqI, G_bigJ, ReqJ)
C
C                 Do energy stuff for this IJ batch.  
C                 Keep 3c integrals & t if we're doing gradients.
C                 Don't keep 4c integrals ever.
C
                  Call RIMP2_E_IJ(TopSpin, ISpin .eq. JSpin,
     $               ILo, IHi, NVir(ISpin), JLo, JHi, NVir(JSpin),
     $               Eig(1, ISpin), Eig(1, JSpin),
     $               Eig(1+NOcc(ISpin), ISpin),
     $               Eig(1+NOcc(JSpin), JSpin), 
     $               G_bigI, DoGrad, G_bigJ, DoGrad,
     $               G_Int, .FALSE., G_T2, DoGrad,
     $               G_PairE(ISpin + JSpin - 1),
     $                 SpinLabel(:Inp_StrLen(SpinLabel)),
     $                 TopN, Int_MB(I_TopN_Ind(ISpin+JSpin-1)),
     $                 Dbl_MB(I_TopN_Val(ISpin+JSpin-1)), ReqI, ReqJ  )
                  If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $                 '       Energy terms done')
C
C                 Do gradient stuff for this IJ batch.  3c integrals,
C                 t no longer needed.
C
                  If ( DoGrad ) then
                     Call RIMP2_G_IJ(TopSpin, ISpin .eq. JSpin,
     $                  ILo, IHi, NVir(ISpin), JLo, JHi, NVir(JSpin),
     $                  Eig(1, ISpin), Eig(1, JSpin),
     $                  Eig(1+NOcc(ISpin), ISpin),
     $                  Eig(1+NOcc(JSpin), JSpin), 
     $                  G_bigI, .FALSE., G_bigJ, .FALSE.,
     $                  G_T2, .FALSE., G_PairE(ISpin + JSpin - 1),
     $               NOcc(ISpin), ISpin, NOcc(JSpin), JSpin, G_P2,
     $               SpinLabel(:Inp_StrLen(SpinLabel)), D_Gam,
     $                    BraKetSame, D_Int(JSpin, 1, 1) )
                  If ( PrPrgRpt ) Call Progress_Report(LuOut,
     $                 '       Gradient terms done')
C
                  EndIf         ! Gradient-only section
C                 
               EndDo            ! J (blocks)
            EndDo               ! JSpin
         EndDo                  ! I (blocks)
      EndDo                     ! ISpin
C
C     ************************
C     * Close & delete files *
C     ************************
C
      Do ISpin = 1, TopSpin
C
C        If we're doing gradients, we still need to keep them,
C        otherwise we can delete them entirely.
C
         If ( DoGrad ) then
            If ( DRA_Close(D_Int(ISpin, 1, 1)) .ne. 0) then
               Write (LuErr, *) 'RIMP2_DRIVER_E:  Problems closing ',
     $              'integral file(s)!'
               RIMP2_DRIVER_E = .FALSE.
               Return
            EndIf
            If ( .NOT. BraKetSame ) then
               If ( DRA_Close(D_Int(ISpin, 2, 1)) .ne. 0) then
                  Write (LuErr, *) 'RIMP2_DRIVER_E:  Problems closing ',
     $                 'integral file(s)!'
                  RIMP2_DRIVER_E = .FALSE.
                  Return
               EndIf
            EndIf
            If ( DRA_Close(D_Gam(ISpin)) .ne. 0) then
               Write (LuErr, *) 'RIMP2_DRIVER_E:  Problems closing ',
     $              'gamma3 file(s)!'
               RIMP2_DRIVER_E = .FALSE.
               Return
            EndIf
         Else       ! .NOT. DoGrad
            If ( DRA_Delete(D_Int(ISpin, 1, 1)) .ne. 0) then
               Write (LuErr, *) 'RIMP2_DRIVER_E:  Problems closing ',
     $              'integral file(s)!'
               RIMP2_DRIVER_E = .FALSE.
               Return
            EndIf
            If ( .NOT. BraKetSame ) then
               If ( DRA_Delete(D_Int(ISpin, 2, 1)) .ne. 0) then
                  Write (LuErr, *) 'RIMP2_DRIVER_E:  Problems closing ',
     $                 'integral file(s)!'
                  RIMP2_DRIVER_E = .FALSE.
                  Return
               EndIf
            EndIf
         EndIf ! DoGrad
      EndDo    ! ISpin
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
     $           NFrzO, NFrzO, NOcc(ISpin), NOcc(JSpin),
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
C     ****************************************
C     * Finished with gradient intermediates *
C     ****************************************
C
      If ( DoGrad ) then
         Do ISpin = 1, TopSpin
            If (.NOT. GA_Destroy( G_P2(ISpin) ) ) Call ErrQuit(
     $         'RIMP2_Driver_E: can''t destroy P2 ', 0)
         EndDo
      EndIf
C
      Call GA_Sync
      Call PStat_Off(PS_Energy)
C
      RIMP2_DRIVER_E = .TRUE.
      Return
      End
