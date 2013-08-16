C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_v_e2 -- Evaluate MP2 energy using V-approximation 3 ctr ints.
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function RIMP2_V_E2( rtdb, D_Int, FNInt, TopSpin, NMO,
     $   NOcc, NVir, NFit, FitBas, Eig, LDEig, Contrib, E2,
     $   Mem_Opt, IBlock, JBlock,actualname)
      Implicit NONE
C
      Character*(*) FNInt
      Integer TopSpin, D_Int(TopSpin), NMO(TopSpin), NFit, FitBas,
     $   IBlock, JBlock, LDEig, rtdb
      Integer NOcc(TopSpin), NVir(TopSpin)
      Double precision Eig(LDEig, TopSpin), Contrib(TopSpin, TopSpin)
      Double precision E2
      Character*(*) Mem_Opt,actualname
C
C ARGUMENTS
C     rtdb     Handle for run-time database [IN]
C     D_Int    Handles for integral files [IN]
C     FNInt    File name for MO basis 3-center integrals [IN]
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
C     Global: 2*NVir*NFit + NVir*NVir
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "numerical_constants.fh"
#include "stdio.fh"
#include "rimp2_timers.fh"
#include "util.fh"
#include "rtdb.fh"
C
C EXTERNAL ROUTINES
      Logical Open_3ci_File, Get_3cInts, Close_3cI_File,
     $  GA_Create_One_By_Shells
      External Open_3cI_File, Get_3cInts, Close_3cI_File,
     $  GA_Create_One_By_Shells
C
      external Sum_Pairs
      Double Precision Sum_Pairs
      Integer LNBlnk
      External LNBlnk
#ifdef BAD_GACCESS
      External rimp2_ddotp
      Double precision rimp2_ddotp
#endif
c$$$C
c$$$      External GA_DASum_Patch, GA_DAMax_Patch
c$$$      Double precision GA_DASum_Patch, GA_DAMax_Patch
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C PARAMETERS
C
      Integer MinChunk             ! Minimum panel size that is reasonable
      Parameter ( MinChunk = -64)  ! for global distribution of data
C
C LOCAL VARIBLES
      Integer I, J, ISpin, JSpin, JTop, JTop2
      Integer Proc
      Double Precision Fact, PairE
      Logical Status, PermSym
      Integer G_bigI, ILo, IHi, NI, G_bigJ, JLo, JHi, NJ, G_bigA, G_BigB
      Integer HH_PairE, IH_PairE
c$$$      Double precision ChkSum, MxElem
      Character*255 String1
      Integer NewIBlk, NewJBlk
      Logical PrInfo
c
      double precision fss, fos ! Spin Component Scaled (SCS) MP2
      logical scs               ! SCS flag
C
      Call GA_Sync
      Call PStat_On(PS_Energy)
C
      PrInfo = Util_Print('information', Print_Low)
C
C     ************************
C     * Check inputs & setup *
C     ************************
C
      If (TopSpin .lt. 1 .OR. TopSpin .gt. 2) then
         Write (LuErr, *) ' Bad value of TopSpin ', TopSpin
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf
      If (LDEig .lt. 1) then
         Write (LuErr, *) ' Bad value of LDEig ', LDEig
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf
      If (NOcc(1) .lt. 0 .OR. NOcc(TopSpin) .lt. 0) then
         Write (LuErr, *) ' Bad value of NOcc ', NOcc
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf
      If (NVir(1) .lt. 0 .OR. NVir(TopSpin) .lt. 0) then
         Write (LuErr, *) ' Bad value of NVir ', NVir
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf
      If (NOcc(1) + NVir(1) .gt. NMO(1)
     $   .OR. NOcc(TopSpin) + NVir(TopSpin) .gt. NMO(TopSpin) ) then
         Write (LuErr, *) ' Bad value of NOcc+NVir ', NOcc, '+', NVir
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf
      If (NFit .lt. 1) then
         Write (LuErr, *) ' Bad value of NFit ', NFit
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf

c     flag and scaling factors for SCS-MP2
      if (.not. rtdb_get(rtdb, 'mp2:scs', mt_log, 1, scs))
     $  scs = .false.       ! flag for SCS
      if (scs) then
         if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $     fss = 1.0d0/3.0d0   ! same spin scaling factor for SCS-MP2
         if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $     fos = 1.2d0         ! opposite spin scaling factor for SCS-MP2
      else
         fss = 1.0d0
         fos = 1.0d0
      endif


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
      Status = .TRUE.
      Status = Status .AND. MA_Push_Get(MT_Int, 2*TopSpin-1,
     $   'Pair energy handles', HH_PairE, IH_PairE)
      If ( .NOT. Status ) then
         Write (LuErr, *) ' Unable to allocate Pair energy handles'
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf

C
      Do ISpin = 1, TopSpin
         Do JSpin = ISpin, TopSpin
            Status = Status .AND.
     $         GA_Create( MT_Dbl, NOcc(ISpin), NOcc(JSpin),
     $         'Pair Energies', MinChunk, MinChunk,
     $         Int_MB(IH_PairE + ISpin + JSpin - 2) )
            If ( .NOT. Status ) then
              Write (LuErr, *) ' Unable to allocate Pair energy handles'
              RIMP2_V_E2 = .FALSE.
              Return
            EndIf
            Call GA_Zero( Int_MB(IH_PairE + ISpin + JSpin - 2 ) )
         EndDo
      EndDo
C
C     ***************************
C     * Determine loop blocking *
C     ***************************
C
      Call RIMP2_E_Mem(TopSpin, NOcc, NVir, FitBas, Mem_Opt,
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
      If ( PrInfo .AND. GA_NodeID() .eq. 0)
     $   Write (LUOut, 9100) IBlock, JBlock
 9100 Format(1X, 'Energy evaluation using block sizes ', 2I5)
C
C     *****************************************************************
C     * Loop over occupied orbitals of both spins for the first index *
C     *****************************************************************
C
      Do ISpin = 1, TopSpin
C
C        Create the global array for the I block
C
         Status = GA_Create_One_By_Shells( MT_Dbl,
     $      Min( NOcc(ISpin), IBlock) * NVir(ISpin), FitBas,
     $      'Integral matrix 1', G_bigI)
C
         If ( .NOT. Status ) then
            Write (LuErr, *) ' Unable to allocate for big I block'
            RIMP2_V_E2 = .FALSE.
            Return
         EndIf
C
C        Open integral file
C
         Status = Open_3ci_file(D_Int(ISpin), FNInt, ISpin,
     $      NVir(ISpin), NOcc(ISpin), G_bigI,actualname)
#ifdef DEBUG
         write(luout,*) ga_nodeid(),' ve2 opened FD ',D_Int(ISpin)
         call util_flush(luout)
#endif
         If ( .NOT. Status ) then
            Write (LuErr, *) 'RIMP2_v_E2:  Problems opening ',
     $         'integral file for spin ', ISpin, '!'
            RIMP2_v_E2 = .FALSE.
            Return
         EndIf
C
         Do ILo = 1, NOcc(ISpin), IBlock
            IHi = Min( NOcc(ISpin), ILo + IBlock - 1)
            NI = IHi - ILo + 1
C
C           Get integrals (ia|s') for this I
C
            Call PStat_On( PS_Rd_3ci)
            If ( .NOT. Get_3cInts(D_Int(ISpin), G_bigI, NVir(ISpin),
     $         ILo, IHi) ) then
               Write( LuErr, *) ' Unable to get integral block ',
     $            I, ISpin
               RIMP2_V_E2 = .FALSE.
               Return
            EndIf
            Call PStat_Off( PS_Rd_3ci)
C
            If ( Util_Print('3c ints', Print_Debug) )
     $         Call GA_Print( G_bigI )
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
C
C              Open integral file if JSpin != ISpin
C              NOTE: This routine gets the distribution of the fitting
C              basis from the GA that's passed in.  We don't create the
C              real target array for these guys until a little later,
C              but we can use G_bigI here because the fitting basis will
C              not change with the spin case.
C
               If ( JSpin .ne. ISpin) then
                  Status = Open_3ci_file(D_Int(JSpin), FNInt, JSpin,
     $               NVir(JSpin), NOcc(JSpin), G_bigI,actualname)
                  If ( .NOT. Status ) then
                     Write (LuErr, *) 'RIMP2_v_E2:  Problems opening ',
     $                  'integral file for spin ', JSpin, ' !'
                     RIMP2_v_E2 = .FALSE.
                     Return
                  EndIf
               EndIf
C
C              Use permutational symmetry if possible
C
               If ( ISpin .eq. JSpin) then
                  JTop = IHi
                  PermSym = .TRUE.
               Else
                  JTop = NOcc(JSpin)
                  PermSym = .FALSE.
               EndIf
C              
               Do JLo = 1, JTop, JBlock
                  JHi = Min( JTop, JLo + JBlock - 1)
                  NJ = JHi - JLo + 1
C
C                 Create GAs for result of multiplication and
C                 second set of integrals (J).  Do the result matrix
C                 first so GA's are treated like a stack to avoid
C                 a big hole when the J is freed and we continue using A.
C
                  If ( .NOT. GA_Create( MT_Dbl, NI*NVir(ISpin),
     $               NJ*NVir(JSpin), 'Temporary matrix big A', MinChunk,
     $               MinChunk, G_bigA) ) then
                     Write (LuErr, *) ' Unable to allocate big A'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
C
                  Status = GA_Create_One_By_Shells( MT_Dbl,
     $               NJ * NVir(JSpin), FitBas,
     $               'Integral matrix 2', G_bigJ)
C
                  If ( .NOT. Status ) then
                     Write (LuErr, *) ' Unable to allocate integral ',
     $                  ' block ', J, JSpin
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
C
C                 Find the integrals (jb|s') for this J.  Try to use
C                 the I patch if possible.  For simplicity, we
C                 must get the _whole_ patch from I because Get_3cInts
C                 doesn't understand GA patches.
C                 
                  If ( ISpin .eq. JSpin .AND. JLo .ge. ILo
     $               .AND. JHi .le. IHi) then
                     Call GA_Sync
#ifdef BAD_GACCESS
                     Call rimp2_copypp('N',
#else
                     Call GA_Copy_Patch('N',
#endif
     $                  G_BigI, (JLo-ILo)*NVir(JSpin)+1,
     $                          (JHi-ILo+1)*NVir(JSpin), 1, NFit,
     $                  G_bigJ, 1, NJ*NVir(JSpin),       1, NFit )
                     Call GA_Sync
                  Else
                     Call PStat_On( PS_Rd_3ci)
                     If ( .NOT. Get_3cInts(D_Int(JSpin), G_bigJ,
     $                  NVir(JSpin), JLo, JHi) ) then
                        Write( LuErr, *)
     $                     ' Unable to get integral block ',
     $                     JLo, JSpin
                        RIMP2_V_E2 = .FALSE.
                        Return
                     EndIf
                     Call PStat_Off( PS_Rd_3ci)
                  EndIf
C
                  If ( Util_Print('3c ints', Print_Debug) )
     $               Call GA_Print( G_bigJ )
C
C                 ******************************************************
C                 * Form a block of four-center two-electron integrals *
C                 ******************************************************
C                 from these three-center things:
C                    (ia|jb) = (ia|s')(s'|jb)
C                 where the |s') = V^{-1/2} |s) has been done
C                 as part of the transformation.
C
                  Call GA_Sync
                  Call PStat_On( PS_form_4ci)
                  Call GA_DGEMM('N', 'T',
     $               NI*NVir(ISpin), NJ* NVir(JSpin), NFit,
     $               FP_One, G_bigI, G_bigJ, FP_Zero, G_bigA)
                  Call PStat_Off( PS_form_4ci)
C
                  If ( Util_Print('progress reports', Print_Default))
     $               then
                     String1 = ' ' ! Note all machines initialize mem.
                     Write (String1, 9210) ILo, IHi, JLo, JHi
                     Call Progress_Report(LuOut,
     $                  String1(:LNBlnk(String1)))
 9210                Format('RIMP2_V_E2: Approx. ints. done (',
     $                  I5, ':', I5, ',', I5, ':', I5, ')')
                  EndIf
C
                  If ( Util_Print('4c ints', Print_Debug) )
     $               Call GA_Print( G_bigA )
C
C                 The J integral block is no longer needed and can
C                 now be freed.  This alows us to reuse the space
C                 to hold the temporary matrix B for the other
C                 half of this iteration.
C
                  If ( .NOT. GA_Destroy( G_bigJ) ) then
                     Write (LuErr, *) ' Unable to free integral block ',
     $                  JLo, JSpin
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
C
C                 Now we have a matrix of (approximate) 4c2e integrals
C                 for this particular IJ.
C
C                 The E(2) energy is a dot product of this matrix
C                 with a linear combination of the matrix and its
C                 transpose. Make the necessary workspace and initialize
C                 it with a copy of A.
C
                  If ( .NOT. GA_Create( MT_Dbl, NI*NVir(ISpin),
     $               NJ*NVir(JSpin), 'Temporary matrix big B', MinChunk,
     $               MinChunk, G_bigB) ) then
                     Write (LuErr, *) ' Unable to allocate big B'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
C
C                 ********************************************
C                 * Form the B matrix, which is related to A *
C                 ********************************************
C                 Term          A                  B
C                 RHF        (IA|JB)     2 (IA|JB) - (IB|JA)
C                 UHF AA     (IA|JB)       (IA|JB) - (IB|JA)
C                 UHF Ab     (IA|jb)       (IA|jb) 
C                 UHF bb     (ia|jb)       (ia|jb) - (ib|ja)
C
C                 Introducing SCS scaling, the above table transforms into
C
C                 Term          A                           B
C                 RHF        (IA|JB)      (fss+fos) (IA|JB) - fss (IB|JA)
C                 UHF AA     (IA|JB)         fss    (IA|JB) - fss (IB|JA)
C                 UHF Ab     (IA|jb)         fos    (IA|jb) 
C                 UHF bb     (ia|jb)         fss    (ia|jb) - fss (ib|ja)
C
C                 Where fss and fos are the same-spin and opposite-spin scaling
C                 factors, respectively. In the conventional case we have 
C                 fss = fos = 1.0, and we recover the previous expressions.
C
C                 So in every case except UHF Ab, we can initialize
C                 B with A transpose, but on a block-wise basis.
C
                  Call GA_Sync
                  Call PStat_On( PS_form_ab)
                  If ( PermSym ) then
                     Call GA_Zero(G_bigB)
C
                     Do I = ILo, IHi
C
C                       Use permutational symmetry if possible
C                       (already know that ISpin .eq. JSpin)
C                       
                        If ( ILo .eq. JLo .AND. IHi .eq. JHi) then
                           JTop2 = I
                        Else
                           JTop2 = JHi
                        EndIf
C
                        Do J = JLo, JTop2
C
C                          Copy from an IJ block of A to the 
C                          coresponding IJ block of B, but with
C                          transposition.  This routine is data parallel
C                          in that nodes which don't have any part of
C                          the patch to be copies will just fall through
C                          and immediately go on to the next loop
C                          iteration.
C
#ifdef BAD_GACCESS
                     Call rimp2_copypp('T',
#else
                           Call GA_Copy_Patch('T',
#endif
     $                        G_bigA, (I-ILo)*NVir(ISpin)+1,
     $                                (I-ILo+1)*NVir(ISpin),
     $                                (J-JLo)*NVir(JSpin)+1,
     $                                (J-JLo+1)*NVir(JSpin),
     $                        G_bigB, (I-ILo)*NVir(ISpin)+1,
     $                                (I-ILo+1)*NVir(ISpin),
     $                                (J-JLo)*NVir(JSpin)+1,
     $                                (J-JLo+1)*NVir(JSpin) )
                        EndDo
                     EndDo
C
C                    Now we must make sure everyone has a consistent
C                    copy of B before moving on.
C
                     Call GA_Sync
C
C                    Scale B and add some multiple of A to it,
C                    according to the SCS scaling above
C
                     Fact = fss + fos * Dble( 2 - TopSpin )
                     Call GA_DAdd( Fact, G_bigA, -fss, G_bigB,
     $                  G_bigB)
                  Else
C
C                    For the mimxed-spin case, its just a straight copy
C
                     Call GA_Copy( G_bigA, G_bigB)
                     if ( fos .ne. 1.0d0 ) call GA_Dscal( G_bigB, fos )
                  EndIf
                  Call PStat_Off( PS_form_ab)
C
                  If ( Util_Print('progress reports', Print_Default))
     $               then
                     String1 = ' ' ! Note all machines initialize mem.
                     Write (String1, 9220) ILo, IHi, JLo, JHi
                     Call Progress_Report(LuOut,
     $                  String1(:LNBlnk(String1)))
 9220                Format('RIMP2_V_E2: Copy to B done (',
     $                  I5, ':', I5, ',', I5, ':', I5, ')')
                  EndIf
C
                  If ( Util_Print('4c ints b', Print_Debug) )
     $               Call GA_Print( G_bigB )
C
C                 ***************************
C                 * Apply denominators to B *
C                 ***************************
C                 Each processor does this on the panel of B it owns.
C
                  Call GA_Sync
                  Call PStat_On( PS_denom)
                  Call RIMP2_Denom( G_bigB, ILo, JLo, NVir(ISpin),
     $               NVir(JSpin), Eig(1, ISpin), Eig(1, JSpin),
     $               Eig(1+NOcc(ISpin), ISpin),
     $               Eig(1+NOcc(JSpin), JSpin) )
                  Call PStat_Off( PS_denom)
C
                  If ( Util_Print('progress reports', Print_Default))
     $               then
                     String1 = ' ' ! Note all machines initialize mem.
                     Write (String1, 9230) ILo, IHi, JLo, JHi
                     Call Progress_Report(LuOut,
     $                  String1(:LNBlnk(String1)))
 9230                Format('RIMP2_V_E2: Denominators done (',
     $                  I5, ':', I5, ',', I5, ':', I5, ')')
                  EndIf
C
                  If ( Util_Print('amplitudes', Print_Debug) )
     $               Call GA_Print( G_bigB )
C
C                 ***********************
C                 * Form IJ pair energy *
C                 ***********************
C                 Another data parallel routine.
C
                  Call GA_Sync
                  Call PStat_On(PS_PrEn)
C
                  Do I = ILo, IHi
C
C                    Use permutational symmetry if possible
C                       
                     If ( (ILo .eq. JLo .AND. IHi .eq. JHi) 
     $                  .AND. PermSym) then
                        JTop2 = I
                     Else
                        JTop2 = JHi
                     EndIf
C
                     Do J = JLo, JTop2
C
#ifdef BAD_GACCESS
                     paire=rimp2_ddotp(
#else
                        PairE = GA_DDot_Patch(
#endif
     $                     G_bigA, 'N', (I-ILo)*NVir(ISpin)+1,
     $                                  (I-ILo+1)*NVir(ISpin),
     $                                  (J-JLo)*NVir(JSpin)+1,
     $                                  (J-JLo+1)*NVir(JSpin),
     $                     G_bigB, 'N', (I-ILo)*NVir(ISpin)+1,
     $                                  (I-ILo+1)*NVir(ISpin),
     $                                  (J-JLo)*NVir(JSpin)+1,
     $                                  (J-JLo+1)*NVir(JSpin) )
#ifndef BAD_GACCESS
                        PairE = PairE/ga_nnodes()
#endif
C
C                       If this turns out to be too much traffic, we
C                       might want to create a local array to contain
C                       only our IJs and accumulate into it.  But
C                       wait until perfomance says we need to.  For
C                       now, we'll just avoid sending zeros across
C                       the machine.
C
                        If ( PairE .ne. FP_Zero )
     $                     Call GA_Acc(
     $                        Int_MB( IH_PairE + ISpin + JSpin - 2),
     $                        I, I, J, J, PairE, 1, FP_One)
                     EndDo
                  EndDo
                  Call PStat_Off(PS_PrEn)
C
                  If ( Util_Print('partial pair energies', Print_Debug))
     $               Call GA_Print(
     $                        Int_MB( IH_PairE + ISpin + JSpin - 2) )
C
C                 Free up these matrices -- the space can be reused for
C                 the J integral matrix, and their sizes may change when
C                 the spin cases change.
C                 
                  If ( .NOT. GA_Destroy( G_bigB) ) then
                     Write (LuErr, *) ' Unable to free big B'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
                  If ( .NOT. GA_Destroy( G_bigA) ) then
                     Write (LuErr, *) ' Unable to free big A'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
C
                  If ( Util_Print('progress reports', Print_Default))
     $               then
                     String1 = ' ' ! Note all machines initialize mem.
                     Write (String1, 9200) ILo, IHi, JLo, JHi
                     Call Progress_Report(LuOut,
     $                  String1(:LNBlnk(String1)))
 9200                Format('RIMP2_V_E2: Completed block (',
     $                  I5, ':', I5, ',', I5, ':', I5, ')')
                  EndIf
C
               EndDo            ! J (blocks)
            EndDo               ! JSpin
         EndDo                  ! I (blocks)
C
C        This may change size with the next spin case
C
         If ( .NOT. GA_Destroy( G_bigI) ) then
            Write (LuErr, *) ' Unable to free big I integral block'
            RIMP2_V_E2 = .FALSE.
            Return
         EndIf
C        
      EndDo                     ! ISpin
C
C     ************************
C     * Close & delete files *
C     ************************
C
      Status = .TRUE.
      Do ISpin = 1, TopSpin
         Status = Status .AND. Close_3cI_File( D_Int(ISpin), 'DELETE',
     ,        actualname)
      EndDo
C
      If ( .NOT. Status ) then
         Write (LuErr, *) 'RIMP2_V_E2:  Problems closing ',
     $      'integral file(s)!'
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf
C
C     ***********************
C     * Print pair energies *
C     ***********************
C
      If ( Util_Print('pair energies', Print_Debug) ) then
         If ( Proc .eq. 0 ) then
            If ( TopSpin .eq. 1) then
               Call Banner(LuOut,
     $            'Pair Correlation Energies (spin adapted)',
     $            '-', .TRUE., .TRUE., .FALSE.)
            Else
               Call Banner(LuOut, 'Pair Correlation Energies',
     $            '-', .TRUE., .TRUE., .FALSE.)
            EndIf
         EndIf
C
         Do ISpin = 1, TopSpin
            Do JSpin = ISpin, TopSpin
               Call GA_Print( Int_MB(IH_PairE+ISpin+JSpin-2) )
            EndDo
         EndDo
      EndIf
C
C     ***********************************
C     * Compute total energy from pairs *
C     ***********************************
C
      Do ISpin = 1, TopSpin
         Do JSpin = ISpin, TopSpin
            Contrib(ISpin, JSpin) =
     $         Sum_Pairs( Int_MB(IH_PairE+ISpin+JSpin-2) ,
     $         ISpin .eq. JSpin)
         EndDo
      EndDo
C
C     ***************************
C     * Finished with the pairs *
C     ***************************
C
      Do ISpin = 1, TopSpin
         Do JSpin = ISpin, TopSpin
            Status = GA_Destroy( Int_MB(IH_PairE+ISpin+JSpin-2) )
         EndDo
      EndDo
C
      Status = Status .AND. MA_Pop_Stack(HH_PairE)
C
C     *****************************************************
C     * For UHF case, pure spin components must be scaled *
C     *****************************************************
C
      If ( TopSpin .eq. 2)
     $   Call dscal(TopSpin, FP_Half, Contrib, TopSpin+1)
C
C     *************************************
C     * And the total from the components *
C     *************************************
C     This means the caller doesn't have to fiddle with anything
C     if only the total correlation energy is important.
C
      E2 = FP_Zero
      Do JSpin = 1, TopSpin
         Do ISpin = 1, TopSpin
            E2 = E2 + Contrib(ISpin, JSpin)
         EndDo
      EndDo
C
      Call GA_Sync
      Call PStat_Off(PS_Energy)
C
      RIMP2_V_E2 = .TRUE.
      Return
      End
#ifdef BAD_GACCESS
      double precision function rimp2_ddotp(
     $                     g_a, tra, ailo,aihi,ajlo,ajhi,
     $                     g_b, trb, bilo,bihi,bjlo,bjhi)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_a, ailo,aihi,ajlo,ajhi,
     $        g_b, bilo,bihi,bjlo,bjhi
      character *1 tra,trb
      double precision ddot
      external ddot
c
      integer n,lda,ldb
      integer k_a,l_a,k_b,l_b
      integer me,iloa,ihia,jloa,jhia
c
      if(tra.ne.'N'.or.trb.ne.'N') call errquit(
     $ 'rimp2ddot not ready ',0,0)
      lda=aihi-ailo+1
      me=ga_nodeid()
      rimp2_ddotp=0d0
c
c     check intersect
c
      call ga_distribution(g_a, me,iloa,ihia,jloa,jhia)
      if( ihia.lt. ailo.or. aihi.lt.iloa) return 
      if( jhia.lt.ajlo.or. ajhi.lt.jloa) return 

      ldb=bihi-bilo+1
      n=(ajhi-ajlo+1)*(aihi-ailo+1)
      if(.not.ma_push_get(MT_DBL,(aihi-ailo+1)*(ajhi-ajlo+1),
     $  'scratch buff', l_a, k_a)) call
     $  errquit('rimp2_ddot: pushget failed',0,0)
      if(.not.ma_push_get(MT_DBL,(bihi-bilo+1)*(bjhi-bjlo+1),
     $  'scratch buff', l_b, k_b)) call
     $  errquit('rimp2_ddot: pushget failed',0,0)

      call ga_get(g_a,ailo,aihi,ajlo,ajhi,dbl_mb(k_a),lda)
      call ga_get(g_b,bilo,bihi,bjlo,bjhi,dbl_mb(k_b),ldb)
      rimp2_ddotp=ddot(n,dbl_mb(k_a),1,dbl_mb(k_b),1)
      if(.not.ma_pop_stack(l_b)) call 
     $  errquit('riddot: popstack failed',0,0)
      if(.not.ma_pop_stack(l_a)) call 
     $  errquit('riddot: popstack failed',0,0)
c
      return
      end
      subroutine rimp2_copypp(trans, g_a, ailo, aihi, ajlo, ajhi,
     *     g_b, bilo, bihi, bjlo, bjhi)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      character *1 trans 
      integer g_a, ailo, aihi, ajlo, ajhi,
     *     g_b, bilo, bihi, bjlo, bjhi
c
      integer n,lda,ldb
      integer k_a,l_a,k_b,l_b
      integer me,iloa,ihia,jloa,jhia
c
      lda=aihi-ailo+1
      ldb=bihi-bilo+1
      n=(ajhi-ajlo+1)*(aihi-ailo+1)
c
c     check intersect
c
      me=ga_nodeid()
      call ga_distribution(g_a, me,iloa,ihia,jloa,jhia)
      if( ihia.lt. ailo.or. aihi.lt.iloa) return 
      if( jhia.lt.ajlo.or. ajhi.lt.jloa) return 
      if(.not.ma_push_get(MT_DBL,(aihi-ailo+1)*(ajhi-ajlo+1),
     $  'scratch buff', l_a, k_a)) call
     $  errquit('rimp2_ddot: pushget failed',0,0)
      if (trans.ne.'N') then
         if(.not.ma_push_get(MT_DBL,(bihi-bilo+1)*(bjhi-bjlo+1),
     $        'scratch buff', l_b, k_b)) call
     $        errquit('rimp2_ddot: pushget failed',0,0)
      endif

      call ga_get(g_a,ailo,aihi,ajlo,ajhi,dbl_mb(k_a),lda)
      if(trans.eq.'N') then
         call ga_put(g_b,bilo,bihi,bjlo,bjhi,dbl_mb(k_a),ldb)
      else
         call util_transpose(dbl_mb(k_a),lda,dbl_mb(k_b),ldb,
     @  aihi-ailo+1,ajhi-ajlo+1)
         call ga_put(g_b,bilo,bihi,bjlo,bjhi,dbl_mb(k_b),ldb)
      endif
      if (trans.ne.'N') then
         if(.not.ma_pop_stack(l_b)) call 
     $        errquit('ricopyp: popstack failed',0,0)
      endif
      if(.not.ma_pop_stack(l_a)) call 
     $  errquit('ricopyp: popstack failed',0,0)
      return
      end
#endif
