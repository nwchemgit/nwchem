C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_v_e2 -- Evaluate MP2 energy using V-approximation 3 ctr ints.
C
C REVISION
C     $Id: rimp2_v_e2.F,v 1.5 1994-08-25 20:28:05 gg502 Exp $
C
C SYNOPSIS
      Logical Function RIMP2_V_E2(FNInt, TopSpin, NMO, NOcc, NVir,
     $   NFit, FitBas, Eig, Contrib, E2, IBlock, JBlock)
      Implicit NONE
C
      Character*(*) FNInt
      Integer TopSpin, NMO, NFit, FitBas, IBlock, JBlock
      Integer NOcc(TopSpin), NVir(TopSpin)
      Double precision Eig(NMO, TopSpin), Contrib(TopSpin, TopSpin), E2
C
C ARGUMENTS
C     FNInt    File name for MO basis 3-center integrals [IN]
C     TopSpin  Number of spins in reference (RHF=1, UHF=1) [IN]
C     NMO      Number of MOs [IN]
C     NOcc     Number of occupied orbitals for each spin [IN]
C     NVir     Number of virtual orbitals for each spin [IN]
C     NFit     Number of functions in fitting basis [IN]
C     FitBas   Handle for fitting basis [IN]
C     Eig      Eigenvalues of reference Fock matrix [IN]
C     Contrib  E(2) broken down by spin case. (The ba, or 2,1 element
C              is always zero) [OUT]
C     E2       The second order MBPT energy [OUT]
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
#include "mafdecls.h"
#include "numerical_constants.fh"
#include "stdio.fh"
#include "rimp2_timers.fh"
C
C EXTERNAL ROUTINES
      Logical Open_3ci_File, Get_3cInts, Close_3cI_File,
     $  GA_Create_One_By_Shells
      Double Precision DSum
      External Open_3cI_File, Get_3cInts, Close_3cI_File, DSum,
     $  GA_Create_One_By_Shells
c$$$C
c$$$      External GA_DASum_Patch, GA_DAMax_Patch
c$$$      Double precision GA_DASum_Patch, GA_DAMax_Patch
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C PARAMETERS
C
      Integer MinChunk             ! Minimum panel size that is reasonable
      Parameter ( MinChunk = 64)   ! for global distribution of data
C
      Integer LUInt
      Parameter (LuInt = 10)
C
C LOCAL VARIBLES
      Integer I, J, ISpin, JSpin, G_J, JTop, JTop2
      Integer A, B, Alo, Ahi, Blo, BHi, G_A, G_B, H_B, LDB, AB
      Integer Proc
      Double Precision DenIJ, Den, Fact, EPair
      Logical Status, PermSym
      Integer Junk
      Integer G_bigI, ILo, IHi, NI, G_bigJ, JLo, JHi, NJ, G_bigA, G_BigB
      Integer IALo, IAHi, JBLo, JBHi, MyILo, MyIHi, MyJLo, MyJHi
      Integer Offset, BBot, BTop, ABot, ATop
c$$$      Double precision ChkSum, MxElem
C
      Call GA_Sync
      Call PStat_On(PS_Energy)
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
      If (NMO .lt. 1) then
         Write (LuErr, *) ' Bad value of NMO ', NMO
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
      If (NOcc(1) + NVir(1) .gt. NMO
     $   .OR. NOcc(TopSpin) + NVir(TopSpin) .gt. NMO ) then
         Write (LuErr, *) ' Bad value of NOcc+NVir ', NOcc, '+', NVir
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf
      If (NFit .lt. 1) then
         Write (LuErr, *) ' Bad value of NFit ', NFit
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf

c$$$      Write (LuOut, *) 'In energy!'
C
C     Get my processor number
C
      Proc = GA_NodeID()
C
C     Clear out the energy contribution matrix
C
      Call DGEZER(TopSpin, TopSpin, Contrib, TopSpin)
c$$$C
c$$$      Call DGEWr(NMO, TopSpin, Eig, NMO, 0, 0, LuOut,
c$$$     $   'SCF eigenvalues', 80, 'E15.6', Junk)
C
C     Label the pair energies we'll print within the loops
C
      If ( Proc .eq. 0 ) then
         If ( TopSpin .eq. 1) then
            Call Banner(LuOut,
     $         'Pair Correlation Energies (spin adapted)',
     $         '-', .TRUE., .TRUE., .FALSE.)
         Else
            Call Banner(LuOut, 'Pair Correlation Energies',
     $         '-', .TRUE., .TRUE., .FALSE.)
         EndIf
      EndIf
C
C     *****************************************************************
C     * Loop over occupied orbitals of both spins for the first index *
C     *****************************************************************
C
      Do ISpin = 1, TopSpin
C
C        Create the global array for the I block
C
         Call PStat_On( PS_Create)
         Status = GA_Create_One_By_Shells( MT_Dbl,
     $      Min( NOcc(ISpin), IBlock) * NVir(ISpin), FitBas,
     $      'Integral matrix 1', G_bigI)
C
         If ( .NOT. Status ) then
            Write (LuErr, *) ' Unable to allocate for big I block'
            RIMP2_V_E2 = .FALSE.
            Return
         EndIf
         Call PStat_Off( PS_Create)
C
C        Open integral file
C
         Status = Open_3ci_file(LUInt+ISpin-1, FNInt, ISpin,
     $      NVir(ISpin), NOcc(ISpin), G_bigI)
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
            If ( .NOT. Get_3cInts(LUInt+ISpin-1, G_bigI, NVir(ISpin),
     $         ILo, IHi) ) then
               Write( LuErr, *) ' Unable to get integral block ',
     $            I, ISpin
               RIMP2_V_E2 = .FALSE.
               Return
            EndIf
            Call PStat_Off( PS_Rd_3ci)
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
                  Status = Open_3ci_file(LUInt+JSpin-1, FNInt, JSpin,
     $               NVir(JSpin), NOcc(JSpin), G_bigI)
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
C                 Create the global array for the J block
C
                  Call PStat_On( PS_Create)
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
                  If ( .NOT. GA_Create( MT_Dbl, NI*NVir(ISpin),
     $               NJ*NVir(JSpin), 'Temporary matrix big A', MinChunk,
     $               MinChunk, G_bigA) ) then
                     Write (LuErr, *) ' Unable to allocate big A'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
                  Call PStat_Off( PS_Create)
C
C                 Find the integrals (jb|s') for this J
C
                  Call PStat_On( PS_Rd_3ci)
                  If ( .NOT. Get_3cInts(LUInt+JSpin-1, G_bigJ,
     $               NVir(JSpin), JLo, JHi) ) then
                     Write( LuErr, *) ' Unable to get integral block ',
     $                  JLo, JSpin
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
                  Call PStat_Off( PS_Rd_3ci)
C
C                 ******************************************************
C                 * Form a block of four-center two-electron integrals *
C                 ******************************************************
C                 from these three-center things:
C                    (ia|jb) = (ia|s')(s'|jb)
C                 where the |s') = V^{-1/2} |s) has been done
C                 as part of the transformation.
C
                  Call PStat_On( PS_form_4ci)
                  Call GA_DGEMM('N', 'T',
     $               NI*NVir(ISpin), NJ* NVir(JSpin), NFit,
     $               FP_One, G_bigI, G_bigJ, FP_Zero, G_bigA)
                  Call PStat_Off( PS_form_4ci)
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
                  Call PStat_On( PS_Create)
                  If ( .NOT. GA_Create( MT_Dbl, NVir(ISpin),
     $               NVir(JSpin), 'Temporary matrix A', MinChunk,
     $               MinChunk, G_A) ) then
                     Write (LuErr, *) ' Unable to allocate A'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
                  If ( .NOT. GA_Create( MT_Dbl, NVir(ISpin),
     $               NVir(JSpin), 'Temporary matrix B', MinChunk,
     $               MinChunk, G_B) ) then
                     Write (LuErr, *) ' Unable to allocate B'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
                  If ( .NOT. GA_Create( MT_Dbl, NI*NVir(ISpin),
     $               NJ*NVir(JSpin), 'Temporary matrix big B', MinChunk,
     $               MinChunk, G_bigB) ) then
                     Write (LuErr, *) ' Unable to allocate big B'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
                  Call PStat_Off( PS_Create)
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
C                 So in every case except UHF Ab, we can initialize
C                 B with A transpose, but on a block-wise basis.
C
                  Call PStat_On( PS_form_ab)
                  If ( ISpin .eq. JSpin) then
                     Call GA_Zero(G_bigB)
C
C                    NOTE: the operations in here are SIMD, so are
C                    not as parallel as one might like. This could
C                    be optimized away by creating a specialized
C                    routine to do the whole thing, but we'll wait
C                    and see if there is a problem.
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
C                          transposition.
C
                           Call GA_Copy_Patch('T',
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
C                    In all of these cases, we also have to scale B
C                    by -1 and add some multiple of A to it,
C                    2 for RHF, 1 for UHF.
C
                     Fact = Dble( 3 - TopSpin)
                     Call GA_DAdd( Fact, G_bigA, -FP_One, G_bigB,
     $                  G_bigB)
                  Else
C
C                    For the mimxed-spin case, its just a straight copy
C
                     Call GA_Copy( G_bigA, G_bigB)
                  EndIf
                  Call PStat_Off( PS_form_ab)
C
C                 ***************************
C                 * Apply denominators to B *
C                 ***************************
C                 Each processor does this on the panel of B it owns.
C
                  Call PStat_On( PS_denom)
                  Call GA_Distribution( G_bigB, Proc, IAlo, IAhi,
     $               JBlo, JBhi)
C
C                 GA_Access doesn't like it if we request a bogus
C                 patch (as might be returned if this node has no data
C                 on it).
C
                  If ( (IAhi-IAlo+1) * (JBHi-JBLo+1) .ne. 0) then
C
C                    Figure out what these correspond to in more 
C                    basic terms
C
                     ALo = Mod( IALo, NVir(ISpin) )
                     If ( ALo .eq. 0) ALo = NVir(ISpin)
                     AHi = Mod( IAHi, NVir(ISpin) )
                     If ( AHi .eq. 0) AHi = NVir(ISpin)
C
                     BLo = Mod( JBLo, NVir(JSpin) )
                     If ( BLo .eq. 0) BLo = NVir(JSpin)
                     BHi = Mod( JBHi, NVir(JSpin) )
                     If ( BHi .eq. 0) BHi = NVir(JSpin)
C
                     MyILo = Int( (IALo-1) / NVir(ISpin) ) + 1
                     MyIHi = Int( (IAHi-1) / NVir(ISpin) ) + 1
C
                     MyJLo = Int( (JBLo-1) / NVir(JSpin) ) + 1
                     MyJHi = Int( (JBHi-1) / NVir(JSpin) ) + 1
C
                     Call GA_Access( G_bigB, IAlo, IAhi, JBlo, JBhi,
     $                        H_B, LDB)
C
                     Offset = -1
                     Do J = MyJLo, MyJHi
C
C                       First and last J might not cover the whole range
C
                        BBot = 1
                        BTop = NVir(JSpin)
                        If ( J .eq. MyJLo) then
                           BBot = BLo
                        ElseIf (J .eq. MyJHi) then
                           BTop = BHi
                        EndIf
C
                        Do B = BBot, BTop
C
                           Do I = MyILo, MyIHi
C                          
C                             Form half of the denominator
C
                              DenIJ = Eig(I+ILo-1, ISpin)
     $                           + Eig(J+JLo-1, JSpin)
C
C                             As for B, A doesn't always cover the whole
C                             range
C
                              ABot = 1
                              ATop = NVir(ISpin)
                              If ( I .eq. MyILo) then
                                 ABot = ALo
                              ElseIf (I .eq. MyIHi) then
                                 ATop = AHi
                              EndIf
C
                              Do A = ABot, ATop
                                 Offset = Offset + 1
C
                                 Den = DenIJ
     $                              - Eig( NOcc(ISpin) + A, ISpin)
     $                              - Eig( NOcc(JSpin) + B, JSpin)
C
                                 Dbl_MB( H_B + Offset )
     $                              =  Dbl_MB( H_B + Offset ) / Den
C
                              EndDo ! A
                           EndDo ! I
                           Offset = Offset + LDB - (IAhi - IAlo + 1)
C
                        EndDo ! B
                     EndDo ! J
C
                     Call GA_Release_Update( G_bigB, IAlo, IAhi,
     $                  JBlo, JBhi)
                  EndIf
                  Call PStat_Off( PS_denom)
C
C                 ***********************
C                 * Form IJ pair energy *
C                 ***********************
C                 Contract the integrals in A with the denominatorified 
C                 integrals in B. Put it into the matrix of 
C                 contributions by spin case because that breakdown is 
C                 sometimes interesting to people. 
C
C                 NOTE: The dot products in this section are SIMD,
C                 which limits the available parallelism.  It might
C                 be necessary to write a special routine to do it.
C
                  Call PStat_On(PS_PrEn)
                  Do I = ILo, IHi
C
C                    Use permutational symmetry if possible
C                       
                     If ( ISpin .eq. JSpin .AND. ILo .eq. JLo
     $                  .AND. IHi .eq. JHi) then
                        JTop2 = I
                     Else
                        JTop2 = JHi
                     EndIf
C                    
                     Do J = JLo, JTop2
                           EPair = GA_DDot_Patch(
     $                        G_bigA, 'N', (I-ILo)*NVir(ISpin)+1,
     $                                     (I-ILo+1)*NVir(ISpin),
     $                                     (J-JLo)*NVir(JSpin)+1,
     $                                     (J-JLo+1)*NVir(JSpin),
     $                        G_bigB, 'N', (I-ILo)*NVir(ISpin)+1,
     $                                     (I-ILo+1)*NVir(ISpin),
     $                                     (J-JLo)*NVir(JSpin)+1,
     $                                     (J-JLo+1)*NVir(JSpin) )
C
                        If ( PermSym .AND. I .ne. J) then
                           Contrib( ISpin, JSpin)
     $                        = Contrib( ISpin, JSpin) + FP_Two * EPair
                        Else
                           Contrib( ISpin, JSpin)
     $                        = Contrib( ISpin, JSpin) + EPair
                        EndIf
C
                        If ( Proc .eq. 0 )
     $                     Write (LuOut, 9000) I, J, ISpin, JSpin, EPair
 9000                   Format( 1X, 'Pair (', I3, ',', I3, ') spin (',
     $                     I2, ',', I2, ') ', F30.15)
C                       
                     EndDo      ! J
                  EndDo         ! I
                  Call PStat_Off(PS_PrEn)
C
C                 Free up these matrices -- the space can be reused for
C                 the J integral matrix, and their sizes may change when
C                 the spin cases change.
C                 
                  If ( .NOT. GA_Destroy( G_B) ) then
                     Write (LuErr, *) ' Unable to free B'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
                  If ( .NOT. GA_Destroy( G_A) ) then
                     Write (LuErr, *) ' Unable to free A'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
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
C     ************************
C     * Close & delete files *
C     ************************
C
      Status = .TRUE.
      Do ISpin = 1, TopSpin
         Status = Status .AND. Close_3cI_File( LUInt+ISpin-1, 'DELETE')
      EndDo
C
      If ( .NOT. Status ) then
         Write (LuErr, *) 'RIMP2_V_E2:  Problems closing ',
     $      'integral file(s)!'
         RIMP2_V_E2 = .FALSE.
         Return
      EndIf
C
      Call GA_Sync
      Call PStat_Off(PS_Energy)
C
      RIMP2_V_E2 = .TRUE.
      Return
      End
