C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rimp2_v_e2 -- Evaluate MP2 energy using V-approximation 3 ctr ints.
C
C REVISION
C     $Id: rimp2_v_e2.F,v 1.2 1994-06-16 22:44:52 gg502 Exp $
C
C SYNOPSIS
      Logical Function RIMP2_V_E2(FNInt, TopSpin, NMO, NOcc, NVir,
     $   NFit, FitBas, Eig, Contrib, E2)
      Implicit NONE
C
      Character*(*) FNInt
      Integer TopSpin, NMO, NFit, FitBas
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
C
C EXTERNAL ROUTINES
      Logical Open_3ci_File, Get_3cInts, Close_3cI_File,
     $  GA_Create_One_By_Shells
      Double Precision DSum
      External Open_3cI_File, Get_3cInts, Close_3cI_File, DSum,
     $  GA_Create_One_By_Shells
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
      Integer I, J, ISpin, JSpin, G_I, G_J, JHi
      Integer A, B, Alo, Ahi, Blo, BHi, G_A, G_B, H_B, LDB, AB
      Integer Proc
      Double Precision DenIJ, Den, Fact, EPair
      Logical Status, PermSym
      Integer Junk
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
         Status = GA_Create_One_By_Shells( MT_Dbl, NVir(ISpin), FitBas,
     $      'Integral matrix 1', G_I)
C
         If ( .NOT. Status ) then
            Write (LuErr, *) ' Unable to allocate for I block'
            RIMP2_V_E2 = .FALSE.
            Return
         EndIf
C
C        Open integral file
C
         Status = Open_3ci_file(LUInt+ISpin-1, FNInt, ISpin,
     $      NVir(ISpin), NOcc(ISpin), G_I)
         If ( .NOT. Status ) then
            Write (LuErr, *) 'RIMP2_v_E2:  Problems opening ',
     $         'integral file for spin ', ISpin, '!'
            RIMP2_v_E2 = .FALSE.
            Return
         EndIf
C
         Do I = 1, NOcc(ISpin)
C
C           Get integrals (ia|s') for this I
C
            If ( .NOT. Get_3cInts(LUInt+ISpin-1, G_I, NVir(ISpin), I, I)
     $         ) then
               Write( LuErr, *) ' Unable to get integral block ',
     $            I, ISpin
               RIMP2_V_E2 = .FALSE.
               Return
            EndIf
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
C              Space for a block of integrals -- size varies
C              with spin case
C
               If ( .NOT. GA_Create( MT_Dbl, NVir(ISpin),
     $            NVir(JSpin), 'Temporary matrix A', MinChunk,
     $            MinChunk, G_A) ) then
                  Write (LuErr, *) ' Unable to allocate A'
                  RIMP2_V_E2 = .FALSE.
                  Return
               EndIf
               Call GA_Distribution( G_A, Proc, Alo, Ahi, Blo, Bhi)
c$$$               Write (LuOut, *) 'A distribution ', Proc, Alo, AHi,
c$$$     $            Blo, BHi
C
C              Open integral file if JSpin != ISpin
C              NOTE: This routine gets the distribution of the fitting
C              basis from the GA that's passed in.  We don't create the
C              real target array for these guys until a little later,
C              but we can use G_I here because the fitting basis will
C              not change with the spin case.
C
               If ( JSpin .ne. ISpin) then
                  Status = Open_3ci_file(LUInt+JSpin-1, FNInt, JSpin,
     $               NVir(JSpin), NOcc(JSpin), G_I)
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
                  JHi = I
                  PermSym = .TRUE.
               Else
                  JHi = NOcc(JSpin)
                  PermSym = .FALSE.
               EndIf
C
               Do J = 1, JHi
C
C                 Form half of the denominator
C
                  DenIJ = Eig(I, ISpin) + Eig(J, JSpin)
C
C                 Create the global array for the J block
C
                  Status = GA_Create_One_By_Shells( MT_Dbl, NVir(JSpin),
     $               FitBas, 'Integral matrix 2', G_J)
C
                  If ( .NOT. Status ) then
                     Write (LuErr, *) ' Unable to allocate integral ',
     $                  ' block ', J, JSpin
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
C
C                 Find the integrals (jb|s') for this J
C
                  If ( .NOT. Get_3cInts(LUInt+JSpin-1, G_J, NVir(JSpin),
     $               J, J) ) then
                     Write( LuErr, *) ' Unable to get integral block ',
     $                  J, JSpin
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
C
C                 ******************************************************
C                 * Form a block of four-center two-electron integrals *
C                 ******************************************************
C                 from these three-center things:
C                    (ia|jb) = (ia|s')(s'|jb)
C                 where the |s') = V^{-1/2} |s) has been done
C                 as part of the transformation.
C
                  Call GA_DGEMM('N', 'T', NVir(ISpin), NVir(JSpin),
     $               NFit, FP_One, G_I, G_J, FP_Zero, G_A)
c$$$C
c$$$                  Call GA_Print(G_I)
c$$$                  Call GA_Print(G_J)
c$$$                  Call GA_Print(G_A)
C
C                 The J integral block is no longer needed and can
C                 now be freed.  This alows us to reuse the space
C                 to hold the temporary matrix B for the other
C                 half of this iteration.
C
c$$$                  Write (LuOut, *) 'About to destroy J'
                  If ( .NOT. GA_Destroy( G_J) ) then
                     Write (LuErr, *) ' Unable to free integral block ',
     $                  J, JSpin
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
                  If ( .NOT. GA_Create( MT_Dbl, NVir(ISpin),
     $               NVir(JSpin), 'Temporary matrix B', MinChunk,
     $               MinChunk, G_B) ) then
                     Write (LuErr, *) ' Unable to allocate B'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
               Call GA_Distribution( G_B, Proc, Alo, Ahi, Blo, Bhi)
c$$$               Write (LuOut, *) 'B distribution ', Proc, Alo, AHi,
c$$$     $            Blo, BHi
C
                  Call GA_Copy( G_A, G_B)
C
C                 *************************************************
C                 * Make B into what is needed by pure spin cases *
C                 *************************************************
C                 For mixed spin UHF B = A, but for others we need
C                 to do some more manipulation.
C
                  If ( ISpin .eq. JSpin) then
C
C                    B <-- 0.5 ( A + A')
C
                     Call GA_Symmetrize( G_B)
c$$$C
c$$$                  Call GA_Print(G_A)
c$$$                  Call GA_Print(G_B)
C
C                    UHF: B <-- 2 A - 2 B = 2 A - (A + A') =   A - A'
C                    RHF: B <-- 3 A - 2 B = 3 A - (A + A') = 2 A - A'
C
                     Fact = Dble( 4 - TopSpin )
                     Call GA_DAdd( Fact, G_A, -FP_Two, G_B, G_B)
                  EndIf
c$$$C
c$$$                  Call GA_Print(G_A)
c$$$                  Call GA_Print(G_B)
C
C                 So now we have the following for a particular IJ:
C                 Term          A                  B
C                 RHF        (IA|JB)     2 (IA|JB) - (IB|JA)
C                 UHF AA     (IA|JB)       (IA|JB) - (IB|JA)
C                 UHF Ab     (IA|jb)       (IA|kb) 
C                 UHF bb     (ia|jb)       (ia|jb) - (ib|ja)
C
C                 ***************************
C                 * Apply denominators to B *
C                 ***************************
C                 Each processor does this on the panel of B it owns.
C
                  Call GA_Distribution( G_B, Proc, Alo, Ahi, Blo, Bhi)
C
C                 GA_Access doesn't like it if we request a bogus
C                 patch (as might be returned if this node has no data
C                 on it).
C
                  If ( (Ahi-Alo+1) * (BHi-BLo+1) .ne. 0) then
                     Call GA_Access( G_B, Alo, Ahi, Blo, Bhi, H_B, LDB)
c$$$C
c$$$                  Call GA_Print(G_B)
c$$$                  Write (LuOut, *) 'DenIJ = ', DenIJ, I, J
C
                     AB = -1
                     Do B = Blo, Bhi
                        Do A = Alo, Ahi
                           AB = AB + 1
C
                           Den = DenIJ - Eig( NOcc(ISpin) + A, ISpin)
     $                        - Eig( NOcc(JSpin) + B, JSpin)
C
c$$$                        Write (LuOut, *) 'Den ', Den, I, J, A, B
                           Dbl_MB( H_B + AB )
     $                        =  Dbl_MB( H_B + AB ) / Den
C
                        EndDo
                        AB = AB + LDB - (Ahi - Alo + 1)
                     EndDo
C
                     Call GA_Release_Update( G_B, Alo, Ahi, Blo, Bhi)
                  EndIf
c$$$C
c$$$                  Call GA_Print(G_A)
c$$$                  Call GA_Print(G_B)
C
C                 ***********************
C                 * Form IJ pair energy *
C                 ***********************
C                 Contract the integrals in A with the denominatorified 
C                 integrals in B. Put it into the matrix of 
C                 contributions by spin case because that breakdown is 
C                 sometimes interesting to people. 
C
                  EPair = GA_DDot( G_A, G_B)
                  If ( PermSym .AND. I .ne. J) then
                     Contrib( ISpin, JSpin)
     $                  = Contrib( ISpin, JSpin) + FP_Two * EPair
                  Else
                     Contrib( ISpin, JSpin)
     $                  = Contrib( ISpin, JSpin) + EPair
                  EndIf
C
                  If ( Proc .eq. 0 )
     $               Write (LuOut, 9000) I, J, ISpin, JSpin, EPair
 9000             Format( 1X, 'Pair (', I3, ',', I3, ') spin (',
     $               I2, ',', I2, ') ', F30.15)
C
C                 Free up the B matrix -- the space can be used by
C                 a J 3ci block before we need another B.
C
c$$$                  Write ( LuOut,*) 'About to destroy B'
                  If ( .NOT. GA_Destroy( G_B) ) then
                     Write (LuErr, *) ' Unable to free B'
                     RIMP2_V_E2 = .FALSE.
                     Return
                  EndIf
C
               EndDo   ! J
C
C              This may change size for the next spin case
C
c$$$                  Write ( LuOut,*) 'About to destroy A'
               If ( .NOT. GA_Destroy( G_A) ) then
                  Write (LuErr, *) ' Unable to free A'
                  RIMP2_V_E2 = .FALSE.
                  Return
               EndIf
C
            EndDo      ! JSpin
         EndDo         ! I
C
C        This may change size with the next spin case
C
c$$$                  Write ( LuOut,*) 'About to destroy I'
         If ( .NOT. GA_Destroy( G_I) ) then
            Write (LuErr, *) ' Unable to free I integral block'
            RIMP2_V_E2 = .FALSE.
            Return
         EndIf
C
      EndDo            ! ISpin
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
      RIMP2_V_E2 = .TRUE.
      Return
      End
