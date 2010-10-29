C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
CNAME
C     RIMP2_Mk_T -- Form T amplitudes for IJ batch in memory
C     
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function RIMP2_Mk_T(ILo, IHi, NVirI, JLo, JHi, NVirJ,
     $   IJSym, IJTri, G_Int, G_T, OccEigI, OccEigJ,
     $   VirEigI, VirEigJ, SpinLabel)
      Implicit NONE
C
      Integer ILo      ![in]  Lowest occ index in batch (row)
      Integer IHi      ![in]  Highest occ index in batch (row)
      Integer NVirI    ![in]  Size of blocks (row)
      Integer JLo      ![in]  Lowest occ index in batch (col)
      Integer JHi      ![in]  Highest occ index in batch (col)
      Integer NVirJ    ![in]  Size of blocks (col)
      Logical IJSym    ![in]  Use I <--> J permutational symmetry if possible
      Logical IJTri    ![in]  G_I/G_T hold only triangular I > J
      Integer G_Int    ![in]  GA handle for integral batch
      Integer G_T      ![out] GA handle for amplitude batch (created)
      Double Precision OccEigI(*)     ![in] Eigenvalues for index I (occ) 
      Double Precision OccEigJ(*)     ![in] Eigenvalues for index J (occ) 
      Double Precision VirEigI(NVirI) ![in] Eigenvalues for index I (vir) 
      Double Precision VirEigJ(NVirJ) ![in] Eigenvalues for index I (vir) 
      Character*(*) SpinLabel ![in] text label for current spin case
C
C DESCRIPTION
C     Given an appropriate batch of (antisymmetrized) integrals, creates
C     a matching set of amplitudes by copying integrals into a new GA
C     (same distribution as integrals) and applying denominators.
C
C     Integrals & amplitudes are expected to be organized as blocks of
C     all AB for a particular IJ.  A batch will contain blocks for
C     multiple Is (row dimension) and Js (column dimension) laid out
C     in rectangular fashion.  IJTri indicates that the minimum (I > J)
C     triangle is provided instead of the full rectangular matrix.
C
C     Where indicated (by IJSym and the fact that the requested ILo:IHi
C     and JLo:JHi ranges match up), only the I > J triangle of the
C     batch will be "denominatored" (the whole thing will be copied).
C
C PRINT CONTROLS
C     'amplitudes' will print G_T at the end of this routine
C
C ROUTINES
C     none
C
C INCLUDED FILES
#include "global.fh"
#include "util.fh"
#include "inp.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Character*255 String1
C
C     Prepare a GA for the amplitudes which matches the integrals
C     in shape and distribution.  That way we can capture anything
C     special the caller has done.
C     
      String1 = 'T2 amplitudes ' // SpinLabel
      RIMP2_Mk_T = GA_Duplicate(G_Int, G_T,
     $   String1(:Inp_StrLen(String1)) )
      If ( .NOT. RIMP2_Mk_T ) Return  ! Perhaps caller can handle failure!
C        
      Call GA_Copy( G_Int, G_T)
C     
C     **********************
C     * Apply denominators *
C     **********************
C     Each processor does this on the panel of T it owns.
C     
      Call GA_Sync
      Call RIMP2_Denom2( IJSym, IJTri, G_T, ILo, IHi, JLo, JHi,
     $   NVirI, NVirJ, OccEigI, OccEigJ, VirEigI, VirEigJ)
      Call GA_Sync
C     
C     Print results if requested
C
      If ( Util_Print('amplitudes', Print_Debug) ) Call GA_Print( G_T )
C
      Return
      End
