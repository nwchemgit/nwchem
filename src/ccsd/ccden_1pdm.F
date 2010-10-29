C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCDen_1PDM -- Form coupled cluster 1-particle density matrix
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine CCDen_1PDM(NOcc, NVir, T1, LDT1, g_T2,
     $   Z1, LDZ1, g_Z2,
     $   Omega1, LDO1, Omega2, LDO2, Omega3, LDO3, Omega4, LDO4,
     $   g_1PDM)
      Implicit NONE
#include "errquit.fh"
      Integer NOcc, NVir, LDT1, LDZ1, LDO1, LDO2, LDO3, LDO4
      Integer g_T2, g_Z2, g_1PDM
      Double Precision T1(LDT1, NVir), Z1(LDT1, NVir)
      Double Precision Omega1(LDO1, NVir), Omega2(LDO2, NVir),
     $   Omega3(LDO3, NVir), Omega4(LDO4, NOcc)
C
C ARGUMENTS
C DESCRIPTION
C
C     The output of this routine is a global array, for convenience 
C     in interfacing with existing NWChem machinery, the inputs are
C     all local, replicated data.  Consequently, the simplest
C     implementation has process 0 doing all the work.
C
C REFERENCES
C     Alistair P. Rendell and Timothy J. Lee, J. Chem. Phys. 94, 
C     6219--6228 (1991). Especially Eq. 43--45
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
C
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer H_Scr1, I_Scr1
C
C     *********************
C     * Initialize Output *
C     *********************
C
      Call GA_Zero(G_1PDM)
C
C     All inputs are replicated, so only process0 needs to do the work.
C
      If ( GA_NodeID() .eq. 0) then
C
C        ****************************************************
C        * Occupied-occupied block                          *
C        * Q(i,j) = -1/2 P(i,j) { [tZ](i,j) + Omega4(i,j) } *
C        *           + 2Delta(i,j)                          *
C        * APR & TJL Eq. 43                                 *
C        ****************************************************
C
         If (.NOT. MA_Push_Get(MT_Dbl, NOcc*NOcc,
     $      'CCDen Scratch 1', H_Scr1, I_Scr1)) Call ErrQuit(
     $      'CCDen_1PDM: unable to allocate scratch', NOcc*NOcc,
     &       MA_ERR)
C
C        Initialize the matrix with 2 on the diagonal, per the last term.
C
         Call dlaset('Full', NOcc, NOcc, 0.0d0, 2.0d0,
     $      Dbl_MB(I_Scr1), NOcc)
C
C        [tZ](i,j) = sum(a) t(i,a) Z(j,a) (APR & TJL Eq. 56)
C        Minus in DGEMM is due to minus in front of whole expression
C
         Call dgemm('N', 'T', NOcc, NOcc, NVir,
     $      -1.0d0, T1, LDT1, Z1, LDZ1, 1.0d0, Dbl_MB(I_Scr1), NOcc)
C
C        Now add in Omega4
C
         Call GEAXPY(NOcc, NOcc, -1.0d0, Omega4, NOcc,
     $      Dbl_MB(I_Scr1), NOcc)
C
C        The 1/2 P(i,j) permutation operation amounts to symmetrization
C        of the block.  We need to do the same thing to the other blocks
C        so we can defer it until later.
C
C        So we've got everything we need to put into the oo block...
C
         Call GA_Acc(g_1PDM, 1, NOcc, 1, NOcc, Dbl_MB(I_Scr1), NOcc,
     $      1.0d0)
C
         If ( .NOT. MA_Pop_Stack(H_Scr1) ) Call ErrQuit(
     $      'CCDen_1PDM: unable to free scratch', 0, MA_ERR)
C
C        ***************************************************
C        * Virtual-virtual block                           *
C        * Q(a,b) = 1/2 P(a,b) { [tZ](a,b) + Omega3(a,b) } *
C        * APR & TJL Eq. 44                                *
C        ***************************************************
C
         If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $      'CCDen Scratch 1', H_Scr1, I_Scr1)) Call ErrQuit(
     $      'CCDen_1PDM: unable to allocate scratch', NVir*NVir,
     &       MA_ERR)
C
C        [tZ](a,b) = sum(i) t(i,a) Z(i,b) (APR & TJL Eq. 57)
C
         Call dgemm('T', 'N', NVir, NVir, NOcc,
     $      1.0d0, T1, LDT1, Z1, LDZ1, 0.0d0, Dbl_MB(I_Scr1), NVir)
C
C        Now add in Omega3
C
         Call GEAXPY(NVir, NVir, 1.0d0, Omega3, NVir,
     $      Dbl_MB(I_Scr1), NVir)
C
C        The 1/2 P(i,j) permutation operation amounts to symmetrization
C        of the block.  We need to do the same thing to the other blocks
C        so we can defer it until later.
C
C        So we've got everything we need to put into the vv block...
C
         Call GA_Acc(g_1PDM, Nocc+1, NOcc+NVir, NOcc+1, NOcc+NVir,
     $      Omega3, NVir, 1.0d0)
C
         If ( .NOT. MA_Pop_Stack(H_Scr1) ) Call ErrQuit(
     $      'CCDen_1PDM: unable to free scratch', 0,
     &       MA_ERR)
C
C        ***************************************************************
C        * Occupied-virtual block                                      *
C        * Q(i,a) = 1/2{ 2t(i,a) + Z(i,a) + 2Omega1(i,a) - Omega2(i,a) *
C        *          - sum(b) Omega3(a,b) t(i,b)                        * 
C        *          - sum(j) Omega4(i,j) t(j,a)                        * 
C        * APR & TJL Eq. 45                                            * 
C        ***************************************************************
C
         If (.NOT. MA_Push_Get(MT_Dbl, NOcc*NVir,
     $      'CCDen Scratch 1', H_Scr1, I_Scr1)) Call ErrQuit(
     $      'CCDen_1PDM: unable to allocate scratch', NOcc*NVir,
     &       MA_ERR)
C
         Call DFill(NOcc*NVir, 0.0d0, Dbl_MB(I_Scr1), 1)
C
C        The first four terms are _really_ easy
C
         Call GEAXPY(NOcc, NVir, 2.0d0, T1, LDT1, Dbl_MB(I_Scr1), NOcc)
         Call GEAXPY(NOcc, NVir, 1.0d0, Z1, LDZ1, Dbl_MB(I_Scr1), NOcc)
         Call GEAXPY(NOcc, NVir, 2.0d0, Omega1, LDO1,
     $      Dbl_MB(I_Scr1), NOcc)
         Call GEAXPY(NOcc, NVir, -1.0d0, Omega2, LDO2,
     $      Dbl_MB(I_Scr1), NOcc)
C
C        sum(b) t(i,b) Omega3(a,b)
C
         Call dgemm('N', 'T', NOcc, NVir, NVir,
     $      -1.0d0, T1, LDT1, Omega3, LDO3,
     $      1.0d0, Dbl_MB(I_Scr1), NOcc)
C
C        sum(j) Omega4(i,j) t(j,a)
C
         Call dgemm('N', 'N', NOcc, NVir, NOcc,
     $      -1.0d0, Omega4, LDO4, T1, LDT1,
     $      1.0d0, Dbl_MB(I_Scr1), NOcc)
C
C        Now, plop it into the 1PDM
C
         Call GA_Acc(g_1PDM, 1, NOcc, NOcc+1, NOcc+NVir,
     $      Dbl_MB(I_Scr1), NOcc, 1.0d0)
C
         If ( .NOT. MA_Pop_Stack(H_Scr1) ) Call ErrQuit(
     $      'CCDen_1PDM: unable to free scratch', 0, MA_ERR)
C
C
C
      EndIf ! End of process 0 work
      Call GA_Sync
C
C     To this point, we have ignored the 1/2P(i,j) (oo block), the
C     1/2P(a,b) (vv block), and the 1/2 (ov block).  Now we can
C     accomplish it in one fell swoop by symmetrizing the matrix.
C     Q <-- 1/2 (Q + Q')
C
C     Note:  Since we filled only one of the two ov blocks, they'll
C     end up as 1/2 (Q(i,a) + 0), which is what we want!  Tricky!
C
      Call GA_Symmetrize(g_1PDM)
C
      Return
      End
