C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCDen_Interm2 -- Form two-index intermediates for density
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine CCDen_Interm2(NOcc, NVir, T1, LDT1, g_T2,
     $   Z1, LDZ1, g_Z2,
     $   Omega1, LDO1, Omega2, LDO2, Omega3, LDO3, Omega4, LDO4)
      Implicit NONE
#include "errquit.fh"
      Integer NOcc, NVir, LDT1, LDZ1, LDO1, LDO2, LDO3, LDO4
      Integer g_T2, g_Z2
      Double Precision T1(LDT1, NVir), Z1(LDT1, NVir)
      Double Precision Omega1(LDO1, NVir), Omega2(LDO2, NVir),
     $   Omega3(LDO3, NVir), Omega4(LDO4, NOcc)
C
C ARGUMENTS
C DESCRIPTION
C     Compute 2-index intermediates used in coupled cluster 1PDM and 
C     2PDM evaluation
C
C     Outputs are small (NVir*NVir at most), so we replicate them.
C
C REFERENCES
C     Alistair P. Rendell and Timothy J. Lee, J. Chem. Phys. 94, 
C     6219--6228 (1991). Especially Eq. 52--55
C
C     Note:  The 2-index [tZ] terms (eq. 56--57) are trivial to compute
C     on the fly and are only used in two places each, so we don't
C     bother computing them now.
C
C INCLUDE FILES
#include "numerical_constants.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer T2_ILo, T2_IHi, T2_JLo, T2_JHi
      Integer H_Scr1, I_Scr1, H_Scr2, I_Scr2
      Integer I, J, M, N
      Logical oprint
C
C     Get print information
C
      oprint = util_print('information', print_low)
C
      Call DGeWr(NOcc, NVir, T1, LDT1, 0, 0, LUOut,
     $   'T1', 80, 'F20.15')
      Call DGeWr(NOcc, NVir, Z1, LDT1, 0, 0, LUOut,
     $   'Z1', 80, 'F20.15')
C
C     **********************
C     * Initialize Outputs *
C     **********************
C
      Call DFill(LDO1*NVir, FP_Zero, Omega1, 1)
      Call DFill(LDO2*NVir, FP_Zero, Omega2, 1)
      Call DFill(LDO3*NVir, FP_Zero, Omega3, 1)
      Call DFill(LDO4*NOcc, FP_Zero, Omega4, 1)
C
C     ***************************
C     * Local patch of T2 & Z2? *
C     ***************************
C     Note:  Not local in the scrictest sense (as required for ga_access)
C
      Call CCSD_T2_MyIJ(GA_NodeID(), NVir, G_T2, T2_ILo, T2_IHi,
     $     T2_JLo, T2_JHi)
      If (.NOT. GA_Compare_Distr(g_T2, g_Z2) ) Call ErrQuit(
     $   'CCDen_Interm2: T2 and Z2 must have same dist.', 0,
     &       GA_ERR)
C
      If (oprint) Write (LuOut, *) 'CCDen_Interm2 ', GA_NodeID(), 
     $   ': T2 region (',
     $   T2_ILo, ':', T2_IHi, ', ', T2_JLo, ':', T2_JHi, ')'
C
C     ********************************************
C     * Omega1(i,a) = sum(m,e) t(i,m,a,e) Z(m,e) *
C     * APR & TJL Eq. 52                         *
C     ********************************************
C     Get the whole vv patch for each t(i,m) and contract it with the
C     appropriate row of Z1 (vector-matrix multiply)
C
      If (.NOT. MA_Push_Get( MT_Dbl, NVir*NVir, 'CC Density scratch 1',
     $   H_Scr1, I_Scr1)) Call ErrQuit(
     $   'CCDen_Interm2: unable to allocate scratch', NVir*NVir, MA_ERR)
C
      Do I = T2_ILo, T2_IHi
         Do M = T2_JLo, T2_JHi
c$$$C
c$$$            If (oprint) Write (LuOut, *) GA_NodeID(), ': Omega1 (', I, ',', M, ')'
c$$$C
            Call GA_Get(g_T2, (I-1)*NVir+1, I*NVir, (M-1)*NVir+1,
     $         M*NVir, Dbl_MB(I_Scr1), NVir)
            Call dgemm('N', 'T', 1, NVir, NVir,
     $         FP_One, Z1(M,1), LDT1, Dbl_MB(I_Scr1), NVir,
     $         FP_One, Omega1(I,1), LDO1)
c$$$C
c$$$            Call DGeWr(1, NVir, Z1(M,1), LDZ1, 0, 0, LUOut,
c$$$     $         'Z1 patch', 80, 'F20.15')
c$$$            Call DGeWr(NVir, NVir, Dbl_MB(I_Scr1), NVir, 0, 0, LuOut,
c$$$     $         'T2 patch', 80, 'F20.15')
c$$$            Call DGeWr(NOcc, NVir, Omega1, LDO1, 0, 0, LUOut,
c$$$     $         'Partial Omega1', 80, 'F20.15')
c$$$C
         EndDo
      EndDo
C
C     Sum the partial results to complete Omega1
C
      Call GA_DGop(Msg_CCDen_1, Omega1, LDO1*NVir, '+')
c$$$C
c$$$      Call DGeWr(NOcc, NVir, Omega1, LDO1, 0, 0, LUOut,
c$$$     $   'Final Omega1', 80, 'F20.15')
C
C     *****************************************************
C     * Omega2(i,a) = sum(m,e) tau(m,i,a,e) Z(m,e)        *
C     *             = Omega1(i,a) + t(m,a) t(i,e) Z(m,e) -  ?  *
C     * APR & TJL Eq. 53                                  *
C     *****************************************************
C
      Do M = T2_ILo, T2_IHi
         Do I = T2_JLo, T2_JHi
c$$$C
c$$$            If (oprint) Write (LuOut, *) GA_NodeID(), ':Omega2 (', M, ',', I, ')'
c$$$C
            Call GA_Get(g_T2, (M-1)*NVir+1, M*NVir, (I-1)*NVir+1,
     $         I*NVir, Dbl_MB(I_Scr1), NVir)
            Call dgemm('N', 'T', 1, NVir, NVir,
     $         FP_One, Z1(M,1), LDT1, Dbl_MB(I_Scr1), NVir,
     $         FP_One, Omega2(I,1), LDO2)
c$$$C
c$$$            Call DGeWr(NVir, NVir, Dbl_MB(I_Scr1), NVir, 0, 0, LuOut,
c$$$     $         'T2 patch', 80, 'F20.15')
c$$$            Call DGeWr(NOcc, NVir, Omega2, LDO2, 0, 0, LUOut,
c$$$     $         'Partial Omega2', 80, 'F20.15')
c$$$C
         EndDo
      EndDo
C
C     Sum the partial results
C
      Call GA_DGop(Msg_CCDen_2, Omega2, LDO1*NVir, '+')
c$$$C
c$$$      Call DGeWr(NOcc, NVir, Omega1, LDO1, 0, 0, LUOut,
c$$$     $   'Summed Partial Omega2', 80, 'F20.15')
C
C     Add in the t1 t1 z1 terms
C
      Call dgemm('T', 'N', NVir, NVir, NOcc,
     $   FP_One, Z1, LDZ1, T1, LDT1,
     $   FP_Zero, Dbl_MB(I_Scr1), NVir)
c$$$C
c$$$      Call DGeWr(NVir, NVir, Dbl_MB(I_Scr1), NVir, 0, 0, LuOut,
c$$$     $   'Z(m,e) T(m,a)', 80, 'F20.15')
C
      Call dgemm('N', 'N', NOcc, NVir, NVir,
     $   FP_One, T1, LDT1, Dbl_MB(I_Scr1), NVir,
     $   FP_One, Omega2, LDO2)
c$$$C
c$$$      Call DGeWr(NOcc, NVir, Omega2, LDO1, 0, 0, LUOut,
c$$$     $   'Final Omega2', 80, 'F20.15')
C
C     **************************************************
C     * Omega3(a,b) = sum(m,n,e) t(m,n,a,e) z(m,n,b,e) *
C     * APR & TJL Eq. 54                               *
C     **************************************************
C
      If (.NOT. MA_Push_Get( MT_Dbl, NVir*NVir, 'CC Density scratch 2',
     $   H_Scr2, I_Scr2)) Call ErrQuit(
     $   'CCDen_Interm2: unable to allocate scratch', NVir*NVir, MA_ERR)
C
      Do M = T2_ILo, T2_IHi
         Do N = T2_JLo, T2_JHi
c$$$C
c$$$            If (oprint) Write (LuOut, *) GA_NodeID(), ': Omega3 (', M, ',', N, ')'
c$$$C
            Call GA_Get(g_T2, (M-1)*NVir+1, M*NVir, (N-1)*NVir+1,
     $         N*NVir, Dbl_MB(I_Scr1), NVir)
            Call GA_Get(g_Z2, (M-1)*NVir+1, M*NVir, (N-1)*NVir+1,
     $         N*NVir, Dbl_MB(I_Scr2), NVir)
            Call dgemm('N', 'T', NVir, NVir, NVir,
     $         FP_One, Dbl_MB(I_Scr1), NVir, Dbl_MB(I_Scr2), NVir,
     $         FP_One, Omega3, LDO3)
c$$$C
c$$$            Call DGeWr(NVir, NVir, Dbl_MB(I_Scr1), NVir, 0, 0, LuOut,
c$$$     $         'T2 patch', 80, 'F20.15')
c$$$            Call DGeWr(NVir, NVir, Dbl_MB(I_Scr2), NVir, 0, 0, LuOut,
c$$$     $         'Z2 patch', 80, 'F20.15')
c$$$            Call DGeWr(NVir, NVir, Omega3, LDO3, 0, 0, LUOut,
c$$$     $         'Partial Omega3', 80, 'F20.15')
c$$$C
         EndDo
      EndDo
C
      If (.NOT. MA_Pop_Stack( H_Scr2 ) ) Call ErrQuit(
     $   'CCDen_Interm2: unable to free scratch', 2, MA_ERR)
      If (.NOT. MA_Pop_Stack( H_Scr1 ) ) Call ErrQuit(
     $   'CCDen_Interm2: unable to free scratch', 1,
     &       MA_ERR)
C
C     Sum the partial results to complete Omega3
C
      Call GA_DGop(Msg_CCDen_3, Omega3, LDO3*NVir, '+')
c$$$C
c$$$      Call DGeWr(NVir, NVir, Omega3, LDO3, 0, 0, LUOut,
c$$$     $   'Final Omega3', 80, 'F20.15')
C
C     **************************************************
C     * Omega4(i,j) = sum(m,e,f) t(i,m,e,f) z(j,m,e,f) *
C     * APR & TJL Eq. 55                               *
C     **************************************************
C     This amounts to dot products between entire bands of t and z.
C
C     Using GA_DDot_Patch, each node computes its contribution
C     to the dot product and returns immediately, so the Omega4
C     we construct this way is must be summed over all nodes.
C
c     Call GA_Sync
      Do I = 1, NOcc
         Do J = 1, NOcc
            Omega4(i,j) = GA_DDot_Patch(
     $         g_T2, 'N', (I-1)*NVir+1, I*NVir, 1, NOcc*NVir,
     $         g_Z2, 'N', (J-1)*NVir+1, J*NVir, 1, NOcc*NVir)
c$$$C
c$$$            If (oprint) Write (LuOut, *) GA_NodeID(), ':Omega4 (', I, ',', J, ') =',
c$$$     $         Omega4(i,j)
c$$$C
         EndDo
      EndDo
c     Call GA_Sync
C
C     Sum the partial results to complete Omega4
C
c     Call GA_DGop(Msg_CCDen_4, Omega4, LDO4*NOcc, '+')
c$$$C
c$$$      Call DGeWr(NOcc, NOcc, Omega4, LDO4, 0, 0, LUOut,
c$$$     $   'Final Omega4', 80, 'F20.15')
C
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C REVISION
C     $Id$
C SYNOPSIS
      Subroutine CCSD_T2_MyIJ(Node, NVir, g_T2, ILo, IHi, JLo, JHi)
      Implicit NONE
      Integer Node, NVir
      Integer g_T2
      Integer ILo, IHi, JLo, JHi
C
C ARGUMENTS
C DESCRIPTION
C     Determines the range of i and j indices local to the specified
C     process for a T2-like global array of size NVir*NOcc x NVir*NOcc,
C     with the NVir dimension running fastest in both directions.
C
C     The range returned is such that the (ILo,JLo) block is the first
C     complete block for this process, and the lower and right edges
C     (IHi,JLo:JHi) and (ILo:IHi,JHi) reside partly (but not necessarily 
C     completely) on this process.
C
C     Consequently, the ranges returned by this routine are NOT
C     suitable for use with ga_access!
C
C INCLUDE FILES
#include "global.fh"
C
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer IALo, IAHi, JBLo, JBHi
C
      Call GA_Distribution(G_T2, Node, IALo, IAHi, JBLo, JBHi)
C
      If ( IALo .eq. 0 .AND. IAHi .eq. -1) then
         ILo = 0
         IHi = -1
      Else
         ILo = IALo / NVir + 1 ! Lowest I that _starts_ on this proc.
         IHi = IAHi / NVir + 1 ! Must watch out for possible exact match
         If ( Mod(IAHi, NVir) .eq. 0) IHi = IHi - 1
      EndIf
C
      If ( JBLo .eq. 0 .AND. JBHi .eq. -1) then
         JLo = 0
         JHi = -1
      Else
         JLo = JBLo / NVir + 1 ! Lowest J that _starts_ on this proc.
         JHi = JBHi / NVir + 1 ! Must watch out for possible exact match
         If ( Mod(JBHi, NVir) .eq. 0) JHi = JHi - 1
      EndIf
C
      Return
      End
