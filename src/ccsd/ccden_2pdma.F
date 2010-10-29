C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCDen_2PDMoooo -- Form coupled cluster 2-particle density matrix
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine CCDen_2PDMoooo(NOcc, NVir, T1, LDT1, g_T2,
     $     g_Z2, g_Q, g_G)
      Implicit NONE
#include "errquit.fh"
      Integer NOcc, NVir, LDT1
      Double Precision T1(LDT1, NVir)
      Integer g_T2, g_Z2, g_Q, g_G
C
C ARGUMENTS
C DESCRIPTION
C     Compute the oooo part of the coupled cluster 2-PDM
C
C REFERENCES
C     Alistair P. Rendell and Timothy J. Lee, J. Chem. Phys. 94, 
C     6219--6228 (1991). Especially Eq. 46
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
C
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer ILo, IHi, JLo, JHi, I, J, K
      Integer h_Qi, i_Qi
      Double precision Qij
      Integer g_Tau
C
C     Clear the result
C
      Call GA_Zero(g_G)
C
C     Find out what portion of G is "local"
C
      Call CCSD_T2_MyIJ(GA_NodeID(), NOcc, g_G, ILo, IHi, JLo, JHi)
C
C     *********************************************************
C     * G(i,j,k,l) = 1/8 P(i,j,k,l) {1/2 [TauZ](i,k,j,l)      *
C     *              + 2 delta(k,l) ( Q(i,j) - delta(i,j) )   *
C     *              -   delta(j,l) ( Q(i,k) - delta(i,k) ) } *
C     * APR & TJL Eq. 46                                      *
C     *********************************************************
C
C     ***********************
C     * 1/2 [TauZ](i,j,k,l) *
C     ***********************
C
      If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $   'CCDen Tau', NVir, NVir, G_Tau) ) Call ErrQuit(
     $   'CCDen_2PDMoooo: unable to allocate Tau',
     $   NVir*NOcc*NVir*NOcc, GA_ERR)
C
      Call CCSD_MkTau_GA(NOcc, NVir, T1, LDT1, g_T2, g_Tau)
      Call CCDen_tZoooo(NOcc, NVir, g_Tau, g_Z2, g_G)
      Call GA_Scale(g_G, 0.50d0)
C
      If (.NOT. GA_Destroy(g_Tau) ) Call ErrQuit(
     $   'CCDen_2PDMoooo: unable to free Tau', 0, GA_ERR)
C
C     Let's do the 1-PDM terms.  Note that the data distribution is 
C     such that each node has complete kl blocks for a particular subset
C     of i and j.
C
C     Note too that the data distribution does not guarantee strict
C     locality -- the edge portion of our patch (determined by
C     CCSD_T2_MyIJ) may in fact be remote.  Thus it is possible this
C     implementation is a bit too naive for good performance.
C     With a little work, it can be made to do fewer comms by using
C     a ga_scatter or other mechanism.
C
C     ****************************************
C     * 2 delta(k,l) ( Q(i,j) - delta(i,j) ) *
C     ****************************************
C     Place the Q(i,j) element all along the diagonal of the block (k=l).
C
      Do I = ILo, IHi
         Do J = JLo, JHi
C
C           Get the appropriate element of Q for this ij block
C
            Call GA_Get(g_Q, I, I, J, J, Qij, 1)
            If ( I .eq. J ) Qij = Qij - 1.0d0
C
C           Now plop it into the diagonal elements of this block
C
            Do K = 1, NOcc
               Call GA_Acc(g_G, (I-1)*NOcc+K, (I-1)*NOcc+K,
     $              (J-1)*NOcc+K, (J-1)*NOcc+K, Qij, 1, 2.0d0)
            EndDo
         EndDo
      EndDo
C
C     **************************************
C     * delta(j,l) ( Q(i,k) - delta(i,k) ) *
C     **************************************
C     Put the Q(i,*) row into one one column (j=l) of the block.
C
      If (.NOT. MA_Push_Get(MT_Dbl, NOcc, 'CCDen Scr',
     $     H_Qi, I_Qi) ) Call ErrQuit(
     $     'CCDen_2PDMoooo: unable to allocate local Qi column', NOcc,
     &       MA_ERR)
C
      Do I = ILo, IHi
C
C        This column of Q goes into one column in each J
C
         Call GA_Get(g_Q, I, I, 1, NOcc, Dbl_MB(I_Qi), 1)
C
C        Handle the delta(i,k) term
C
         Dbl_MB(I_Qi + I - 1) = Dbl_MB(I_Qi + I - 1) - 1.0d0
C
         Do J = JLo, JHi
            Call GA_Acc(g_G, (I-1)*NOcc+1, I*NOcc,
     $           (J-1)*NOcc+J, (J-1)*NOcc+J, Dbl_MB(i_Qi), 1,
     $           -1.0d0)
         EndDo
      EndDo
C
      If ( .NOT. MA_Pop_Stack(H_Qi) ) Call ErrQuit(
     $     'CCDen_2PDMoooo: unable to free Qi', NOcc, MA_ERR)
C
C
C
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCDen_tZoooo -- Form [tZ] intermediate for 2PDM
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine CCDen_tZoooo(NOcc, NVir, g_t, g_Z, g_tZ)
      Implicit NONE
#include "errquit.fh"
      Integer NOcc, NVir
      Integer g_t, g_Z, g_tZ
C
C ARGUMENTS
C DESCRIPTION
C     Computes the [tZ](i,j,k,l) intermediate
C
C     In paper, the term is defined as (Eq. 58)
C     [tZ](i,j,k,l) = sum(a>=b) t(i,j,a,b)(+) Z(k,l,a,b)(+)
C                   + sum(a>b)  t(i,j,a,b)(-) Z(k,l,a,b)(-)
C     with (Eq. 66--68)
C     Z(k,l,a,b)(+/-) = 1/2 [ Z(k,l,a,b) +/- Z(l,k,a,b) ]
C     t(i,j,a,b)(+/-) = 1/2 [ t(i,j,a,b) +/- t(j,i,a,b) ]
C     t(i,j,a,a)(+)   = 1/2 [ t(i,j,a,a) ]
C
C     This can be simplified to
C     [tZ](i,j,k,l) = 1/4 sum(a) [ t(i,j,a,a) {Z(k,l,a,a) + Z(l,k,a,a)} ]
C                   + 1/2 sum(a>b) [ t(i,j,a,b) Z(k,l,a,b) 
C                                  + t(j,i,a,b) Z(l,k,a,b) ]
C
C     Note that [tZ](i,j,k,l) = [tZ](j,i,l,k)
C
C REFERENCES
C     Alistair P. Rendell and Timothy J. Lee, J. Chem. Phys. 94, 
C     6219--6228 (1991). Especially Eq. 58, 66-68.
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
C
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I, J, K, L, ILo, IHi, JLo, JHi, KLo, KHi, LLo, LHi, A, B
      Integer Sorta_Node, Node
      Integer H_Scr1, I_Scr1, H_Scr2, I_Scr2
      Double precision tZ
      Logical oprint
C
C Get print information
C
      oprint = util_print('information',print_low)
C
C     Find out what portion of t is "local"
C
      Call CCSD_T2_MyIJ(GA_NodeID(), NVir, g_t, ILo, IHi, JLo, JHi)
C
      If (oprint)  Write (LuOut, *) 'CCDen_tZoooo ', GA_NodeID(),
     $     ': T2 region (', ILo, ':', IHi, ', ', JLo, ':', JHi, ')'
C
C     Allocate scratch space
C
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scratch 1', H_Scr1, I_Scr1)) Call ErrQuit(
     $   'CCDen_tZoooo: unable to allocate scratch', NVir*NVir, MA_ERR)
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scratch 1', H_Scr2, I_Scr2)) Call ErrQuit(
     $   'CCDen_tZoooo: unable to allocate scratch', NVir*NVir, MA_ERR)
C
C     Basically, we get blocks of t and Z and dot-product them to
C     produce an element [tZ].  For the local chunk of t, contract
C     with all blocks of Z.  The comms can probably be handled smarter.
C     Should check if this is a problem or not.
C
C
C     Loop over processors.  Start with ours and work around the
C     list to avoid everyone pounding on the same processor at once.
C
      Do Sorta_Node = GA_NodeID(), GA_NNodes() + GA_NodeID() - 1
         Node = Mod(Sorta_Node, GA_NNodes() )
C
C        Figure out what Z patches are on the node and loop over them
C
         Call CCSD_T2_MyIJ(Node, NVir, g_Z, KLo, KHi, LLo, LHi)
         If (oprint) Write (LuOut, *) 'CCDen_tZoooo ', GA_NodeID(),
     $        ': Z2 region (', KLo, ':', KHi, ', ', LLo, ':', LHi, ')'
C
         Do K = KLo, KHi
            Do L = LLo, LHi
               Call GA_Get(g_Z, (K-1)*NVir+1, K*NVir, (L-1)*NVir+1,
     $         L*NVir, Dbl_MB(I_Scr1), NVir)
C
C              Contract this Z block with each local t block
C
               Do I = ILo, IHi
                  Do J = JLo, JHi
                     Call GA_Get(g_t, (I-1)*NVir+1, I*NVir,
     $                    (J-1)*NVir+1, J*NVir, Dbl_MB(I_Scr2), NVir)
C
                     tZ = 0.0d0
                     Do A = 1, NVir
                        Do B = 1, A-1
                           tZ = tZ + Dbl_MB(i_Scr1 + (B-1)*NVir + A - 1)
     $                             * Dbl_MB(i_Scr2 + (B-1)*NVir + A - 1)
                        EndDo
C
C                       Now handle diagonal.  Note factor of 1/2 here.
C                       Additional 1/2 applied to final result below.
C
                        tZ = tZ + 0.5d0
     $                       * Dbl_MB(i_Scr1 + (A-1)*NVir + A - 1)
     $                       * Dbl_MB(i_Scr2 + (A-1)*NVir + A - 1)
                     EndDo
C
C                    Don't forget to apply factor of 1/2
C
                     tZ = 0.5d0 * tZ
C
C                    Put this in [tZ](i,j,k,l) and [tZ](j,i,l,k)
C
                     Call GA_Acc(g_tZ, (I-1)*NOcc+K, (I-1)*NOcc+K,
     $                    (J-1)*NOcc+L, (J-1)*NOcc+L, tZ, 1, 1.0d0)
                     Call GA_Acc(g_tZ, (J-1)*NOcc+L, (J-1)*NOcc+L,
     $                    (I-1)*NOcc+K, (I-1)*NOcc+K, tZ, 1, 1.0d0)
                     
                  EndDo
               EndDo
C
            EndDo
         EndDo
C
      EndDo ! Nodes
C
C     Clean up scratch
C
      If (.NOT. MA_Pop_Stack(H_Scr2) ) Call ErrQuit(
     $   'CCDen_tZoooo: unable to free scratch', 0, MA_ERR)
      If (.NOT. MA_Pop_Stack(H_Scr1) ) Call ErrQuit(
     $   'CCDen_tZoooo: unable to free scratch', 0, MA_ERR)
C
C
C
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCDen_tZvvvv -- Form [tZ] intermediate for 2PDM
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine CCDen_tZvvvv(NOcc, NVir, g_t, g_Z, g_tZ)
      Implicit NONE
#include "errquit.fh"
      Integer NOcc, NVir
      Integer g_t, g_Z, g_tZ
C
C ARGUMENTS
C DESCRIPTION
C     Computes the [tZ](a,b,c,d) intermediate
C
C     In paper, the term is defined as (Eq. 59)
C     [tZ](a,b,c,d) = sum(i>=j) t(i,j,a,b)(+) Z(i,j,c,d)(+)
C                   + sum(i>j)  t(i,j,a,b)(-) Z(i,j,c,d)(-)
C     with (Eq. 69--71)
C     Z(i,j,c,d)(+/-) = 1/2 [ Z(i,j,c,d) +/- Z(i,j,d,c) ]
C     t(i,j,a,b)(+/-) = 1/2 [ t(i,j,a,b) +/- t(j,i,a,b) ]
C     t(i,i,a,b)(+)   = 1/2 [ t(i,i,a,b) ]
C
C     This can be simplified to
C     [tZ](a,b,c,d) = 1/4 sum(i) [ t(i,i,a,b) {Z(i,j,c,d) + Z(i,j,d,c)} ]
C                   + 1/2 sum(i>j) [ t(i,j,a,b) Z(i,j,c,d) 
C                                  + t(j,i,a,b) Z(i,j,d,c) ]
C
C     Note that [tZ](a,b,c,d) = [tZ](b,a,d,c)
C
C REFERENCES
C     Alistair P. Rendell and Timothy J. Lee, J. Chem. Phys. 94, 
C     6219--6228 (1991). Especially Eq. 59, 69--71
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
C
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I, J, K, L, ILo, IHi, JLo, JHi, KLo, KHi, LLo, LHi, A, B
      Integer Sorta_Node, Node
      Integer H_Scr1, I_Scr1, H_Scr2, I_Scr2
      Double precision tZ
      Logical oprint
C
C     Get print information
C
      oprint = util_print('information',print_low)
C
C     Find out what portion of t is "local"
C
      Call CCSD_T2_MyIJ(GA_NodeID(), NVir, g_t, ILo, IHi, JLo, JHi)
C
      If (oprint) Write (LuOut, *) 'CCDen_tZvvvv ', GA_NodeID(),
     $     ': T2 region (', ILo, ':', IHi, ', ', JLo, ':', JHi, ')'
C
C     Allocate scratch space
C
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scratch 1', H_Scr1, I_Scr1)) Call ErrQuit(
     $   'CCDen_tZvvvv: unable to allocate scratch', NVir*NVir,
     &       MA_ERR)
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scratch 1', H_Scr2, I_Scr2)) Call ErrQuit(
     $   'CCDen_tZvvvv: unable to allocate scratch', NVir*NVir,
     &       MA_ERR)
C
C     Basically, we get blocks of t and Z and dot-product them to
C     produce an element [tZ].  For the local chunk of t, contract
C     with all blocks of Z.  The comms can probably be handled smarter.
C     Should check if this is a problem or not.
C
C
C     Loop over processors.  Start with ours and work around the
C     list to avoid everyone pounding on the same processor at once.
C
      Do Sorta_Node = GA_NodeID(), GA_NNodes() + GA_NodeID() - 1
         Node = Mod(Sorta_Node, GA_NNodes() )
C
C        Figure out what Z patches are on the node and loop over them
C
         Call CCSD_T2_MyIJ(Node, NVir, g_Z, KLo, KHi, LLo, LHi)
         If (oprint) Write (LuOut, *) 'CCDen_tZvvvv ', GA_NodeID(),
     $        ': Z2 region (', KLo, ':', KHi, ', ', LLo, ':', LHi, ')'
C
         Do K = KLo, KHi
            Do L = LLo, LHi
               Call GA_Get(g_Z, (K-1)*NVir+1, K*NVir, (L-1)*NVir+1,
     $         L*NVir, Dbl_MB(I_Scr1), NVir)
C
C              Contract this Z block with each local t block
C
               Do I = ILo, IHi
                  Do J = JLo, JHi
                     Call GA_Get(g_t, (I-1)*NVir+1, I*NVir,
     $                    (J-1)*NVir+1, J*NVir, Dbl_MB(I_Scr2), NVir)
C
                     tZ = 0.0d0
                     Do A = 1, NVir
                        Do B = 1, A-1
                           tZ = tZ + Dbl_MB(i_Scr1 + (B-1)*NVir + A - 1)
     $                             * Dbl_MB(i_Scr2 + (B-1)*NVir + A - 1)
                        EndDo
C
C                       Now handle diagonal.  Note factor of 1/2 here.
C                       Additional 1/2 applied to final result below.
C
                        tZ = tZ + 0.5d0
     $                       * Dbl_MB(i_Scr1 + (A-1)*NVir + A - 1)
     $                       * Dbl_MB(i_Scr2 + (A-1)*NVir + A - 1)
                     EndDo
C
C                    Don't forget to apply factor of 1/2
C
                     tZ = 0.5d0 * tZ
C
C                    Put this in [tZ](i,j,k,l) and [tZ](j,i,l,k)
C
                     Call GA_Acc(g_tZ, (I-1)*NOcc+K, (I-1)*NOcc+K,
     $                    (J-1)*NOcc+L, (J-1)*NOcc+L, tZ, 1, 1.0d0)
                     Call GA_Acc(g_tZ, (J-1)*NOcc+L, (J-1)*NOcc+L,
     $                    (I-1)*NOcc+K, (I-1)*NOcc+K, tZ, 1, 1.0d0)
                     
                  EndDo
               EndDo
C
            EndDo
         EndDo
C
      EndDo ! Nodes
C
C     Clean up scratch
C
      If (.NOT. MA_Pop_Stack(H_Scr2) ) Call ErrQuit(
     $   'CCDen_tZvvvv: unable to free scratch', 0, MA_ERR)
      If (.NOT. MA_Pop_Stack(H_Scr1) ) Call ErrQuit(
     $   'CCDen_tZvvvv: unable to free scratch', 0, MA_ERR)
C
C
C
      Return
      End
