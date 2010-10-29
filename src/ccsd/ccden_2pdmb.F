C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCDen_2PDMccoo -- Form coupled cluster 2-particle density matrix
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine CCDen_2PDMvvoo(NOcc, NVir, T1, LDT1, Z1, LDZ1, 
     $     g_Z2, g_Omega1, g_Omega2, g_Q, g_G)
      Implicit NONE
#include "errquit.fh"
      Integer NOcc, NVir, LDT1, LDZ1
      Double Precision T1(LDT1, NVir), Z1(LDZ1, NVir)
      Integer g_Z2, g_Omega1, g_Omega2, g_Q, g_G
C
C ARGUMENTS
C DESCRIPTION
C     Compute the vvoo part of the coupled cluster 2-PDM.
C
C NOTES
C     Output (g_G) is stored in the same way as g_Z2, despite 
C     the labeling as "vvoo"!
C
C REFERENCES
C     Alistair P. Rendell and Timothy J. Lee, J. Chem. Phys. 94, 
C     6219--6228 (1991). Especially Eq. 48
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
C
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer ILo, IHi, JLo, JHi, I, J, L
      Integer H_Scr1, I_Scr1, H_Scr2, I_Scr2, H_Scr3, I_Scr3,
     $     H_Qvv, I_Qvv
C
C     Clear out target
C
      Call GA_Zero(g_G)
C
C     Find out what portion of G is "local"
C
      Call CCSD_T2_MyIJ(GA_NodeID(), NVir, g_G, ILo, IHi, JLo, JHi)
C
C     ********************************************************************
C     * G(a,b,i,j) = 1/8 P(i,j) P(a,b) {-Z(i,a) t(j,b) - Omega3(i,j,a,b) *
C     *              - sum(l,d) t(l,a) t(i,d) Z(l,j,b,d)                 *
C     *              + 2 delta(i,j) Q(a,b)}                              *
C     ********************************************************************
C     
C     ******************
C     * -Z(i,a) t(j,b) *
C     ******************
C     Local memory for Z1 * T1
C
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scr 1', H_Scr1, I_Scr1)) Call ErrQuit(
     $   'CCDen_2PDMvvoo: unable to allocate scratch 1', NVir*NVir,
     &       MA_ERR)
C
C     Loop through local portion, adding Z1 * T1 updates
C
      Do I = ILo, IHi
         Do J = JLo, JHi
C
C        Form Z(i,a) * t(j,b)
C
            Call dgemm('T', 'N', NVir, NVir, 1,
     $           1.0d0, Z1(I, 1), LDZ1, T1(J, 1), LDT1,
     $           0.0d0, Dbl_MB(i_Scr1), NVir)
            Call GA_Acc(g_G, (I-1)*NVir+1, I*NVir,
     $           (J-1)*NVir+1, J*NVir, Dbl_MB(i_Scr1), 1, -1.0d0)
         EndDo
      EndDo
C
C     Clean up scratch
C
      If (.NOT. MA_Pop_Stack(H_Scr1) ) Call ErrQuit(
     $   'CCDen_2PDMvvoo: unable to free scratch 1', 0, MA_ERR)
C
C     ******************************************************
C     * -Omega3 = -1/2{Omega1 + Omega2} (APR & TJL Eq. 62) *
C     ******************************************************
C
      Call GA_Add(-0.5d0, g_Omega1, 1.0d0, g_G, g_G)
      Call GA_Add(-0.5d0, g_Omega2, 1.0d0, g_G, g_G)
C
C     **************************************
C     * -sum(l,d) t(l,a) t(i,d) Z(l,j,b,d) *
C     **************************************
C     We'll do this by forming t(l,a)t(i,d) products, contracting
C     with Z, and accumulating into the local G.
C
C     It seems like there ought to be a more efficient way to do this
C     by doing one contraction at a time, but as far as I can tell,
C     it would involve two GA intermediates (size vooo) and some
C     personalized all-to-all comms.  Off the top of my head, it seems
C     like a lot of bother and I'm not convinced there is so much to
C     gain.  If this term turns out to be a performance problem, it might
C     be worth looking into.
C
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scr 1', H_Scr1, I_Scr1)) Call ErrQuit(
     $   'CCDen_2PDMvvoo: unable to allocate scratch 1', NVir*NVir,
     &       MA_ERR)
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scr 2', H_Scr2, I_Scr2)) Call ErrQuit(
     $   'CCDen_2PDMvvoo: unable to allocate scratch 2', NVir*NVir,
     &       MA_ERR)
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scr 3', H_Scr3, I_Scr3)) Call ErrQuit(
     $   'CCDen_2PDMvvoo: unable to allocate scratch 3', NVir*NVir,
     &       MA_ERR)
C
      Do J = JLo, JHi
         Do L = 1, NOcc
C
C           Get a patch Z(l,j,*,*).  This is not generally local. Each
C           node ends up getting JLo:JHi full columns of Z2.
C
            Call GA_Get(g_Z2, (L-1)*NVir+1, L*NVir,
     $           (J-1)*NVir+1, J*NVir, Dbl_MB(i_Scr1), 1)
C
            Do I = ILo, IHi
C
C              Form t(l,a) * t(i,d)
C
               Call dgemm('T', 'N', NVir, NVir, 1,
     $              1.0d0, T1(L, 1), LDT1, T1(I, 1), LDT1,
     $              0.0d0, Dbl_MB(i_Scr2), NVir)
C
C              Contract [t(l,a)t(i,d)] with Z(l,j,b,d) and plunk it into
C              G(i,j,a,b).  This patch of G should be "local".
C
               Call dgemm('N', 'T', NVir, NVir, NVir,
     $              1.0d0, Dbl_MB(i_Scr2), NVir, Dbl_MB(i_Scr1), NVir,
     $              0.0d0, Dbl_MB(i_Scr3), NVir)
C
               Call GA_Acc(g_G, (I-1)*NVir+1, I*NVir,
     $              (J-1)*NVir+1, J*NVir, Dbl_MB(i_Scr3), 1, -1.0d0)
            EndDo
         EndDo
      EndDo
C
      If (.NOT. MA_Pop_Stack(H_Scr3) ) Call ErrQuit(
     $   'CCDen_2PDMvvoo: unable to free scratch 3', 0, MA_ERR)
      If (.NOT. MA_Pop_Stack(H_Scr2) ) Call ErrQuit(
     $   'CCDen_2PDMvvoo: unable to free scratch 2', 0, MA_ERR)
      If (.NOT. MA_Pop_Stack(H_Scr1) ) Call ErrQuit(
     $   'CCDen_2PDMvvoo: unable to free scratch 1', 0, MA_ERR)
C
C     ***********************
C     * 2 delta(i,j) Q(a,b) *
C     ***********************
C     We could use GA_Add_Patch, but that is collective and we'd much
C     rather do things in parallel.
C
      If ( .NOT. MA_Push_Get(MT_Dbl, NVir*NVir, 'CCDen Scr Q', H_Qvv,
     $     I_Qvv) ) Call ErrQuit(
     $     'CCDen_2PDMvvoo: unable to allocate local Q(v,v)',
     $     NVir*NVir, MA_ERR)
C
      Call GA_Get(g_Q, NOcc+1, NOcc+NVir, NOcc+1, NVir,
     $     Dbl_MB(I_Qvv), 1)
C
C     Loop over the (i,j) diagonal elements we "own"
C
      Do I = Max(ILo, JLo), Min(IHi, JHi)
         Call GA_Acc(g_G, (I-1)*NVir+1, I*NVir, (I-1)*NVir+1, I*NVir,
     $        Dbl_MB(I_Qvv), 1, 2.0d0)
      EndDo
C
      If ( .NOT. MA_Pop_Stack(H_Qvv) ) Call ErrQuit(
     $     'CCDen_2PDMvvoo: unable to free Q(v,v)', NVir*NVir,
     &       MA_ERR)
C
C
C
      Return
      End
