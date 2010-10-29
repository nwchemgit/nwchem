C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     CCSD_MkTau_GA -- Form Tau = T2 + T1T1 in a GA
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine CCSD_MkTau_GA(NOcc, NVir, T1, LDT1, g_T2, g_Tau)
      Implicit NONE
#include "errquit.fh"
      Integer NOcc, NVir, LDT1
      Double Precision T1(LDT1, NVir)
      Integer g_T2, g_Tau
C
C ARGUMENTS
C DESCRIPTION
C     Compute Tau(i,j,a,b) = t(i,j,a,b) + t(i,a) t(j,b)
C
C INCLUDE FILES
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
C
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I, J, ILo, IHi, JLo, JHi
      Integer H_Scr1, I_Scr1
C
C     Copy T2 into Tau
C
      Call GA_Copy(g_T2, g_Tau)
C
C     Find out what portion of Tau is "local"
C
      Call CCSD_T2_MyIJ(GA_NodeID(), NVir, g_Tau, ILo, IHi, JLo, JHi)
C
C     Local memory for T1 * T1
C
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCSD Scratch 1', H_Scr1, I_Scr1)) Call ErrQuit(
     $   'CCSD_MkTau_GA: unable to allocate scratch', NVir*NVir, MA_ERR)
C
C     Loop through local portion, adding T1 * T1 updates
C
      Do I = ILo, IHi
         Do J = JLo, JHi
C
C        Form t(i,a) * t(j,b)
C
            Call dgemm('T', 'N', NVir, NVir, 1,
     $           1.0d0, T1(I, 1), LDT1, T1(J, 1), LDT1,
     $           0.0d0, Dbl_MB(i_Scr1), NVir)
            Call GA_Acc(g_Tau, (I-1)*NVir+1, I*NVir,
     $           (J-1)*NVir+1, J*NVir, Dbl_MB(i_Scr1), 1, 1.0d0)
         EndDo
      EndDo
C
C     Clean up scratch
C
      If (.NOT. MA_Pop_Stack(H_Scr1) ) Call ErrQuit(
     $   'CCSD_MkTau_GA: unable to free scratch', 0, MA_ERR)
C
C
C
      Return
      End
