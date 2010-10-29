C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     ccden_interm4 -- Form four-index intermediates for density
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine CCDen_Interm4(NOcc, NVir,
     $   g_T2, g_Z2, g_Omega1, g_Omega2)
      Implicit NONE
#include "errquit.fh"
      Integer NOcc, NVir
      Integer g_T2, g_Z2, g_Omega1, g_Omega2
C
C ARGUMENTS
C DESCRIPTION
C     Compute 4-index intermediates used in coupled cluster 2PDM
C     evaluation
C
C     Note that Omega3 (Eq. 62) is a simple combination of Omega1 and
C     Omega2. Omega2 is used in only one place in G, and most uses of
C     Omega3 are fairly simple.  Therefore we'll compute just Omega1
C     and Omega2, and leave to later replacing Omega2 with Omega3.
C
C REFERENCES
C     Alistair P. Rendell and Timothy J. Lee, J. Chem. Phys. 94, 
C     6219--6228 (1991). Especially Eq. 60--62
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
      Integer ILo, IHi, JLo, JHi, J, K, L
      Integer G_Scr1, H_Scr2, I_Scr2
      Logical oprint
C
C     Get print information
C
      oprint = util_print('information', print_low)
C
      Call GA_Print(g_T2)
      Call GA_Print(g_Z2)
C
C     ************************
C     * Allocate temporaries *
C     ************************
C
      If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $   'CCDen Scratch', NVir, NVir, G_Scr1) ) Call ErrQuit(
     $   'CCDen_Interm4: unable to allocate scratch',
     $   NVir*NOcc*NVir*NOcc, GA_ERR)
C
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scratch 2', H_Scr2, I_Scr2)) Call ErrQuit(
     $   'CCDen_Interm4: unable to allocate scratch', NVir*NVir, MA_ERR)
C
C     Figure out which patch of Scr1 this processor has
C
      Call CCSD_T2_MyIJ(GA_NodeID(), NVir, G_Scr1, ILo, IHi, JLo, JHi)
C
      If (oprint) Write (LuOut, *) 'CCDen_Interm4 ', GA_NodeID(),
     $     ': Scr1 region (', ILo, ':', IHi, ', ', JLo, ':', JHi, ')'
C
C     ********************************************************************
C     * Omega1(j,k,a,c) = sum(l,d) [2t(j,l,a,d) - t(l,j,a,d)] Z(k,l,c,d) *
C     * APR & TJL Eq. 60                                                 *
C     ********************************************************************
C     If necessary, this can probably be done with less memory.
C
C     First form 2t(j,l,a,d) - t(l,j,a,d) in the scratch GA
C
C     Implemented in a way that required no new GA coding and should be
C     efficient enough.  A "data parallel" (non-collective) version
C     of ga_add_patch could do this in one shot, though.
C
      Call GA_Copy(g_T2, g_Scr1)
      Call GA_Scale(g_Scr1, 2.0d0)
C
C     Get the T2 patches that contribute to the local part of Scr1
C     and accumulate them into their proper places.
C
      Call GA_Sync
      Do J = ILo, IHi
         Do L = JLo, JHi
            Call GA_Get(g_T2, (L-1)*NVir+1, L*NVir, (J-1)*NVir+1,
     $         J*NVir, Dbl_MB(I_Scr2), NVir)
            Call GA_Acc(g_Scr1, (J-1)*NVir+1, J*NVir, (L-1)*NVir+1,
     $         L*NVir, Dbl_MB(I_Scr2), NVir, -1.0d0)
         EndDo
      EndDo
      Call GA_Sync
C
C     Now contract the intermediate with Z2
C
      Call GA_DGEMM('T', 'N', NOcc*NVir, NOcc*NVir, NOcc*NVir,
     $     1.0d0, g_Scr1, g_Z2, 0.0d0, g_Omega1)
C
C     ********************************************************************
C     * Omega2(j,k,a,c) = sum(l,d) [2Z(l,k,c,d) + Z(k,l,c,d)] t(l,j,a,d) *
C     * APR & TJL Eq. 61                                                 *
C     ********************************************************************
C     If necessary, this can probably be done with less memory.
C
C     First form 2Z(l,k,c,d) + Z(k,l,c,d) in the scratch GA
C
C     Implemented in a way that required no new GA coding and should be
C     efficient enough.  A "data parallel" (non-collective) version
C     of ga_add_patch could do this in one shot, though.
C
      Call GA_Copy(g_Z2, g_Scr1)
C
C     Get the T2 patches that contribute to the local part of Scr1
C     and accumulate them into their proper places.
C
      Call GA_Sync
      Do K = ILo, IHi
         Do L = JLo, JHi
            Call GA_Get(g_T2, (L-1)*NVir+1, L*NVir, (K-1)*NVir+1,
     $         K*NVir, Dbl_MB(I_Scr2), NVir)
            Call GA_Acc(g_Scr1, (K-1)*NVir+1, K*NVir, (L-1)*NVir+1,
     $         L*NVir, Dbl_MB(I_Scr2), NVir, 2.0d0)
         EndDo
      EndDo
      Call GA_Sync
C
C     Now contract the intermediate with Z2
C
      Call GA_DGEMM('T', 'N', NOcc*NVir, NOcc*NVir, NOcc*NVir,
     $     1.0d0, g_Scr1, g_t2, 0.0d0, g_Omega2)
C
C     ************************
C     * Clean up temporaries *
C     ************************
C
      If (.NOT. MA_Pop_Stack(H_Scr2) ) Call ErrQuit(
     $   'CCDen_Interm4: unable to free scratch', 0, MA_ERR)
      If (.NOT. GA_Destroy(g_Scr1) ) Call ErrQuit(
     $   'CCDen_Interm4: unable to free scratch', 0, MA_ERR)
C
      Return
      End
