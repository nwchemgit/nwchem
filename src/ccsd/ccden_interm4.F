C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     ccden_interm4 -- Form four-index intermediates for density
C
C REVISION
C     $Id: ccden_interm4.F,v 2.1 1997-03-24 02:01:47 gg502 Exp $
C
C SYNOPSIS
      Subroutine CCDen_Interm4(NOcc, NVir,
     $   g_T2, g_Z2, g_Omega1, g_Omega2)
      Implicit NONE
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
C EXTERNAL ROUTINES
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES

C
C     Temporary GA, same size as T2
C
      If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $   'CCDen Scratch', NVir, NVir, G_Scr1) ) Call ErrQuit(
     $   'CCDen_Interm4: unable to allocate scratch',
     $   NVir*NOcc*NVir*NOcc)
C
C     Figure out which patch of Scr1 this processor has
C
      Call CCSD_T2_MyIJ(NVir, G_Scr1, ILo, IHi, JLo, JHi)
C
      Write (LuOut, *) GA_NodeID(), ': Scr1 region (',
     $   ILo, ':', IHi, ', ', JLo, ':', JHi, ')'
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
      If (.NOT. MA_Push_Get(MT_Dbl, NVir*NVir,
     $   'CCDen Scratch 2', H_Scr2, I_Scr2)) Call ErrQuit(
     $   'CCDen_Interm4: unable to allocate scratch', NVir*NVir)
C
C     Get the T2 patches that contribute to the local part of Scr1
C     and accumulate them into their proper places.
C
      Do J = ILo, IHi
         Do L = JLo, JHi
            Call GA_Get(g_T2, (L-1)*NVir+1, L*NVir, (J-1)*NVir+1,
     $         J*NVir, Dbl_MB(I_Scr2), NVir)
            Call GA_Acc(g_Scr1, (J-1)*NVir+1, J*NVir, (L-1)*NVir+1,
     $         L*NVir, Dbl_MB(I_Scr2), NVir, -1.0d0)
         EndDo
      EndDo
C
      Return
      End
