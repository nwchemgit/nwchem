C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Gamma3 -- Form 3-index Gamma from 3e2c and 4-index Gamma
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_Gamma3(IJSym, IJTri, ILo, IHi, NVirI, JLo, JHi,
     $   NVirJ, G_I, G_J, G_T, D_Gam3, ReadJ, D_J)
      Implicit NONE
#include "errquit.fh"
C
      Logical IJSym
      Logical IJTri
      Integer ILo
      Integer IHi
      Integer NVirI
      Integer JLo
      Integer JHi
      Integer NVirJ
      Integer G_I ![in] Used to size temporary GA for gamma3
      Integer G_J ![in] Integrals to contract with amplitudes
      Integer G_T ![in] Amplitudes
      Integer D_Gam3
      Logical ReadJ ![in] Read J integrals from disk (trashes G_J!)
      Integer D_J ![in] DRA handle to get J integrals from
C
C DESCRIPTION
C     Evaluate G(D,ia) = (D|jb) t(ij,ab)
C
C     Non-zero spin cases:
C     i j   a b
C     a a   a a   G(D,IA) = (D|JB) t(IJ,AB)
C     a b   a b   G(D,IA) = (D|jb) t(Ij,Ab)
C
C     b a   b a   G(D,ia) = (D|JB) t(iJ,aB) -->  (D|JB) t(Ji,Ba)
C           -or-  G(D,jb) = (D|IA) t(Ij,Ab)
C     b b   b b   G(D,ia) = (D|jb) t(ij,ab)
C
C     Note:  Currently assumes that outer loops cover all IJ, even when 
C     there is permutational symmetry.  Such symmetry could be used
C     here by treating pure-spin t like mixed spin, and contracting
C     with both batches of integrals, producing both G(D,ia) and
C     G(D,jb).
C
C PRINT CONTROLS
C     'partial gamma3' will print the G_Gam3 array at the end of the routine.
C
C INCLUDED FILES
#include "global.fh"
#include "dra.fh"
#include "util.fh"
#include "stdio.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer LT1, NI, NJ, Type, NFit, Req
      Integer G_Gam3, ReqJ
C
C     Determine basic dimensions
C
      Call GA_Inquire( G_J, Type, LT1, NFit)
      NI = IHi - ILo +1
      NJ = JHi - JLo +1
C
C     If necessary read integrals into G_J
C
      If ( ReadJ ) then
         If ( DRA_Read_Section(.FALSE.,
     $        G_J, 1,               NVirJ*NJ,  1, NFit,
     $        D_J, 1+NVirJ*(JLo-1), NVirJ*JHi, 1, NFit, ReqJ) .ne. 0)
     $        Call ErrQuit( 'RIMP2_Gamma3: Error reading J integrals',
     $        JLo, DISK_ERR)
         If ( DRA_Wait(ReqJ) .ne. 0) Call ErrQuit(
     $        'RIMP2_Gamma3: Error from DRA_Wait on J integrals', 2,
     &       CALC_ERR)
      EndIf
C
C     Create temporary GA to receive this portion of gamma3
C
      If ( .NOT. GA_Duplicate(G_I, G_Gam3, 'MO 3-index gamma') )
     $     Call ErrQuit('RIMP2_Gamma3: ga_duplicate failed', 0, GA_ERR)
C
C     Since there may be more than one contribution to each patch of 
C     gamma3, we need to _accumulate_ into the DRA.  At present, there
C     is no explicit accumulate function for DRAs, so we'll have to 
C     do it ourselves.  A DRA_Accum_Section routine would save two
C     memory copies.
C
      If ( DRA_Read_Section(.FALSE.,
     $     G_Gam3, 1,               NVirI*NI,  1, NFit,
     $     D_Gam3, 1+NVirI*(ILo-1), NVirI*IHi, 1, NFit, Req) .ne. 0)
     $   Call ErrQuit( 'RIMP2_Gamma3: DRA_Read_Section error.', ILo,
     &       CALC_ERR)
C
      If ( DRA_Wait(Req) .ne. 0) Call ErrQuit(
     $   'RIMP2_Gamma3: Error from DRA_Wait', ILo, CALC_ERR)
C
C     Compute Gam3(ia,D) = t(ia,jb) * (jb|D)
C     Note: t indices here reflect the storage format of t rather than 
C     the traditional notation.
C
      Call GA_DGEMM('N', 'N', NI*NVirI, NFit, NJ*NVirJ,
     $     1.0d0, G_T, G_J, 1.0d0, G_Gam3)
C
      If ( DRA_Write_Section(.FALSE.,
     $     G_Gam3, 1,               NVirI*NI,  1, NFit,
     $     D_Gam3, 1+NVirI*(ILo-1), NVirI*IHi, 1, NFit, Req) .ne. 0)
     $   Call ErrQuit( 'RIMP2_Gamma3: DRA_Write_Section error.', ILo,
     &       CALC_ERR)
C
      If ( Util_Print('partial gamma3', Print_Debug)) then
         Write (LuOut, *) 'After RIMP2_Gamma3...'
         Call GA_Print(g_Gam3)
      EndIf
C
      If ( DRA_Wait(Req) .ne. 0) Call ErrQuit(
     $     'RIMP2_Gamma3: Error from DRA_Wait', ILo, CALC_ERR)
C     
C     Free temporary GA
C
      If ( .NOT. GA_Destroy(G_Gam3) ) Call ErrQuit(
     $     'RIMP2_Gamma3: ga_destroy failed', 0, GA_ERR)
C      
      Return
      End
