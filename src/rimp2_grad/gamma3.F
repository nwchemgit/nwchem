C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_Gamma3 -- Form 3-index Gamma from 3e2c and 4-index Gamma
C
C REVISION
C     $Id: gamma3.F,v 1.1 1998-07-06 21:10:00 bernhold Exp $
C
C SYNOPSIS
      Subroutine RIMP2_Gamma3(IJSym, IJTri, ILo, IHi, NVirI, JLo, JHi,
     $   NVirJ, Mixed, G_I, G_J, G_T, D_Gam3)
      Implicit NONE
C
      Logical IJSym
      Logical IJTri
      Integer ILo
      Integer IHi
      Integer NVirI
      Integer JLo
      Integer JHi
      Integer NVirJ
      Logical Mixed
      Integer G_I
      Integer G_J
      Integer G_T
      Integer D_Gam3
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
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer LT1, NI, NJ, Type, NFit, Req
      Integer G_Gam3
C
C     Determine basic dimensions
C
      Call GA_Inquire( G_J, Type, LT1, NFit)
      NI = IHi - ILo +1
      NJ = JHi - JLo +1
C
C     Create temporary GA to receive this portion of gamma3
C
      If ( .NOT. GA_Duplicate(G_I, G_Gam3, 'MO 3-index gamma') )
     $     Call ErrQuit('RIMP2_Gamma3: ga_duplicate failed', 0)
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
     $   Call ErrQuit( 'RIMP2_Gamma3: DRA_Read_Section error.', ILo)
C
      If ( DRA_Wait(Req) .ne. 0) Call ErrQuit(
     $   'RIMP2_Gamma3: Error from DRA_Wait', ILo)
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
     $   Call ErrQuit( 'RIMP2_Gamma3: DRA_Write_Section error.', ILo)
C
      If ( Util_Print('partial gamma3', Print_Debug))
     $     Call GA_Print(g_Gam3)
C
      If ( DRA_Wait(Req) .ne. 0) Call ErrQuit(
     $     'RIMP2_Gamma3: Error from DRA_Wait', ILo)
C     
C     Free temporary GA
C
      If ( .NOT. GA_Destroy(G_Gam3) ) Call ErrQuit(
     $     'RIMP2_Gamma3: ga_destroy failed', 0)
C      
      Return
      End
