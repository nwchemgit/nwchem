      SubRoutine DMat(Geom,Basis,g_Vec,g_Dns,nOcc,nBF)
*
* $Id$
*
      Implicit None
#include "errquit.fh"
*
************************************************************************
*     Include standard input and memory checking common structures.
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
************************************************************************
*     
*---- Arguments of a DMat call
      Integer Geom,Basis
      Integer g_Vec,g_Dns,nBF
      Double Precision nOcc(nBF)
*---- Local variables
      Integer iProc,iLo,iHi,jLo,jHi,lthi,lthj
      Integer g_Scr,k_Scr,l_Scr,iCol
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      Double Precision Scale
      Logical LResult,DeBug
*
************************************************************************
*     
      DeBug = .False.
      iProc = ga_nodeid()
*
*---- Allocate the global array
      g_Scr = ga_create_atom_blocked(Geom, Basis, 'Scr')
      Call GA_Copy(g_Vec,g_Scr)
      Call GA_Sync
*
*---- Find the distribution of the g_Scr
      Call GA_Distribution(g_Scr,iProc,iLo,iHi,jLo,jHi)
      lthi = iHi - iLo + 1
      lthj = jHi - jLo + 1
      If (DeBug) Then
         Write(*,*)' DMat: iProc,iLo,iHi,lthi ',iProc,iLo,iHi,lthi
         Write(*,*)' DMat: iProc,jLo,jHi,lthj ',iProc,jLo,jHi,lthj
      End If
      If (lthi.eq.0 . or . lthj.eq.0) GoTo 9999
      If (lthi.lt.0 . or . lthj.lt.0)
     &   Call ErrQuit('DMat: Wrong distribution !!!',0, UNKNOWN_ERR)
*
*---- Allocate the local vector
      LResult = MA_Push_Get(MT_Dbl,lthi,'Vector',l_Scr,k_Scr)
      Call ErrMem(LResult,'DMat: Failed to allocate Vector',0)
*
*---- Scale each column of g_Scr by corresponding occupation number
      Do iCol = jLo, jHi
         Call GA_Get(g_Scr,iLo,iHi,iCol,iCol,Dbl_mb(k_Scr),lthi)
         Scale = nOcc(iCol)
         Call dscal(lthi,Scale,Dbl_mb(k_Scr),1)
         Call GA_Put(g_Scr,iLo,iHi,iCol,iCol,Dbl_mb(k_Scr),lthi)
      End Do
*
*---- Deallocate vector
      LResult = MA_Pop_Stack(l_Scr)
      Call ErrMem(LResult,'DMat: Failed to deallocate Vector',0)
*
 9999 Continue
      Call GA_Sync
      Call GA_Dgemm('N','T',nBF,nBF,nBF,1.d0,g_Vec,g_Scr,0.d0,g_Dns)
*
      if (.not. ga_destroy(g_Scr)) call errquit
     &   ('dmat: could not destroy g_Scr', 0, GA_ERR)
*
*---- Normal termination
      Return
      End
