      SubRoutine DMat(Geom,Basis,g_Vec,g_Dns,nOcc,nBF)
      Implicit None
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
      Integer g_Scr,k_Scr,l_Scr,iCol
      Double Precision Scale
      Logical LResult
*
************************************************************************
*     
      Call GACrea(Geom,Basis,nBF,nBF,'Scr',g_Scr,'atom')
      Call GA_Copy(g_Vec,g_Scr)
*
      LResult = MA_Push_Get(MT_Dbl,nBF,'Vector',l_Scr,k_Scr)
      Call ErrMem(LResult,'DMat: Failed to allocate Vector',0)
*
*---- PARALLELIZE IN THE FUTURE !!!
      Do iCol = 1, nBF
         Call GA_Get(g_Scr,1,nBF,iCol,iCol,Dbl_mb(k_Scr),nBF)
         Scale = nOcc(iCol)
         Call DScal(nBF,Scale,Dbl_mb(k_Scr),1)
         Call GA_Sync
         Call GA_Put(g_Scr,1,nBF,iCol,iCol,Dbl_mb(k_Scr),nBF)
      End Do
*
      LResult = MA_Pop_Stack(l_Scr)
      Call ErrMem(LResult,'DMat: Failed to deallocate Vector',0)
*
      Call GA_Sync
      Call GA_Dgemm('N','T',nBF,nBF,nBF,1.d0,g_Vec,g_Scr,0.d0,g_Dns)
*
      Call GADest(g_Scr)
*
*---- Normal termination
      Return
      End
