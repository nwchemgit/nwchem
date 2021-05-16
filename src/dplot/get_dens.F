      SubRoutine Get_Dens(Geom,Basis,nSet,nBF,g_Vec,g_Dns,
     &                    Spin,File_Vec,nOrb,k_Orb,iProc)
*
* $Id$
*
      Implicit None
#include "errquit.fh"
*     
************************************************************************
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "bas.fh"
#include "context.fh"
************************************************************************
*
*---- Arguments of a Get_Dens call
      Integer Geom,Basis,nSet,nBF,g_Vec(2),g_Dns(2)
      Integer nOrb,k_Orb,iProc
      Character*8 Spin
      Character*256 File_Vec
*---- Local variables
      Integer iSet,nOcc(2),k_Occ(2),k_Ene(2),l_Occ(2),l_Ene(2)
      Integer iBF
      Logical DeBug,LResult
      Logical  movecs_read
      External movecs_read
*
************************************************************************
*
      DeBug = .False.
*
*---- Allocate memory for orbital energies and occupation numbers
      LResult = MA_Push_Get(MT_Dbl,nBF,'Ene1',l_Ene(1),k_Ene(1))
      Call ErrMem(LResult,'Get_Dens: Failed to allocate orb. en.',0)
      LResult = MA_Push_Get(MT_Dbl,nBF,'Occ1',l_Occ(1),k_Occ(1))
      Call ErrMem(LResult,'Get_Dens: Failed to allocate occ. num.',0)
      If (nSet.eq.2) Then
         LResult = MA_Push_Get(MT_Dbl,nBF,'Ene2',l_Ene(2),k_Ene(2))
         Call ErrMem(LResult,'Get_Dens: Failed to allocate orb. en.',0)
         LResult = MA_Push_Get(MT_Dbl,nBF,'Occ2',l_Occ(2),k_Occ(2))
         Call ErrMem(LResult,'Get_Dens: Failed to allocate occ. num.',0)
      End If
*
*---- Generate the density matrix
      Do iSet = 1, nSet
*------- Read the orbital file
         If (.not.movecs_read(File_Vec,iSet,Dbl_mb(k_Occ(iSet)),
     &      Dbl_mb(k_Ene(iSet)),g_Vec(iSet)))
     &      Call ErrQuit('Get_Dens: movecs_read failed',0, DISK_ERR)
*------- Get the number of occ. orbitals (for DFT_DensM, no longer used)
         nOcc(iSet) = 0
         Do iBF = 0, nBF-1
            If (Dbl_mb(k_Occ(iSet)+iBF).gt.1.d-8)
     &         nOcc(iSet)=nOcc(iSet)+1
         End Do
         If (Debug . and . iProc.eq.0) 
     &      Write(LuOut,*)' nOcc       = ',nOcc(iSet)
*------- Change the orbital occupancies in the restricted case
         If (nSet.eq.1 . and . Spin.ne.'TOTAL')
     &      Call ChangeOcc(Dbl_mb(k_Occ(iSet)),nBF,Spin)
*------- Eliminate orbitals from the calculation of the density
         If (nOrb.ne.999999) Call Eliminate(Int_mb(k_Orb),nOrb,
     &                       Dbl_mb(k_Occ(iSet)),nBF,iProc)
*------- Finally, get the density in AO basis: D(AO) = C x D(MO) x C(T)
         Call DMat(Geom,Basis,
     &             g_Vec(iSet),g_Dns(iSet),Dbl_mb(k_Occ(iSet)),nBF)
      End Do
      Call GA_Sync
*
*---- Form the relevant density matrix
      If (nSet.eq.2) Then
         If (Spin.eq.'TOTAL') Then
            Call GA_dAdd(1.d0,g_Dns(1),1.d0,g_Dns(2),g_Dns(1))
         Else If (Spin.eq.'SPINDENS') Then
            Call GA_dAdd(1.d0,g_Dns(1),-1.d0,g_Dns(2),g_Dns(1))
         Else If (Spin.eq.'ALPHA') Then
         Else If (Spin.eq.'BETA') Then
            Call GA_Copy(g_Dns(2),g_Dns(1))
         End If
      End If
*
*---- Deallocate memory
      If (nSet.eq.2) Then
         LResult = MA_Pop_Stack(l_Occ(2))
         Call ErrMem(LResult,
     &        'Get_Dens: Failed to deallocate occ. num.',0)
         LResult = MA_Pop_Stack(l_Ene(2))
         Call ErrMem(LResult,
     &        'Get_Dens: Failed to deallocate orb. en.',0)
      End If
      LResult = MA_Pop_Stack(l_Occ(1))
      Call ErrMem(LResult,'Get_Dens: Failed to deallocate occ. num.',0)
      LResult = MA_Pop_Stack(l_Ene(1))
      Call ErrMem(LResult,'Get_Dens: Failed to deallocate orb. en.',0)
*
      Call GA_Sync
*
*---- Normal termination
      Return
      End
*
*-----------------------------------------------------------------------
*
      SubRoutine Eliminate(OrbNo,nOrb,OccNo,nBF,iProc)
      Implicit None
#include "errquit.fh"
*---- Arguments of a Eliminate call
      Integer nOrb,nBF,iProc
      Integer OrbNo(nOrb)
      Double Precision OccNo(nBF)
*---- Local Variables
      Integer i,iCount
      Logical DeBug
*
      DeBug = .False.
      If (DeBug . and . iProc.eq.0)
     &   Write(*,'(30(1x,f3.1))') (OccNo(i),i=1,nBF)
*
      iCount = 1
      Do i = 1, nBF
         If (i.eq.OrbNo(iCount)) Then
            iCount = iCount + 1
         Else
            OccNo(i) = 0.d0
         End If
      End Do
      If (iCount-1.ne.nOrb) Call ErrQuit('Eliminate: failed',0,
     &       UNKNOWN_ERR)
*
      If (DeBug . and . iProc.eq.0)
     &   Write(*,'(30(1x,f3.1))') (OccNo(i),i=1,nBF)
*
      Return
      End
*
*-----------------------------------------------------------------------
*
      SubRoutine ChangeOcc(OccNo,nBF,Spin)
      Implicit None
#include "errquit.fh"
*
************************************************************************
*
*---- Arguments of a ChangeOcc call
      Integer nBF
      Double Precision OccNo(nBF)
      Character*8 Spin
*---- Local variables
      Integer iBF
*
************************************************************************
*
      Do iBF = 1, nBF
         If ( Abs(OccNo(iBF)-2.d0) . lt . 1.d-8 ) Then
            If (Spin.eq.'SPINDENS') Then
               OccNo(iBF) = 0.d0
            Else If (Spin.eq.'ALPHA') Then
               OccNo(iBF) = 1.d0
            Else If (Spin.eq.'BETA') Then
               OccNo(iBF) = 1.d0
            End If
         Else If ( Abs(OccNo(iBF)-1.d0) . lt . 1.d-8 ) Then
            If (Spin.eq.'BETA') Then
               OccNo(iBF) = 0.d0
            End If
         Else If ( Abs(OccNo(iBF) ) . lt . 1.d-8 ) Then
         Else
            Call ErrQuit('Cannot handle non-integer occupations yet',0,
     &       CAPMIS_ERR)
         End If
      End Do
*
      Return
      End
