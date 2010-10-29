      SubRoutine Get_Charges(rtdb,Geom,AO_Bas_Han,nBF,nGrid,nAtom,
     &                       g_Dns,Grid,Values,g_Vls,
     &                       XYZ,Charge,Tags,AMt,BMt,Dst,
     &                       STD,Tot_Charge)
*
* $Id$
*
      Implicit None
*
************************************************************************
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
************************************************************************
*
*---- Arguments of a Get_Charges call
      Integer rtdb,Geom,AO_Bas_Han
      Integer nBF,nGrid,nAtom,g_Dns,g_Vls
      Double Precision Grid(3,nGrid),Values(nGrid)
      Double Precision XYZ(3,nAtom),Charge(nAtom)
      Double Precision AMt(nAtom+1,nAtom+1),BMt(nAtom+1)
      Double Precision STD,Tot_Charge
      Double Precision Dst(nGrid,nAtom)
      Character*16 Tags
*---- Local variables
      Integer iAtom,jAtom,iCart
      Integer iPoint,iGrid,iProc,nProc,nGrPP,lthPP,iOffG
      Integer l_Pot,k_Pot,l_Dns,k_Dns,l_Pivot,k_Pivot,Info_var
      Double Precision Sum,SumA,SumB
      Logical LResult,DeBug
*
      Double Precision ddot
      External         ddot
*
************************************************************************
*
      DeBug = .False.
      iProc = ga_nodeid()
      nProc = ga_nnodes()
*
*---- Get distances (qdist computes squares of distances)
      Do iGrid = 1, nGrid
         Do iAtom = 1, nAtom
            Dst(iGrid,iAtom) = Sqrt( Dst(iGrid,iAtom) )
         End Do
      End Do
*
      If (DeBug . and . iProc.eq.0) Then
         Write(*,*)' Cartesian coordinates'
         Do iAtom = 1, nAtom
            Write(*,'(10f12.6)')(XYZ(iCart,iAtom),iCart=1,3)
         End Do
         Write(*,*)' Grid Points'
         Do iGrid = 1, nGrid
            Write(*,'(10f12.6)')(Grid(iCart,iGrid),iCart=1,3)
         End Do
         Write(*,*)' Distances'
         Do iGrid = 1, nGrid
            Write(*,'(10f12.6)')(Dst(iGrid,iAtom),iAtom=1,nAtom)
         End Do
      End If
*
*---- Share the grid points among processes
      lthPP = nGrid/nProc
      If (lthPP*nProc.ne.nGrid) Then
         lthPP = lthPP + 1
         nGrPP = lthPP
         If (iProc.eq.nProc-1) nGrPP = nGrid - (nProc-1)*lthPP
      Else
         nGrPP = lthPP
      End If
      iOffG = iProc*lthPP+1
      If (nGrPP.eq.0) GoTo 9999
*
      If (DeBug . and . iProc.eq.0) Then
         Write(LuOut,*)
         Write(LuOut,*)'        nBF = ',nBF
         Write(LuOut,*)'      nGrid = ',nGrid
         Write(LuOut,*)'      nProc = ',nProc
         Write(LuOut,*)'      iProc = ',iProc
         Write(LuOut,*)'      nGrPP = ',nGrPP
         Write(LuOut,*)
      End If
*
*---- Initialize integrals
      call int_init(rtdb,1,AO_Bas_Han)
*
*---- Allocate memory for the integrals and the density matrix
      LResult = MA_Push_Get(MT_Dbl,nBF*nBF,'Pot',l_Pot,k_Pot)
      Call ErrMem(LResult,'Get_Charges: Failed to allocate Pot',0)
      LResult = MA_Push_Get(MT_Dbl,nBF*nBF,'Dns',l_Dns,k_Dns)
      Call ErrMem(LResult,'Get_Charges: Failed to allocate Pot',0)
*
*---- Get the density matrix
      Call GA_Get(g_Dns,1,nBF,1,nBF,Dbl_mb(k_Dns),nBF)
*
      If (DeBug . and . iProc.eq.0) Then
         Write(*,*)' Charges:'
         Write(*,'(5f10.5)')(Charge(iAtom),iAtom = 1, nAtom)
      End If
*
*---- Loop over the points assigned to the given node
      Do iPoint = 1, nGrPP
         iGrid = iOffG + iPoint - 1
*------- Get the electronic potential at a given point
         Call Pot_At_Pos(Grid(1,iGrid),1.d0,Dbl_mb(k_Pot),nBF,
     &                   AO_Bas_Han,AO_Bas_Han)
         Values(iGrid) = ddot(nBF*nBF,Dbl_mb(k_Dns),1,Dbl_mb(k_Pot),1) 
*------- Get the nuclear potential
         Sum = 0.d0
         Do iAtom = 1, nAtom
            Sum = Sum + Charge(iAtom)/Dst(iGrid,iAtom)
         End Do
         Values(iGrid) = Values(iGrid) + Sum
      End Do
*
*---- Deallocate memory
      LResult = MA_Pop_Stack(l_Dns)
      Call ErrMem(LResult,'Get_Charges: Failed to deallocate Dns',0)
      LResult = MA_Pop_Stack(l_Pot)
      Call ErrMem(LResult,'Get_Charges: Failed to deallocate Pot',0)
*
*---- Put values into the global array and synchronize
      Call GA_Put(g_Vls,iOffG,iOffG+nGrPP-1,1,1,Values(iOffG),nGrPP)
 9999 Continue
      Call GA_Sync
*
*---- Sequential part
      If (iProc.eq.0) Then
         Call GA_Get(g_Vls,1,nGrid,1,1,Values,nGrid)
         Call dfill((nAtom+1)**2,0.d0,AMt,1)
         Call dfill(nAtom+1,     0.d0,BMt,1)
*------- Set up the A-matrix and B-vector
         Do iAtom = 1, nAtom
            SumB = 0.d0
            Do iGrid = 1, nGrid
               SumB = SumB + Values(iGrid)/Dst(iGrid,iAtom)
            End Do
            BMt(iAtom) = SumB
            Do jAtom = 1, iAtom
               SumA = 0.d0
               Do iGrid = 1, nGrid
                  SumA = SumA + 1.d0/Dst(iGrid,iAtom)/Dst(iGrid,jAtom)
               End Do
               AMt(iAtom,jAtom) = SumA
               AMt(jAtom,iAtom) = SumA
            End Do
            AMt(nAtom+1,iAtom) = 1.d0
            AMt(iAtom,nAtom+1) = 1.d0
         End Do
         AMt(nAtom+1,nAtom+1) = 0.d0
         BMt(nAtom+1) = Tot_Charge
*------- DeBug printout
         If (DeBug) Then
            Write(*,*)' A-matrix'
            Do iAtom = 1, nAtom + 1
               Write(*,'(10f12.6)')(AMt(iAtom,jAtom),jAtom=1,nAtom + 1)
            End Do
            Write(*,*)' B-matrix'
            Do iAtom = 1, nAtom + 1
               Write(*,'(10f12.6)')BMt(iAtom)
            End Do
         End If
*------- Solve a set of linear equations to get charges
         Info_var = 0
         LResult = MA_Push_Get(MT_Int,nAtom+1,'Pivot',l_Pivot,k_Pivot)
         Call ErrMem(LResult,'Get_Charges: Failed to allocate Pivot',0)
         Call dgesv(nAtom+1,1,AMt,nAtom+1,Int_mb(k_Pivot),
     &              BMt,nAtom+1,Info_var)
         LResult = MA_Pop_Stack(l_Pivot)
         Call ErrMem(LResult,
     &               'Get_Charges: Failed to deallocate Pivot',0)
*------- Overwrite nuclear charges (for the final printout in DPlot)
         Do iAtom = 1, nAtom
            Charge(iAtom)= BMt(iAtom)
         End Do
*------- Find the deviations
         STD = 0.d0
         Do iGrid = 1, nGrid
            Sum = 0.d0
            Do iAtom = 1, nAtom
               Sum = Sum + BMt(iAtom)/Dst(iGrid,iAtom)
            End Do
*           Write(*,'(3(2x,e12.6))')Values(iGrid),Sum,Values(iGrid)-Sum
            Sum = Values(iGrid)-Sum
            STD = STD + Sum*Sum
         End Do
         STD = Sqrt(STD/Dble(nGrid-1))
      End If
*
*---- Normal termination
      Return
      End
*
*-----------------------------------------------------------------------
*
      SubRoutine Pot_At_Pos(Pos,Charg,Pot,nBF,ibas,jbas)
      Implicit None
#include "errquit.fh"
*
************************************************************************
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
************************************************************************
*
*---- Arguments of a Pot_At_Pos call
      Double Precision Pos(3),Charg ! [input] point & associated charge
      Integer ibas, jbas            ! [input] bra and ket basis sets
      Integer nBF                   ! [input] # of basis functions
      Double Precision Pot(nBF,nBF) ! [output] calculated potential
*---- Local variables
      Integer nshell_i, nshell_j
      Integer ishell, jshell, mem1, max1e
      Integer ilo, ihi, jlo, jhi, idim, jdim, i, j, iFrom
      Integer l_buf, l_scr
      Integer k_buf, k_scr
      Integer iPos, jPos, nInt
*     
*---- Get info about the basis sets
      If (.not. bas_numcont(ibas, nshell_i)) call errquit
     &     ('Pot_At_Pos: bas_numcont failed for ibas', ibas, BASIS_ERR)
      If (.not. bas_numcont(jbas, nshell_j)) call errquit
     &     ('Pot_At_Pos: bas_numcont failed for jbas', jbas, BASIS_ERR)
*     
*---- Allocate necessary local temporary arrays on the stack
      call int_mem_1e(max1e, mem1)
      if (.not. MA_push_get(MT_DBL,max1e,'buf',l_buf,k_buf))
     &     call errquit('Pot_At_Pos: ma failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem1,'scr',l_scr,k_scr))
     &     call errquit('Pot_At_Pos: ma failed', mem1, MA_ERR)
*     
*---- Get the integrals (use Int_1er from NWDFT module)
      Do jShell = 1, nShell_j
         Do iShell = 1, jShell
            If (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) 
     &         call errquit('int_1e_ga: bas_cn2bfr ?', ibas, BASIS_ERR)
            If (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     &         call errquit('int_1e_ga: bas_cn2bfr ?', jbas, BASIS_ERR)
            idim = ihi - ilo + 1
            jdim = jhi - jlo + 1
            nInt = idim*jdim
*---------- Generate the integrals
            Call Int_1er(Pos,1,Charg,
     &                   jbas,jshell,ibas,ishell, 
     &                   mem1,dbl_mb(k_scr),max1e,dbl_mb(k_buf))
*---------- Put the blockinto the Pot array
            iFrom = k_buf-1
            Do j = jlo, jhi
               Do i = ilo, ihi
                  iFrom = iFrom + 1
                  Pot(j,i) = dbl_mb(iFrom)
               End Do
            End Do
         End Do
      End Do
*
      Do iPos = 1, nBF
         Do jPos = 1, iPos - 1
            Pot(jPos,iPos) = Pot(iPos,jPos)
         End Do
      End Do
*     
*---- Deallocate memory     
      If (.not. MA_Pop_Stack(l_scr))
     &   Call ErrQuit('Pot_At_Pos: pop failed', 0,0)
      if (.not. MA_Pop_Stack(l_buf))
     &   Call ErrQuit('Pot_At_Pos: pop failed', 0,0)
*     
      call int_terminate
      Return
      End
