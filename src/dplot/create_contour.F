      SubRoutine Create_Contour(rtdb,Geom,AO_Bas_Han,nBF,nGrid,
     &                          nSet,nSetTri,g_Dns,Grid,Values,g_Vls)
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
*---- NWChem handles
      Integer rtdb,Geom,AO_Bas_Han
*---- Arguments of a Create_Contour call
      Integer nBF,nGrid,nSet,nSetTri,g_Dns,g_Vls,Lth_Tit
      Double Precision Grid(3,nGrid),Values(nGrid)
*---- Geometry variables
      Integer nAtom, iAtom, nxyz, l_xyz, k_xyz
      Integer l_Charge, k_Charge, l_Tags, k_Tags
*---- Used by QDist
      Integer l_rq, k_rq, l_rq0, k_rq0
*---- Basis set variables
      Integer nXPrim, nBF_Max, nCoef_Max
      Double Precision Acc_AO_Gauss
*---- Used by nBF_to_mBF
      Integer l_iniz,k_iniz,l_ifin,k_ifin,l_Expo,k_Expo,l_Coef,k_Coef
      Integer mBF
*---- Used by Int_Eval_GBSets
      Integer lth1, lth2, lth3, nXAng, nScr4
      Integer l_Scr1,k_Scr1,l_Scr2,k_Scr2,l_Scr3,k_Scr3,l_Scr4,k_Scr4
      Integer l_Hess,k_Hess,nXCont,Max_Elem_Ang_Scr,Max_Component
*---- Used by xc_rhogen
      Double Precision TOLL
      Integer max_at_bf, ilo, ihi, max_at_bf2, lthB, lthZ, nelem, lthDel
      Logical ldelley
      Integer l_FMat,k_FMat,l_PMat,k_PMat,l_BMat,k_BMat,l_ZMat,k_ZMat,
     &        l_DelRho,k_DelRho,l_gamma,k_gamma,l_bp,k_bp
*---- Other variables
      Integer i,j,iOff,iProc,nProc,nGrPP,lthPP,iOffG
      Logical LResult,DeBug,GRAD
*
      Logical  Bas_nContr_CN_Max
      External Bas_nContr_CN_Max
*
************************************************************************
*
      DeBug = .False.
      GRAD  = .False.
      iProc = ga_nodeid()
      nProc = ga_nnodes()
*
      Call dcopy(nGrid,0.0d0,0,Values,1)
*---- Share the grid points among processes
      lthPP = nGrid/nProc
      If (lthPP*nProc.ne.nGrid) Then
         lthPP = lthPP + 1
         nGrPP  = lthPP
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
*---- Evaluate distances between each center and the grid points (use
*     QDist from the NWDFT module)
*
      LResult = Geom_Ncent(Geom,nAtom)
      nxyz    = 3*nAtom
      LResult = MA_Push_Get(MT_Dbl,nxyz,'xyz',l_xyz,k_xyz)
      Call ErrMem(LResult,'C_C: Failed to allocate xyz',0)
      LResult = MA_Push_Get(MT_Byte,nAtom*16,'Tags',l_Tags,k_Tags)
      Call ErrMem(LResult,'C_C: Failed to allocate Tags',0)
      LResult = MA_Push_Get(MT_Dbl,nAtom,'Charge',l_Charge,k_Charge)
      Call ErrMem(LResult,'C_C: Failed to allocate Charge',0)
*
      If (.not.Geom_Cart_Get(Geom,nAtom,Byte_mb(k_Tags),
     &   Dbl_mb(k_xyz),Dbl_mb(k_charge)))
     &   Call ErrQuit('C_C: Geom_Cart_Get failed',0)
*
      If (DeBug . and . iProc.eq.0) Then
         Write(LuOut,*)' nAtom      = ',nAtom
         Write(LuOut,*)' Atomic labels and coordinates:'
         Do iAtom = 1, nAtom
            iOff = (iAtom - 1)*3
            Write(LuOut,'(A16,3f12.6)')Byte_mb(k_Tags-1 + iAtom),
     &                                 (Dbl_mb(k_xyz+iOff+j), j = 0, 2)
         End Do
      End If
*
      LResult = MA_Pop_Stack(l_Charge)
      Call ErrMem(LResult,'C_C: Failed to deallocate Charge',0)
      LResult = MA_Pop_Stack(l_Tags)
      Call ErrMem(LResult,'C_C: Failed to deallocate Tags',0)
*
      LResult = MA_Push_Get(MT_Dbl,nGrPP*nAtom,'rq', l_rq, k_rq)
      Call ErrMem(LResult,'C_C: Failed to allocate rq',0)
      LResult = MA_Push_Get(MT_Dbl,      nAtom,'rq0',l_rq0,k_rq0)
      Call ErrMem(LResult,'C_C: Failed to allocate rq0',0)
*
      Call QDist(Dbl_mb(k_rq0),Dbl_mb(k_rq),Grid(1,iOffG),Dbl_mb(k_xyz),
     &           nGrPP,nAtom)
*
*---- Include only basis functions that are "non-zero" for at least one 
*     grid point (use nBF_to_mBF from the NWDFT module)
*
      If (.not.Bas_NPrim_CN_Max(AO_Bas_Han,nXPrim))
     &   Call ErrQuit('C_C: Bas_NPrim_CN_Max failed',0)
      If (.not.Bas_nBF_CN_Max(AO_Bas_Han,nBF_Max))
     &   Call ErrQuit('C_C: Bas_nBF_CN_Max failed',0)
      If (DeBug . and . iProc.eq.0) Then
         Write(LuOut,*)' nXPrim     = ',nXPrim
         Write(LuOut,*)' nBF_Max    = ',nBF_Max
      End If
*
      nCoef_Max = nBF_Max * nXPrim
      LResult = MA_Push_Get(MT_Int,nAtom,'iniz',l_iniz,k_iniz)
      Call ErrMem(LResult,'C_C: Failed to allocate iniz',0)
      LResult = MA_Push_Get(MT_Int,nAtom,'ifin',l_ifin,k_ifin)
      Call ErrMem(LResult,'C_C: Failed to allocate ifin',0)
      LResult = MA_Push_Get(MT_Dbl,nXPrim,'Expo',l_Expo,k_Expo)
      Call ErrMem(LResult,'C_C: Failed to allocate Expo',0)
      LResult = MA_Push_Get(MT_Dbl,nCoef_Max,'Coef',l_Coef,k_Coef)
      Call ErrMem(LResult,'C_C: Failed to allocate Coef',0)
*
      Call dcopy(nXPrim,0.0d0,0,Dbl_mb(k_Expo),1)
      Call dcopy(nCoef_Max,0.0d0,0,Dbl_mb(k_Coef),1)
      Acc_AO_Gauss = 20.d0
      mBF = nBF
      Call nBF_to_mBF(Geom,AO_Bas_Han,
     &                Dbl_mb(k_rq0),mBF,nAtom,
     &                Int_mb(k_iniz),Int_mb(k_ifin),
     &                Dbl_mb(k_Expo),nXPrim,Acc_AO_Gauss)
      If (DeBug . and . iProc.eq.0) Write(LuOut,*)' mBF        = ',mBF
*
*---- Compute values of the basis functions ant the grid points (use
*     Int_Eval_GBSets from the NWDFT module)
*
      lth1 = nGrPP*mBF
      If (lth1.eq.0) lth1 = 1
      lth2 = 1
      If (GRAD) lth2 = nGrPP*mBF*3
      If (lth2.eq.0) lth2 = 1
      lth3 = mBF
      If (lth3.eq.0) lth3 = 1
      LResult = MA_Push_Get(MT_Dbl,lth1,'Scr1',l_Scr1,k_Scr1)
      Call ErrMem(LResult,'C_C: Failed to allocate Scr1',0)
      LResult = MA_Push_Get(MT_Dbl,lth2,'Scr2',l_Scr2,k_Scr2)
      Call ErrMem(LResult,'C_C: Failed to allocate Scr2',0)
      LResult = MA_Push_Get(MT_Int,1   ,'Hess',l_Hess,k_Hess)
      Call ErrMem(LResult,'C_C: Failed to allocate Hess',0)
      LResult = MA_Push_Get(MT_Int,lth3,'Scr3',l_Scr3,k_Scr3)
      Call ErrMem(LResult,'C_C: Failed to allocate Scr3',0)
*
      LResult = Bas_High_Angular(AO_Bas_Han,nXAng)
      LResult = Bas_nContr_CN_Max(AO_bas_han,nXCont)
      If (DeBug . and . iProc.eq.0)Write(LuOut,*)' nXAng      = ',nXAng
      If (DeBug . and . iProc.eq.0)Write(LuOut,*)' nXCont     = ',nXCont
*
*---- Allocate an extra scratch space
      Max_Component = ((nXAng+2)*(nXAng+1))/2
      If ( nXang .eq. 0 ) Then
         Max_Elem_Ang_Scr = Max_Component *  3
      ElseIf ( nXang .le. 3 ) Then
         Max_Elem_Ang_Scr = Max_Component * 9
      Else
         Max_Elem_Ang_Scr = Max_Component * 28
      Endif 
      nScr4=2*nGrPP*nXCont+ Max(2*nGrPP*nXPrim,nGrPP*Max_Elem_Ang_Scr)+1
      If (GRAD) nScr4 = 3*nGrPP*nXCont+
     &          Max(3*nGrPP*nXPrim,nGrPP*Max_Elem_Ang_Scr)+1
*
      LResult = MA_Push_Get(MT_Dbl,nScr4,'Scr4',l_Scr4,k_Scr4)
      Call ErrMem(LResult,'C_C: Failed to allocate Scr4',0)
*
      Call New_Eval_GBSets(Geom,AO_Bas_Han,GRAD,.false.,
     &                     Dbl_mb(k_Scr1),Dbl_mb(k_Scr2),
     &                     Dbl_mb(k_Hess),Int_mb(k_Scr3),
     &                     0,Dbl_mb(k_rq0),Dbl_mb(k_rq),
     &                     Grid(1,iOffG),Dbl_mb(k_xyz),
     &                     nGrPP,nGrPP,mBF,nAtom,
     &                     Dbl_mb(k_Scr4),nScr4,
     &                     Dbl_mb(k_Expo),nXPrim,
     &                     Dbl_mb(k_Coef),nCoef_Max,
     &                     Acc_AO_Gauss)
*
*---- Matrix Scr1(nGrid,nBF) contains HOPEFULY!!! the values of each
*     basis function at a given set of points. The value of the density
*     at a given point is now obtained as a matrix-vector producd
*     Rho(i) = Vi x D x Vi(T), where Vi is the row of the Scr1 matrix.
*
*---- Get the values of the density at the grid points (use xc_rhogen
*     from the NWDFT module)
*
      TOLL=1.D-15
      max_at_bf = 0
      Do iAtom = 1, nAtom
         LResult = bas_ce2bfr(AO_Bas_Han,iAtom,ilo,ihi)
         max_at_bf = max(max_at_bf,ihi-ilo+1)
      End Do
      max_at_bf2 = max_at_bf*max_at_bf
      ldelley=.FALSE.
*     If (.not.rtdb_get(rtdb,'dft:delley',mt_log,1,ldelley))
*    &    ldelley=.FALSE.
*
      lthB = max_at_bf*nGrPP
      lthZ=1
      If (GRAD) lthZ = 3*nGrPP*max_at_bf
*
      LResult = MA_Push_Get(MT_dbl,max_at_bf2,'DENSF',l_FMat,k_FMat)
      Call ErrMem(LResult,'C_C: Failed to allocate DENSF',0)
      LResult = MA_Push_Get(MT_dbl,max_at_bf2,'DENS',l_PMat,k_PMat)
      Call ErrMem(LResult,'C_C: Failed to allocate DENS',0)
      LResult = MA_Push_Get(MT_dbl,lthB,'BMat',l_BMat,k_BMat)
      Call ErrMem(LResult,'C_C: Failed to allocate BMat',0)
      LResult = MA_Push_Get(MT_dbl,lthZ,'ZMat',l_ZMat,k_ZMat)
      Call ErrMem(LResult,'C_C: Failed to allocate ZMat',0)
*
      lthDel=3*nGrPP
      LResult = MA_Push_Get(MT_dbl,lthDel,'DelRho',l_DelRho,k_DelRho)
      Call ErrMem(LResult,'C_C: Failed to allocate DelRho',0)
      nelem=1
      if(GRAD) nelem=nGrPP
      LResult = MA_Push_Get(MT_Dbl,nelem,'gamma',l_gamma,k_gamma)
      Call ErrMem(LResult,'C_C: Failed to allocate Gamma',0)
      LResult = MA_Push_Get(MT_dbl,nGrPP*nAtom,'bp',l_bp,k_bp)
      Call ErrMem(LResult,'C_C: Failed to allocate bp',0)
      Call DFill(nGrPP*nAtom,0.d0,Dbl_mb(k_bp),1)
*
      call rhogen(TOLL,AO_Bas_Han,g_Dns,
     &               max_at_bf,nAtom,1,nGrPP,nBF,mBF,
     &               ldelley,GRAD,
     &               Dbl_mb(k_FMat),Dbl_mb(k_PMat),
     &               Dbl_mb(k_BMat),Dbl_mb(k_ZMat),
     &               Dbl_mb(k_Scr1),Dbl_mb(k_Scr2),Int_mb(k_Scr3),
     &               Int_mb(k_iniz),Int_mb(k_ifin),
     &               Values(iOffG),Dbl_mb(k_DelRho),Dbl_mb(k_gamma),
     &               Dbl_mb(k_bp))
*
*---- Deallocate Scratch memory
      LResult = MA_Pop_Stack(l_bp)
      Call ErrMem(LResult,'C_C: Failed to deallocate bp',0)
      LResult = MA_Pop_Stack(l_gamma)
      Call ErrMem(LResult,'C_C: Failed to deallocate gamma',0)
      LResult = MA_Pop_Stack(l_DelRho)
      Call ErrMem(LResult,'C_C: Failed to deallocate DelRho',0)
      LResult = MA_Pop_Stack(l_ZMat)
      Call ErrMem(LResult,'C_C: Failed to deallocate ZMat',0)
      LResult = MA_Pop_Stack(l_BMat)
      Call ErrMem(LResult,'C_C: Failed to deallocate BMat',0)
      LResult = MA_Pop_Stack(l_PMat)
      Call ErrMem(LResult,'C_C: Failed to deallocate DENS',0)
      LResult = MA_Pop_Stack(l_FMat)
      Call ErrMem(LResult,'C_C: Failed to deallocate DENSF',0)
      LResult = MA_Pop_Stack(l_Scr4)
      Call ErrMem(LResult,'C_C: Failed to deallocate Scr4',0)
      LResult = MA_Pop_Stack(l_Scr3)
      Call ErrMem(LResult,'C_C: Failed to deallocate Scr3',0)
      LResult = MA_Pop_Stack(l_Hess)
      Call ErrMem(LResult,'C_C: Failed to deallocate Hess',0)
      LResult = MA_Pop_Stack(l_Scr2)
      Call ErrMem(LResult,'C_C: Failed to deallocate Scr2',0)
      LResult = MA_Pop_Stack(l_Scr1)
      Call ErrMem(LResult,'C_C: Failed to deallocate Scr1',0)
      LResult = MA_Pop_Stack(l_Coef)
      Call ErrMem(LResult,'C_C: Failed to deallocate Coef',0)
      LResult = MA_Pop_Stack(l_Expo)
      Call ErrMem(LResult,'C_C: Failed to deallocate Expo',0)
      LResult = MA_Pop_Stack(l_ifin)
      Call ErrMem(LResult,'C_C: Failed to deallocate ifin',0)
      LResult = MA_Pop_Stack(l_iniz)
      Call ErrMem(LResult,'C_C: Failed to deallocate iniz',0)
      LResult = MA_Pop_Stack(l_rq0)
      Call ErrMem(LResult,'C_C: Failed to deallocate rq0',0)
      LResult = MA_Pop_Stack(l_rq)
      Call ErrMem(LResult,'C_C: Failed to deallocate rq',0)
      LResult = MA_Pop_Stack(l_xyz)
      Call ErrMem(LResult,'C_C: Failed to deallocate xyz',0)
*
*---- Put values into the global array and synchronize
      Call GA_Put(g_Vls,iOffG,iOffG+nGrPP-1,1,1,Values(iOffG),nGrPP)
 9999 Continue
      Call GA_Sync
*
*---- Normal termination
      Return
*
*---- Error exit
  111 Continue
      Call ErrQuit('Create_Contour: Failed to open output file',0)
      End
