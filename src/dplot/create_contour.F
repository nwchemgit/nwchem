      SubRoutine Create_Contour(rtdb,geom,AO_bas_han,nbf,nGrid,
     &                          nSet,nSetTri,g_Dns,g_Vec,
     &                          Grid,Values,g_Vls,
     &                          What,Spin,Orb_No)
*
* $Id: create_contour.F,v 1.16 2002-03-04 18:54:21 edo Exp $
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
*---- NWChem handles
      Integer rtdb,Geom,AO_Bas_Han,nbf,nshells_ao

*---- Arguments of a Create_Contour call
      Integer nGrid,nSet,nSetTri,g_Dns,g_Vec(2),g_Vls,Orb_No
      Double Precision Grid(3,nGrid),Values(nGrid)
      Character*8 What,Spin
*---- Geometry variables
      Integer nAtom, iAtom, nxyz, l_xyz, k_xyz
      Integer l_Charge, k_Charge, l_Tags, k_Tags
*---- Used by QDist
      Integer l_rq, k_rq, l_rq0, k_rq0
*---- Basis set variables
      Integer nbf_ao_mxprim, nBF_Max, nCoef_Max
      Double Precision Acc_AO_Gauss
*---- Used by nBF_to_mBF
      Integer l_iniz,k_iniz,l_ifin,k_ifin,l_Expo,k_Expo,l_Coef,k_Coef
      Integer mBF
*---- Used by Int_Eval_GBSets
      Integer lth1, lth2, lth3, nbf_ao_mxang, nScr4
      Integer l_Scr1,k_Scr1,l_Scr2,k_Scr2,l_Scr3,k_Scr3,l_Scr4,k_Scr4
      Integer l_Hess,k_Hess,nbf_ao_mxcont,Max_Elem_Ang_Scr,Max_Component
*---- Used by xc_rhogen
      Double Precision TOLL
      Integer max_at_bf, ilo, ihi, max_at_bf2, lthB, lthZ, nelem, lthDel
      Integer l_FMat,k_FMat,l_PMat,k_PMat,l_BMat,k_BMat,
     &        l_DelRho,k_DelRho,l_gamma,k_gamma,
     &        l_MOVec,k_MOVec
*---- Other variables
      Integer j,iOff,iProc,nProc,nGrPP,lthPP,iOffG,
     .     k_ibf,l_ibf,lrchi_atom, irchi_atom,
     .       lrdelchi_atom, irdelchi_atom,k_cetobfr,l_cetofbr,
     .     k_rdat,l_rdat,lcntoce, icntoce, 
     .     lcntobfr, icntobfr,
     .     iscr,lscr,nscr,nbf_ao_mxnbf_ce,nbf_ao_mxnbf_ce2
cedo      double precision rdens_atom(natom*natom)
      Logical LResult,DeBug,GRAD
*
************************************************************************
*
      DeBug = .False.
      GRAD  = .False.
      iProc = ga_nodeid()
      nProc = ga_nnodes()
*
***      Call dcopy(nGrid,0.0d0,0,Values,1)
      call dfill(ngrid, 0.0d0, values, 1)

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
      LResult = MA_alloc_get(MT_Byte,nAtom*16,'CTags',l_Tags,k_Tags)
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
      if(.not.MA_Pop_Stack(l_Charge))
     .     Call errquit('C_C: Failed to deallocate Charge',0)
!      if(.not.MA_Pop_Stack(l_Tags)) 
!     .     call errquit('C_C: Failed to deallocate CTags',0)
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
      if (.not. bas_nbf_ce_max(AO_bas_han,nbf_ao_mxnbf_ce))
     &     call errquit('crecon: basis set error:', 86)
      If (.not.Bas_NPrim_CN_Max(AO_Bas_Han,nbf_ao_mxprim))
     &   Call ErrQuit('C_C: Bas_NPrim_CN_Max failed',0)
      If (.not.Bas_nBF_CN_Max(AO_Bas_Han,nBF_Max))
     &   Call ErrQuit('C_C: Bas_nBF_CN_Max failed',0)
      If (DeBug . and . iProc.eq.0) Then
         Write(LuOut,*)' nbf_ao_mxprim     = ',nbf_ao_mxprim
         Write(LuOut,*)' nBF_Max    = ',nBF_Max
      End If
*
      nCoef_Max = nBF_Max * nbf_ao_mxprim
      LResult = MA_Push_Get(MT_Int,nAtom,'iniz',l_iniz,k_iniz)
      Call ErrMem(LResult,'C_C: Failed to allocate iniz',0)
      LResult = MA_Push_Get(MT_Int,nAtom,'ifin',l_ifin,k_ifin)
      Call ErrMem(LResult,'C_C: Failed to allocate ifin',0)
      LResult = MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'Expo',l_Expo,k_Expo)
      Call ErrMem(LResult,'C_C: Failed to allocate Expo',0)
      LResult = MA_Push_Get(MT_Dbl,nCoef_Max,'Coef',l_Coef,k_Coef)
      Call ErrMem(LResult,'C_C: Failed to allocate Coef',0)
      LResult = MA_Push_Get(MT_Int,nbf,'ibf',l_ibf,k_ibf)
      Call ErrMem(LResult,'C_C: Failed to allocate ifin',0)
*
      Call dfill(nbf_ao_mxprim,0.0d0,Dbl_mb(k_Expo),1)
      Call dfill(nCoef_Max,0.0d0,Dbl_mb(k_Coef),1)
      Acc_AO_Gauss = 1.d20
      If (What.eq.'VIEW') Acc_AO_Gauss = 1.d20
      mBF = nBF
      Call oldnBF_to_mBF(AO_Bas_Han,dbl_mb(k_rq0),
     &                mBF,nAtom,
     &                Int_mb(k_iniz),Int_mb(k_ifin),
     .     int_mb(k_ibf),
     &                Dbl_mb(k_Expo),nbf_ao_mxprim,Acc_AO_Gauss)
      If (What.eq.'VIEW' . and . nBF.ne.mBF)
     &   Call ErrQuit('Create_Contour: nBF.ne.mBF for VIEW',mbf)
      If (DeBug . and . iProc.eq.0) Write(LuOut,*)' mBF        = ',mBF
*
*---- Compute values of the basis functions and the grid points (use
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
      Call ErrMem(LResult,'C_C: Failed to allocate Scr1',lth1)
      LResult = MA_Push_Get(MT_Dbl,lth2,'Scr2',l_Scr2,k_Scr2)
      Call ErrMem(LResult,'C_C: Failed to allocate Scr2',lth2)
      LResult = MA_Push_Get(MT_Int,1   ,'Hess',l_Hess,k_Hess)
      Call ErrMem(LResult,'C_C: Failed to allocate Hess',1)
      LResult = MA_Push_Get(MT_Int,lth3,'Scr3',l_Scr3,k_Scr3)
      Call ErrMem(LResult,'C_C: Failed to allocate Scr3',lth3)
*
      LResult = Bas_High_Angular(AO_Bas_Han,nbf_ao_mxang)
      LResult = Bas_nContr_CN_Max(AO_bas_han,nbf_ao_mxcont)
!      If (DeBug . and . iProc.eq.0)Write(LuOut,*)' nbf_ao_mxang      = ',nbf_ao_mxang
!      If (DeBug . and . iProc.eq.0)Write(LuOut,*)' nbf_ao_mxcont     = ',nbf_ao_mxcont
*
*---- Allocate an extra scratch space
      Max_Component = ((nbf_ao_mxang+2)*(nbf_ao_mxang+1))/2
      If ( nbf_ao_mxang .eq. 0 ) Then
         Max_Elem_Ang_Scr = Max_Component *  3
      ElseIf ( nbf_ao_mxang .le. 3 ) Then
         Max_Elem_Ang_Scr = Max_Component *  9
      Else
         Max_Elem_Ang_Scr = Max_Component * 28
      Endif 
      nscr4 = 3*ngrpp*nbf_ao_mxcont + 
     &     max(3*ngrpp*nbf_ao_mxprim,
     &     ngrpp*max_elem_ang_scr) + 1
*
      LResult = MA_Push_Get(MT_Dbl,nScr4,'Scr4',l_Scr4,k_Scr4)
      Call ErrMem(LResult,'C_C: Failed to allocate Scr4',0)
*
      Call New_Eval_GBSets(AO_Bas_Han,GRAD,.false.,
     &                     Dbl_mb(k_Scr1),Dbl_mb(k_Scr2),
     &                     Dbl_mb(k_Hess),Int_mb(k_ibf),
     &                     0,Dbl_mb(k_rq0),Dbl_mb(k_rq),
     &                     Grid(1,iOffG),Dbl_mb(k_xyz),
     &                     nGrPP,nGrPP,mBF,nAtom,
     &                     Dbl_mb(k_Scr4),nScr4,
     &                     Dbl_mb(k_Expo),nbf_ao_mxprim,
     &                     Dbl_mb(k_Coef),nCoef_Max,
     &                     Int_mb(k_iniz),Int_mb(k_ifin))

*---- Matrix Scr1(nGrid,nBF) contains HOPEFULY!!! the values of each
*     basis function at a given set of points. The value of the density
*     at a given point (i) is now obtained as a matrix-vector producd
*     Rho(i) = Vi x D x Vi(T), and the values of the orbitals are
*     obtained as Phi(i) = Vi x C, where Vi is the row of the Scr1
*     matrix, D is the density matrix and C is the Phi-th column of the
*     moleculat orbitals.
*
      If (What.eq.'DENSITY') Then
*------- Get the values of the density at the grid points (use xc_rhogen
*        from the NWDFT module)
         TOLL=1.D-15
         max_at_bf = 0
         Do iAtom = 1, nAtom
            LResult = bas_ce2bfr(AO_Bas_Han,iAtom,ilo,ihi)
            max_at_bf = max(max_at_bf,ihi-ilo+1)
         End Do
         max_at_bf2 = max_at_bf*max_at_bf
*
         lthB = nbf_ao_mxnbf_ce*nGrPP
         lthZ=1
         If (GRAD) lthZ = 3*nGrPP*max_at_bf
         nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
*
         LResult = MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENSF',
     .        l_FMat,k_FMat)
         Call ErrMem(LResult,'C_C: Failed to allocate DENSF',
     .       nbf_ao_mxnbf_ce )
         LResult = MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENS',
     .        l_PMat,k_PMat)
         Call ErrMem(LResult,'C_C: Failed to allocate DENS',0)
         LResult = MA_Push_Get(MT_dbl,lthB,'BMat',l_BMat,k_BMat)
         Call ErrMem(LResult,'C_C: Failed to allocate BMat',0)
         call dfill(nbf_ao_mxnbf_ce2, 0.0d0, dbl_mb(k_Fmat), 1)
         call dfill(nbf_ao_mxnbf_ce2, 0.0d0,dbl_mb(k_Pmat), 1)
         call dfill(nbf_ao_mxnbf_ce*ngrpp, 0.0d0,dbl_mb(k_Bmat), 1)
         
*
         lthDel=3*nGrPP
         LResult = MA_Push_Get(MT_dbl,lthDel,'DelRho',l_DelRho,k_DelRho)
         Call ErrMem(LResult,'C_C: Failed to allocate DelRho',0)
         nelem=1
         if(GRAD) nelem=nGrPP
         LResult = MA_Push_Get(MT_Dbl,nelem,'gamma',l_gamma,k_gamma)
         Call ErrMem(LResult,'C_C: Failed to allocate Gamma',0)
         LResult = MA_Push_Get(MT_int,2*nAtom,'cet',l_cetofbr,
     .        k_cetobfr)
         Call ErrMem(LResult,'C_C: Failed to allocate cetof',0)
         LResult = MA_Push_Get(MT_dbl,natom*nAtom,'cet',l_rdat,
     .        k_rdat)
         Call ErrMem(LResult,'C_C: Failed to allocate cetof',0)
         if(.not.bas_numcont(ao_bas_han,nshells_ao))
     .        call errquit(' stop in basnuc',11)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &   call errquit('cre_cont:push_get failed', 13)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit('cre_cont:push_get failed', 13)
         call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &                int_mb(k_cetobfr), natom, nshells_ao)
      call dfill(natom*natom, 0.0d0, dbl_mb(k_rdat), 1)
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &   call errquit('cre_cont: cannot allocate scr',0)
      call util_ga_mat_reduce(nbf, natom, int_mb(k_cetobfr), g_dns, 
     &                        1, dbl_mb(k_rdat), 'rms', 
     &                        dbl_mb(iscr), nbf_ao_mxnbf_ce)
      if (.not.ma_pop_stack(lscr))
     &   call errquit('cre_cont: cannot pop stack',0)
            if (.not.MA_Push_Get(MT_Dbl,natom,'rchi_atom',
     &         lrchi_atom,irchi_atom))
     &         call errquit('cre_cont: cannot allocate rchi_atom',0)
            call dfill(natom, 0.0d0, dbl_mb(irchi_atom), 1)
c
            if (.not.MA_Push_Get(MT_Dbl,natom,'rdelchi_atom',
     &         lrdelchi_atom,irdelchi_atom))call errquit
     &         ('cre_cont: cannot allocate rdelchi_atom',0)
            call dfill(natom, 0.0d0, dbl_mb(irdelchi_atom), 1)
            call util_rmsatbf(ngrpp, natom,
     .           int_mb(k_iniz),int_mb(k_ifin),
     ,           dbl_mb(k_scr1),dbl_mb(irchi_atom))
*
!            call dplot_nobq(natom,byte_mb(k_tags),
!     .           int_mb(k_iniz),int_mb(k_ifin),
!     ,           dbl_mb(irchi_atom),dbl_mb(k_rdat))
*

         call xc_rhogen(TOLL,AO_Bas_Han,g_Dns,
     &                  nbf_ao_mxnbf_ce,nAtom,1,nGrPP,nBF,mBF,
     &                  .false.,GRAD,
     &                  Dbl_mb(k_FMat),Dbl_mb(k_PMat),
     &                  Dbl_mb(k_BMat),
     &                  Dbl_mb(k_Scr1),Dbl_mb(k_Scr2),Int_mb(k_ibf),
     &                  Int_mb(k_iniz),Int_mb(k_ifin),
     &                  Values(iOffG),Dbl_mb(k_DelRho),
     &              dbl_mb(irchi_atom), dbl_mb(irdelchi_atom),
     &              dbl_mb(k_rdat), int_mb(k_cetobfr),1d0)
*------- Deallocate Scratch memory needed by xc_rhogen
         if(.not.ma_chop_stack(l_FMat)) call errquit(
     .        'CreateCont: Failed to deallocate DENSF',0)
      Else
*------- Get the values of the orbital at the grid points
         LResult = MA_Push_Get(MT_Dbl,nBF,'MOVec',l_MOVec,k_MOVec)
         Call ErrMem(LResult,'C_C: Failed to allocate MOVec',0)
*
         Call PhiGen(Dbl_mb(k_Scr1),Dbl_mb(k_MOVec),Values(iOffG),
     &               nGrPP,nBF,g_Vec(1),nSet,Spin,Orb_No)
*
         LResult = MA_Pop_Stack(l_MOVec)
         Call ErrMem(LResult,'C_C: Failed to deallocate MOVec',0)
      End If
*
*---- Deallocate Scratch memory
      if(.not.ma_free_heap(l_tags)) call errquit(' cc: heap',0)
      if(.not.MA_chop_Stack(l_xyz))
     .  Call errquit('C_C: Failed to deallocate xyz',0)
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
*
*-----------------------------------------------------------------------
*
      SubRoutine PhiGen(Chi,Vector,Values,
     &                  nPoint,nBF,g_Vec,nSet,Spin,Orb_No)
      Implicit None
*     
************************************************************************
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "bas.fh"
#include "context.fh"
************************************************************************
*
*---- Arguments of a PhiGen call
      Integer nPoint,nBF,Orb_No,nSet
      Integer g_Vec(nSet)
      Character*8 Spin
      Double Precision Chi(nPoint,nBF),Vector(nBF),Values(nPoint)
*---- Local variables
      Integer jFrom,g_From,iPoint
      Logical DeBug
      Double Precision ddot
      External ddot
*
************************************************************************
*
      DeBug = .False.
*
*---- Pick up the relevant orbital vector
      jFrom = Orb_No
      g_From = g_Vec(1)
      If (nSet.eq.2 . and . Spin.eq.'BETA') g_From = g_Vec(2)
      Call GA_Get(g_From,1,nBF,jFrom,jFrom,Vector,nBF)
*
*---- Loop over the grid points to get the orbital values
      Do iPoint = 1, nPoint
         Values(iPoint) = ddot(nBF,Chi(iPoint,1),nPoint,Vector,1)
*        Write(*,'(d10.4)')Values(iPoint)
      End Do
*
      Call GA_Sync
*
*---- Normal termination
      Return
      End
      Subroutine oldnbf_to_mbf(basis_hand, rq0, mbf, mcenters, 
     &                      iniz, ifin, ibf, zprim, nzprim, acc_gauss)
c
C$Id: create_contour.F,v 1.16 2002-03-04 18:54:21 edo Exp $
c
      implicit none
c      
      integer basis_hand
      integer mcenters ! [input]
      integer nzprim ! [input]
      double precision acc_gauss ! [input]
      integer mbf  ! [input/output]
c
#include "bas.fh"
c
c     Distance Squared between Sampling Points and Centers
c
      double precision rq0(mcenters)
      integer  iniz(mcenters),ifin(mcenters)
      double precision zprim(nzprim)
      integer ncontrset, n1, icset, ictr, nprimo, ncontr, 
     &        isphere, nshbf, l, nang, iprimo
      integer n,npt,m
      integer ibf(*) ! [output]
      double precision zmin
      double precision alpha,logeps,acc_loc
      double precision gaussian_range,r_arg
      gaussian_range(n,logeps,alpha) = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)
c
      acc_loc=-acc_gauss
      if (.not.bas_numcont(basis_hand, ncontrset))
     &   call errquit('Exiting in nbf_to_mbf.',1)
c
      call ifill(mcenters,0,iniz,1)
      call ifill(mcenters,0,ifin,1)
c
      n1 = 0
      npt=0
c
      do 400 icset = 1,ncontrset
         if (.not.bas_cn2ce(basis_hand, icset, ictr))
     &      call errquit('Exiting in nbf_to_mbf.',2)
c         
c        get info about current contraction set
c         
         if (.not.bas_continfo(basis_hand, icset,
     &      l, nprimo, ncontr, isphere))
     &      call errquit('Exiting in nbf_to_mbf.',3)
c         
c        angular momentum
c         
c         
c        cartesian/spherical harmonic
c         
         if (isphere .eq. 0)then !  cartesian set
            nang = (l+1)*(l+2)/2
         elseif (isphere .eq. 1)then !  spherical harmonic
            nang = 2*l+1
         else 
            call errquit('nbf_to_mbf: illegal isphere value', 6)
         endif 
c         
c        get exponents and contraction coefficients for this contraction set
c         
         if (.not.bas_get_exponent(basis_hand, icset, zprim))
     &      call errquit('Exiting in nbf_to_mbf.',7)
c
c        Determine the minimum Gaussian exponent.
c
         zmin = 1.D+06
         do 300 iprimo = 1,nprimo
            zmin = min(zprim(iprimo),zmin)
  300    continue
c
c        Only include those basis functions that are "non-zero" for at least one 
c        point in the sampling set.
c
          if(isphere.eq.1) then
            nshbf=ncontr*(2*l+1)
          else
            nshbf=ncontr*(((l+1)*(l+2))/2)
          endif
          r_arg=rq0(ictr)**2
        if (r_arg.lt.gaussian_range(l,acc_loc,zmin)) then

          do m=1,nshbf
            ibf(n1+m ) = npt+m
          enddo
            if (iniz(ictr).eq.0) iniz(ictr)=n1+1
            n1 = n1 + nshbf
         endif
c
         ifin(ictr)= n1
c
         npt=npt+nshbf
  400 continue
c
      mbf = n1
      return 
      end
c
      subroutine dplot_nobq(mcenters,tags,iniz,ifin,rchi)
      implicit none
#include "inp.fh"
      integer mcenters
      integer  iniz(mcenters),ifin(mcenters)
      double precision rchi(*)
      
c
c     get rid of  bqs in iniz/ifin
c         
      character*16 tags(*),tag
      integer iptr,nleft,i,j
      iptr=1
      nleft=mcenters
      do i=1,mcenters
         tag=tags(i)
         if (inp_compare(.true.,tag(1:2),'bq')) then
            
            do j=iptr+1,nleft
               rchi(j-1)=rchi(j)
               iniz(j-1)=iniz(j)
               ifin(j-1)=ifin(j)
            enddo
            nleft=nleft-1
         else
            iptr=iptr+1
         endif
      enddo
      mcenters=nleft

      return
      end
