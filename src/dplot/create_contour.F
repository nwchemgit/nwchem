      SubRoutine Create_Contour(rtdb,geom,AO_bas_han,nbf,nGrid,
     &                          nSet,nSetTri,g_Dns,g_Vec,
     &                          Grid,Values,g_Vls,
     &                          What,Spin,Orb_No)
*
* $Id: create_contour.F,v 1.23 2003-01-11 05:29:33 edo Exp $
*
      Implicit None
*
************************************************************************
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
************************************************************************
*
!#define DEBUG 1
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
      Integer nbf_ao_mxprim, nBF_Max
      Double Precision Acc_AO_Gauss
*---- Used by nBF_to_mBF
      Integer l_iniz,k_iniz,l_ifin,k_ifin,l_Expo,k_Expo
      Integer mBF
*---- Used by Int_Eval_GBSets
      Integer lth1
      Integer l_Scr1,k_Scr1
*---- Used by xc_rhogen
      Double Precision TOLL
      Integer   lthB
      Integer l_FMat,k_FMat,l_PMat,k_PMat,l_BMat,k_BMat,
     &        l_MOVec,k_MOVec
*---- Other variables
      Integer j,iOff,iProc,nProc,nGrPP,lthPP,iOffG,
     .     k_ibf,l_ibf,lrchi_atom, irchi_atom,
     .       k_cetobfr,l_cetofbr,
     .     k_rdat,l_rdat,lcntoce, icntoce, 
     .     lcntobfr, icntobfr,
     .     iscr,lscr,nscr,nbf_ao_mxnbf_ce,nbf_ao_mxnbf_ce2
      integer lbas_cent_info, ibas_cent_info,
     &     lbas_cset_info, ibas_cset_info,
     &     ldocset, idocset,i
      integer avail,nstep,ii,ioffg_s,ngrpp_s,needed,jj,
     .     ngrpp_end

            

      Logical LResult,DeBug,GRAD
*
************************************************************************
*
      DeBug = util_print('dplot', print_high)
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
      ngrpp_end=ioffg+ngrpp
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
c
c     get basis set information
c
      if (.not. bas_nbf_ce_max(AO_bas_han,nbf_ao_mxnbf_ce))
     &     call errquit('crecon: basis set error:', 86)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      If (.not.Bas_NPrim_CN_Max(AO_Bas_Han,nbf_ao_mxprim))
     &   Call ErrQuit('C_C: Bas_NPrim_CN_Max failed',0)
      If (.not.Bas_nBF_CN_Max(AO_Bas_Han,nBF_Max))
     &   Call ErrQuit('C_C: Bas_nBF_CN_Max failed',0)
      if(.not.bas_numcont(ao_bas_han,nshells_ao))
     .     call errquit(' stop in basnuc',11)
      If (DeBug . and . iProc.eq.0) Then
         Write(LuOut,*)' nbf_ao_mxprim     = ',nbf_ao_mxprim
         Write(LuOut,*)' nBF_Max    = ',nBF_Max
      End If
*
      LResult = MA_Push_Get(MT_Int,nAtom,'iniz',l_iniz,k_iniz)
      Call ErrMem(LResult,'C_C: Failed to allocate iniz',0)
      LResult = MA_Push_Get(MT_Int,nAtom,'ifin',l_ifin,k_ifin)
      Call ErrMem(LResult,'C_C: Failed to allocate ifin',0)
      LResult = MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'Expo',l_Expo,k_Expo)
      Call ErrMem(LResult,'C_C: Failed to allocate Expo',0)
      LResult = MA_Push_Get(MT_Int,nbf,'ibf',l_ibf,k_ibf)
      Call ErrMem(LResult,'C_C: Failed to allocate ifin',0)
*
      Acc_AO_Gauss = 20d0
      mBF = nBF
      if (.not.MA_Push_Get(mt_log, nshells_ao, 'docset',
     &     ldocset, idocset))
     &     call errquit('grid_quadv0a: cannot allocate docset',
     .     nshells_ao)
      If (DeBug . and . iProc.eq.0) Write(LuOut,*)' mBF        = ',mBF

      if (.not.MA_Push_Get(mt_int, 3*natom, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('grid_quadv0a: cannot allocate bas_cent_info',0)
      if (.not.MA_Push_Get(mt_int, 6*nshells_ao, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('grid_quadv0a: cannot allocate bas_cset_info',0)
      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), natom)

      LResult = MA_Push_Get(MT_int,2*nAtom,'cet',l_cetofbr,
     .     k_cetobfr)
      Call ErrMem(LResult,'C_C: Failed to allocate cetof',0)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &   call errquit('cre_cont:push_get failed', 13)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit('cre_cont:push_get failed', 13)
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &     int_mb(k_cetobfr), natom, nshells_ao)

      If (What.eq.'DENSITY') Then
         if(.not.MA_Push_Get(MT_dbl,natom*natom,'cet',l_rdat,
     .        k_rdat)) call errquit('cc: cant alloc rdat',0)
         call dfill(natom*natom, 0.0d0, dbl_mb(k_rdat), 1)
         nscr = nbf_ao_mxnbf_ce2
         if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &        call errquit('cre_cont: cannot allocate scr',0)
         call util_ga_mat_reduce(nbf, natom, int_mb(k_cetobfr), g_dns, 
     &        1, dbl_mb(k_rdat), 'rms', 
     &        dbl_mb(iscr), nbf_ao_mxnbf_ce)
         if(.not.ma_pop_stack(lscr)) call errquit(
     .        'CreateCont: Failed to deallocate DENSF',0)
*     
         LResult = MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENSF',
     .        l_FMat,k_FMat)
         Call ErrMem(LResult,'C_C: Failed to allocate DENSF',
     .        nbf_ao_mxnbf_ce )
         LResult = MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENS',
     .        l_PMat,k_PMat)
         Call ErrMem(LResult,'C_C: Failed to allocate DENS',0)
         if (.not.MA_Push_Get(MT_Dbl,natom,'rchi_atom',
     &        lrchi_atom,irchi_atom))
     &        call errquit('cre_cont: cannot allocate rchi_atom',0)
         else
            if(.not.MA_Push_Get(MT_Dbl,nBF,'MOVec',l_MOVec,k_MOVec))
     M           call errquit('C_C: Failed to allocate MOVec',0)
         endif
      LResult = MA_Push_Get(MT_Dbl,      nAtom,'rq0',l_rq0,k_rq0)
      Call ErrMem(LResult,'C_C: Failed to allocate rq0',0)
c
c     make it in pieces that use 60% of  MA
c
      avail = MA_inquire_avail(mt_dbl)*50/100
      needed=ngrpp*(nbf*10/9+natom*10/9+nbf_ao_mxnbf_ce)
      nstep=int(dble(needed)/dble(avail))
      nstep=max(nstep,1)
#ifdef DEBUG
!      nstep=5
#endif
      ngrpp_s=nint(dble(ngrpp)/(nstep))
      lth1=ngrpp_s*nbf*10/9
      if(debug.and.iproc.eq.0) then
         write(luout,*) ' no. doubles avail ',avail,' nstep ',nstep
         write(luout,*) ' length ',lth1
      endif
      LResult = MA_Push_Get(MT_Dbl,ngrpp_s*natom*10/9,'rq', l_rq, k_rq)
      Call ErrMem(LResult,'C_C: Failed to allocate rq',0)
      if(.not.MA_Push_Get(MT_Dbl,lth1,'Scr1',l_Scr1,k_Scr1))
     .   Call errquit('C_C: Failed to allocate Scr1',lth1)
      lthB = nbf_ao_mxnbf_ce*ngrpp_s
      LResult = MA_Push_Get(MT_dbl,lthB,'BMat',l_BMat,k_BMat)
      Call ErrMem(LResult,'C_C: Failed to allocate BMat',0)
C
      do 123 ii=1,nstep
         ioffg_s=ioffg+(ii-1)*ngrpp_s
         if(ii.eq.nstep) ngrpp_s = ngrpp_end - ioffg_s
         If (DeBug ) Then
            write(luout,14)  iproc,ii,
     ,           ioffg_s,ngrpp_s,ngrpp_end
 14         format(' me= ',i3,' step# ',i3,' offset: ',i7,
     .           ' stride ',i7,' total  ',i7)
         endif
         Call QDist(Dbl_mb(k_rq0),Dbl_mb(k_rq),
     .        Grid(1,ioffg_s),Dbl_mb(k_xyz),
     &        ngrpp_s,nAtom)
         Call oldnBF_to_mBF(AO_Bas_Han,dbl_mb(k_rq0),
     &        mBF,nAtom,
     &        Int_mb(k_iniz),Int_mb(k_ifin),
     .        int_mb(k_ibf),
     &        Dbl_mb(k_Expo),nbf_ao_mxprim,
     .        log_mb(idocset),Acc_AO_Gauss)
         if(mbf.eq.0) goto 123
#ifdef DEBUG
         do jj=1,ngrpp_s
         write(0,2) ioffg_s+jj-1,grid(1,ioffg_s+jj-1),
     .           dbl_mb(k_rq+jj-1)
 2       format(I3,' galue ',2E20.12)
         enddo
#endif
         call xc_eval_basis(ao_bas_han, 0, dbl_mb(k_scr1),
     &        0,0,0,
     ,        Dbl_mb(k_rq), Grid(1,iOffG_s),Dbl_mb(k_xyz),
     ,        ngrpp_s, natom, int_mb(k_iniz), log_mb(idocset), 
     .        int_mb(ibas_cent_info), int_mb(ibas_cset_info))
#ifdef DEBUG
         do jj=1,ngrpp_s
         write(0,3) ioffg_s+jj-1,
     .           dbl_mb(k_scr1+jj-1),
     .           grid(1,ioffg_s+jj-1),grid(2,ioffg_s+jj-1),
     .           grid(3,ioffg_s+jj-1)
         enddo
 3       format(I3,' gralue ',4(2X,E18.12))
#endif

*---- Matrix Scr1(nGrid,nBF) contains HOPEFULY!!! the values of basis
*     each function at a given set of points. The value of the density
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
*
         call dcopy(nbf_ao_mxnbf_ce2, 0.0d0,0,dbl_mb(k_Fmat), 1)
         call dcopy(nbf_ao_mxnbf_ce2,0.0d0,0,dbl_mb(k_Pmat), 1)
         call dcopy(nbf_ao_mxnbf_ce*ngrpp_s,0.0d0,0,dbl_mb(k_Bmat),1)
         call dcopy(natom, 0.0d0, 0,dbl_mb(irchi_atom), 1)
c
         call util_rmsatbf(ngrpp_s, natom,
     .        int_mb(k_iniz),int_mb(k_ifin),
     ,        dbl_mb(k_scr1),dbl_mb(irchi_atom))

         call xc_rhogen(TOLL,AO_Bas_Han,g_Dns,
     &                  nbf_ao_mxnbf_ce,nAtom,1,ngrpp_s,nBF,mBF,
     &                  .false.,.false.,
     &                  Dbl_mb(k_FMat),Dbl_mb(k_PMat),
     &                  Dbl_mb(k_BMat),
     &                  Dbl_mb(k_Scr1),0,Int_mb(k_ibf),
     &                  Int_mb(k_iniz),Int_mb(k_ifin),
     &                  Values(iOffG_s),0,
     &              dbl_mb(irchi_atom), 0,
     &              dbl_mb(k_rdat), int_mb(k_cetobfr),1d0)
*------- Deallocate Scratch memory needed by xc_rhogen
      Else
*------- Get the values of the orbital at the grid points
*
         Call PhiGen(Dbl_mb(k_Scr1),Dbl_mb(k_MOVec),Values(iOffg_s),
     &               ngrpp_s,nBF,g_Vec(1),nSet,Spin,Orb_No,
     &        mbf,int_mb(k_ibf))
      End If
 123  continue
*
  111 Continue
      if(.not.MA_chop_Stack(l_xyz))
     .  Call errquit('C_C: Failed to deallocate xyz',0)
*
*---- Put values into the global array and synchronize
      Call GA_Put(g_Vls,iOffG,iOffG+nGrPP-1,1,1,Values(iOffG),nGrPP)
 9999 Continue
      if(.not.ma_free_heap(l_tags)) then
         call Ma_summarize_allocated_blocks()
         call errquit(' cc: heap',0)
      endif
      Call GA_Sync
*
*---- Normal termination
      Return
*
*---- Error exit
      Call ErrQuit('Create_Contour: Failed to open output file',0)
      End
*
*-----------------------------------------------------------------------
*
      SubRoutine PhiGen(Chi,Vector,Values,
     &                  nPoint,nBF,g_Vec,nSet,Spin,Orb_No,
     ,     mbf,ibf)
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
      integer mbf,ibf(*)
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
c
c     compress vectors
c
      do ipoint=1,mbf
         vector(ipoint)=vector(ibf(ipoint))
      enddo
*
*---- Loop over the grid points to get the orbital values
      Do iPoint = 1, nPoint
         Values(iPoint) = ddot(mBF,Chi(iPoint,1),nPoint,Vector,1)
*        Write(*,'(d10.4)')Values(iPoint)
      End Do
*
      Call GA_Sync
*
*---- Normal termination
      Return
      End
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
