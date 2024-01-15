c
c     qmd_coup
c
      subroutine qmd_coup(rtdb,geom,geom_p,ipol,nstates,dt,nmo,
     $                    nocc,g_vecs,g_vecs_p,coup,swap_occ,swap_vir,
     $                    en,tda,do_tdks,ks_spin)
c
      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
#include "util.fh"
c
#include "qmd_common.fh"
c
      integer rtdb             ! RTDB handle
      integer geom             ! geometry handle
      integer geom_p           ! geometry handle
      integer nstates           ! number of nstates
      integer g_vecs(2)        ! MOvecs handles
      integer g_vecs_p(2)      ! MOvecs handles
      integer ipol
      integer nocc(2)
      integer nmo(2)
      integer swap_occ(nstates)
      integer swap_vir(nstates)
      integer ks_spin
      double precision coup(nstates,nstates)  ! NAC handle
      double precision dt      ! nuclear time step
      double precision en(nstates) !state energies
      logical do_tdks
      logical tda
c
      integer basis(2)
      integer ispin
      integer g_mo(2)           ! MO overlap handles
      integer g_ao              ! AO overlap handle
      integer g_tmp
      integer nbf
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical int_normalize
      external int_normalize
c
      character*30 pname
c
c     preliminaries
      pname='qmd_coup: '
c
c     Create basis set objects
c     previous geometry
      if (.not.bas_create(basis(1),'prevbasis'))
     &     call errquit(pname//'failed to create',0,BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom_p,basis(1),'ao basis'))
     &     call errquit(pname//'failed to load',0,BASIS_ERR)
      if (.not.int_normalize(rtdb,basis(1)))
     &     call errquit(pname//'failed to normalize',0,BASIS_ERR)
c     current geometry
      if (.not.bas_create(basis(2),'currbasis'))
     &     call errquit(pname//'failed to create',0,BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis(2),'ao basis'))
     &     call errquit(pname//'failed to load',0,BASIS_ERR)
      if (.not.int_normalize(rtdb,basis(2)))
     &     call errquit(pname//'failed to normalize',0,BASIS_ERR)
c
c     get number of basis functions
      if (.not.bas_numbf(basis(2),nbf))
     &     call errquit(pname//'failed to get nbf',0,BASIS_ERR)
c
c     create GAs
c     AO overlap
      g_ao=ga_create_atom_blocked(geom,basis(2),'AO ovrlp')
c     MO overlap
      do ispin=1,ipol
       if(.not.ga_create(mt_dbl,nmo(ispin),nmo(ispin),'MO ovrlp',0,0,
     &                   g_mo(ispin)))
     &     call errquit(pname//'failed to create MO ovrlp',0,GA_ERR)
      end do
c     scratch space
      if(.not.ga_create(mt_dbl,nbf,nmo(1),'tmp',0,0,g_tmp))
     &    call errquit(pname//'failed to create tmp',0,GA_ERR)
c
c     initialize integrals
      call int_init(rtdb,2,basis)
c
c     make sure GA is empty
      call ga_zero(g_ao)
c
c     calculate AO overlap between basis sets
      call int_1e_ga(basis(2),basis(1),g_ao,'overlap',.false.)
c
c     clean up integrals
      call int_terminate()
c
c     calculate MO overlaps: S'(t,t-dt) = C^+(t)S(t,t-dt)C(t-dt)
      do ispin=1,ipol
        call ga_dgemm('n','n',nbf,nmo(ispin),nbf,1.d0,g_ao,
     &                g_vecs_p(ispin),0.d0,g_tmp)
        call ga_dgemm('t','n',nmo(ispin),nmo(ispin),nbf,1.d0,
     &                g_vecs(ispin),g_tmp,0.d0,g_mo(ispin))
      end do
c
      if (do_tdks) then
        call qmd_tdks_coup(nstates,dt,nmo(ks_spin),g_mo(ks_spin),
     $                     coup,swap_occ,swap_vir)
      else
        call qmd_tddft_coup(rtdb,nstates,dt,ipol,nmo,nocc,g_mo,en,tda,
     $                      coup)
      end if
c
c     load geometry for current step into handle for previous step in
c     preparation for next step
      if (.not. geom_rtdb_load(rtdb,geom_p,'geometry'))
     &     call errquit(pname//'geom_rtdb_load',geom,RTDB_ERR)
c
c     copy current MOvecs to previous MOvecs GA
      do ispin=1,ipol
        call ga_copy(g_vecs(ispin),g_vecs_p(ispin))
      end do
c
c     destroy basis set objects
      if(.not.bas_destroy(basis(1)))
     &    call errquit(pname//'bas_destroy',0,BASIS_ERR)
      if(.not.bas_destroy(basis(2)))
     &    call errquit(pname//'bas_destroy',0,BASIS_ERR)
c
c     free memory
      if(.not.ga_destroy(g_ao))
     + call errquit(pname//'ga_destroy',g_ao,GA_ERR)
      if(.not.ga_destroy(g_tmp))
     + call errquit(pname//'ga_destroy',g_tmp,GA_ERR)
      do ispin=1,ipol
       if(.not.ga_destroy(g_mo(ispin)))
     +  call errquit(pname//'ga_destroy',g_mo(ispin),GA_ERR)
      end do
c
      call ga_sync()
c
      return
      end
c
c     qmd_tddft_coup
c
      subroutine qmd_tddft_coup(rtdb,nstates,dt,ipol,nmo,nocc,g_mo,en,
     $                          tda,coup)
c
      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
#include "util.fh"
c
#include "qmd_common.fh"
c
      integer rtdb
      integer nstates           ! number of nstates
      integer g_mo(2)          ! MO overlap handles
      integer nocc(2)
      integer nmo(2)
      integer ipol
      double precision coup(nstates,nstates)  ! NAC handle
      double precision dt      ! nuclear time step
      double precision en(nstates)
      logical tda
c
      integer g_cix1(2)
      integer g_cix2(2)
      integer g_ciy1(2)
      integer g_ciy2(2)
      integer g_citmp(2)
      integer nvir(2)
      integer istate, jstate
      integer imo, jmo, kmo
      integer amo, bmo
      integer ici, jci
      integer ispin
      integer dims(3)
      integer chun(3)
      integer cilo1(3)
      integer cihi1(3)
      integer cild1(2)
      integer cilo2(3)
      integer cihi2(3)
      integer cild2(2)
      integer molo1(2)
      integer mohi1(2)
      integer molo2(2)
      integer mohi2(2)
      integer ndims
      integer inntsize, ok
      double precision scr_mo1, scr_mo2
      double precision scr_ci1, scr_ci2
      double precision scr_co
c
      logical oprint_ham
      character*30 pname
c
c     preliminaries
      pname='qmd_tddft_coup: '
c
c     print
c
      oprint_ham = util_print('ham', print_high).and.ga_nodeid().eq.0
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      nvir=(nmo-nocc)
c
c     allocate local arrays
c     GAs
      chun=0
      do ispin=1,ipol
       dims(1)=nocc(ispin)
       dims(2)=nvir(ispin)
       dims(3)=nstates-1
c      CI X vectors at previous time step (X+Y if RPA)
       if (.not.nga_create(mt_dbl,3,dims,'CIX1',chun,g_cix1(ispin)))
     &     call errquit(pname//'failed to create CIX1',0,GA_ERR)
c      CI X vectors at current time step (X+Y if RPA)
       if (.not.nga_create(mt_dbl,3,dims,'CIX2',chun,g_cix2(ispin)))
     &     call errquit(pname//'failed to create CIX2',0,GA_ERR)
       if (.not.tda) then
c       CI X-Y vectors at previous time step
        if (.not.nga_create(mt_dbl,3,dims,'CIY1',chun,g_ciy1(ispin)))
     &      call errquit(pname//'failed to create CIY1',0,GA_ERR)
c       CI X-Y vectors at current time step
        if (.not.nga_create(mt_dbl,3,dims,'CIY2',chun,g_ciy2(ispin)))
     &      call errquit(pname//'failed to create CIY2',0,GA_ERR)
       end if
c      Workspace
       if(.not.ga_create(mt_dbl,nmo(ispin),nmo(ispin),'tmp',0,0,
     &                   g_citmp(ispin)))
     &     call errquit(pname//'failed to create CIovl',0,GA_ERR)
      end do
c
c
      do ispin=1,ipol
       call ga_zero(g_cix1(ispin))
       call ga_zero(g_cix2(ispin))
       if (.not.tda) then
        call ga_zero(g_ciy1(ispin))
        call ga_zero(g_ciy2(ispin))
       end if
      end do
c
c     read CI vectors into memory
      call qmd_civecs_read(rtdb,g_cix1,g_cix2,g_ciy1,g_ciy2,tda,
     $                 nocc,nvir,nstates-1,ipol)
c
      coup(:,:)=0.d0
c
c     Calculate coupings
c
c     Derivatives approximated as:
c      <Xi(t)|(dXj(t)/dt)> ~ [<Xi(t-dt)|Xj(t)> - <Xi(t)|Xj(t-dt)>]/(2dt)
c
c     MO overlap is S(t,t-dt); S(t-dt,t)=(S(t,t-dt))^+
      do ispin=1,ipol
       do jstate=1,nstates-1
c
c       ###Ground-Excited Coupling###
c       [CI^{J}_{jb}(t)S_{jb}(t-dt,t)-CI^{J}_{jb}(t)S_{jb}(t,t-dt)]
c
c       CI^{J}_{jb}(t),j=>occ, b=>vir
        cilo2(3)=jstate; cihi2(3)=jstate
        cilo2(1)=1; cihi2(1)=nocc(ispin)
        cilo2(2)=1; cihi2(2)=nvir(ispin)
c
c       S_{jb}(t-dt,t)=S_{bj}(t,t-dt), j=>occ, b=>vir
        molo2(1)=nocc(ispin)+1; mohi2(1)=nmo(ispin)
        molo2(2)=1; mohi2(2)=nocc(ispin)
c
        scr_co=nga_ddot_patch(g_cix2(ispin),'n',cilo2,cihi2,
     $                        g_mo(ispin),'t',molo2,mohi2)
c
        coup(1,jstate+1)=coup(1,jstate+1)+scr_co
c
c       CI^{J}_{jb}(t),j=>occ, b=>vir
        cilo1(3)=jstate; cihi1(3)=jstate
        cilo1(1)=1; cihi1(1)=nocc(ispin)
        cilo1(2)=1; cihi1(2)=nvir(ispin)
c
c       S_{jb}(t,t-dt), j=>occ, b=>vir
        molo1(1)=1; mohi1(1)=nocc(ispin)
        molo1(2)=nocc(ispin)+1; mohi1(2)=nmo(ispin)
c
        scr_co=nga_ddot_patch(g_cix2(ispin),'n',cilo1,cihi1,
     $                          g_mo(ispin),'n',molo1,mohi1)
c
        coup(1,jstate+1)=coup(1,jstate+1)-scr_co
c
        coup(jstate+1,1)=-coup(1,jstate+1)
c
        do istate=jstate+1,nstates-1
c
c        ###Excited-Excited Coupling###
c        ***CI derivative contribution***
c        [CI_{ia}(t-dt)CI_{ia}(t)-CI_{ia}(t)CI_{ia}(t-dt)]
c
c        CI^{I}_{ia}(t-dt),i=>occ, a=>vir
         cilo1(3)=istate; cihi1(3)=istate
         cilo1(1)=1; cihi1(1)=nocc(ispin)
         cilo1(2)=1; cihi1(2)=nvir(ispin)
c
c        CI^{J}_{ia}(t),i=>occ, a=>vir
         cilo2(3)=jstate; cihi2(3)=jstate
         cilo2(1)=1; cihi2(1)=nocc(ispin)
         cilo2(2)=1; cihi2(2)=nvir(ispin)
c
         scr_co=nga_ddot_patch(g_cix1(ispin),'n',cilo1,cihi1,
     $                         g_cix2(ispin),'n',cilo2,cihi2)
c
         if (.not.tda) then
c
c         [(X^{I}+Y^{I})(X^{J}+Y^{J})+(X^{I}-Y^{I})(X^{J}-Y^{J})]/2
c          =(X^{I}X^{J}-Y^{I}Y^{J})
c
          scr_co=scr_co+nga_ddot_patch(g_ciy1(ispin),'n',cilo1,cihi1,
     $                                 g_ciy2(ispin),'n',cilo2,cihi2)
          scr_co=scr_co/2.d0
c
         end if
c
         coup(istate+1,jstate+1)=coup(istate+1,jstate+1)+scr_co
c
c        CI^{I}_{ia}(t),i=>occ, a=>vir
         cilo2(3)=istate; cihi2(3)=istate
         cilo2(1)=1; cihi2(1)=nocc(ispin)
         cilo2(2)=1; cihi2(2)=nvir(ispin)
c
c        CI^{J}_{ia}(t-dt),i=>occ, a=>vir
         cilo1(3)=jstate; cihi1(3)=jstate
         cilo1(1)=1; cihi1(1)=nocc(ispin)
         cilo1(2)=1; cihi1(2)=nvir(ispin)
c
         scr_co=nga_ddot_patch(g_cix2(ispin),'n',cilo2,cihi2,
     $                         g_cix1(ispin),'n',cilo1,cihi1)
         if (.not.tda) then
c
          scr_co=scr_co+nga_ddot_patch(g_ciy2(ispin),'n',cilo2,cihi2,
     $                                 g_ciy1(ispin),'n',cilo1,cihi1)
          scr_co=scr_co/2.d0
c
         end if
c
         coup(istate+1,jstate+1)=coup(istate+1,jstate+1)-scr_co
c
c        ***First MO derivative contribution***
c        [CI_{ia}(t)CI_{ib}(t)S_{ab}(t-dt,t)-CI_{ia}(t)CI_{ib}(t)S_{ab}(t,t-dt)]
c
c        (CI^{I}_{ia}(t))^+,i=>occ, a=>vir
         cilo1(3)=istate; cihi1(3)=istate
         cilo1(1)=1; cihi1(1)=nvir(ispin)
         cilo1(2)=1; cihi1(2)=nocc(ispin)
c
c        CI^{J}_{ib}(t),i=>occ, b=>vir
         cilo2(3)=jstate; cihi2(3)=jstate
         cilo2(1)=1; cihi2(1)=nocc(ispin)
         cilo2(2)=1; cihi2(2)=nvir(ispin)
c
         molo2(1)=1; mohi2(1)=nvir(ispin)
         molo2(2)=1; mohi2(2)=nvir(ispin)
c
         call nga_matmul_patch('t','n',1.d0,0.d0,g_cix2(ispin),cilo1,
     $                         cihi1,g_cix2(ispin),cilo2,cihi2,
     $                         g_citmp(ispin),molo2,mohi2)
c
         if (.not.tda) then
          call nga_matmul_patch('t','n',1.d0,1.d0,g_ciy2(ispin),cilo1,
     $                          cihi1,g_ciy2(ispin),cilo2,cihi2,
     $                          g_citmp(ispin),molo2,mohi2)
          call ga_scale(g_citmp(ispin),5.d-1)
         end if
c
c        S_{ab}(t-dt,t)=S_{ba}(t,t-dt), a=>vir, b=>vir
         molo2(1)=nocc(ispin)+1; mohi2(1)=nmo(ispin)
         molo2(2)=nocc(ispin)+1; mohi2(2)=nmo(ispin)
c
         scr_co=ga_ddot_patch(g_mo(ispin),'t',molo2(1),mohi2(1),
     $                        molo2(2),mohi2(2),g_citmp(ispin),'n',1,
     $                        nvir(ispin),1,nvir(ispin))
c
         coup(istate+1,jstate+1)=coup(istate+1,jstate+1)+scr_co
c
c        (CI^{I}_{ia}(t))^+,i=>occ, a=>vir
         cilo2(3)=istate; cihi2(3)=istate
         cilo2(1)=1; cihi2(1)=nvir(ispin)
         cilo2(2)=1; cihi2(2)=nocc(ispin)
c
c        CI^{J}_{ib}(t),i=>occ, b=>vir
         cilo1(3)=jstate; cihi1(3)=jstate
         cilo1(1)=1; cihi1(1)=nocc(ispin)
         cilo1(2)=1; cihi1(2)=nvir(ispin)
c
         molo2(1)=1; mohi2(1)=nvir(ispin)
         molo2(2)=1; mohi2(2)=nvir(ispin)
c
         call nga_matmul_patch('t','n',1.d0,0.d0,g_cix2(ispin),cilo2,
     $                         cihi2,g_cix2(ispin),cilo1,cihi1,
     $                         g_citmp(ispin),molo2,mohi2)
c
         if (.not.tda) then
          call nga_matmul_patch('t','n',1.d0,1.d0,g_ciy2(ispin),cilo2,
     $                          cihi2,g_ciy2(ispin),cilo1,cihi1,
     $                          g_citmp(ispin),molo2,mohi2)
          call ga_scale(g_citmp(ispin),5.d-1)
         end if
c
c        S_{ab}(t,t-dt), a=>vir, b=>vir
         molo1(1)=nocc(ispin)+1; mohi1(1)=nmo(ispin)
         molo1(2)=nocc(ispin)+1; mohi1(2)=nmo(ispin)
c
         scr_co=ga_ddot_patch(g_mo(ispin),'n',molo1(1),mohi1(1),
     $                        molo1(2),mohi1(2),g_citmp(ispin),'n',1,
     $                        nvir(ispin),1,nvir(ispin))
c
         coup(istate+1,jstate+1)=coup(istate+1,jstate+1)-scr_co
c
c        ***Second MO derivative contribution***
c        -[CI_{ia}(t)CI_{ja}(t)S_{ji}(t-dt,t)-CI_{ia}(t)CI_{ja}(t)S_{ji}(t,t-dt)]
c
c        CI^{I}_{ia}(t),i=>occ, a=>vir
         cilo1(3)=istate; cihi1(3)=istate
         cilo1(1)=1; cihi1(1)=nocc(ispin)
         cilo1(2)=1; cihi1(2)=nvir(ispin)
c
c        (CI^{J}_{ja}(t))^+,j=>occ, a=>vir
         cilo2(3)=jstate; cihi2(3)=jstate
         cilo2(1)=1; cihi2(1)=nvir(ispin)
         cilo2(2)=1; cihi2(2)=nocc(ispin)
c
         molo2(1)=1; mohi2(1)=nocc(ispin)
         molo2(2)=1; mohi2(2)=nocc(ispin)
c
         call nga_matmul_patch('n','t',1.d0,0.d0,g_cix2(ispin),cilo1,
     $                         cihi1,g_cix2(ispin),cilo2,cihi2,
     $                         g_citmp(ispin),molo2,mohi2)
c
         if (.not.tda) then
          call nga_matmul_patch('n','t',1.d0,1.d0,g_ciy2(ispin),cilo1,
     $                          cihi1,g_ciy2(ispin),cilo2,cihi2,
     $                          g_citmp(ispin),molo2,mohi2)
          call ga_scale(g_citmp(ispin),5.d-1)
         end if
c
c        S_{ji}(t-dt,t)=S_{ij}(t,t-dt), j=>occ, i=>occ
         molo2(1)=1; mohi2(1)=nocc(ispin)
         molo2(2)=1; mohi2(2)=nocc(ispin)
c
         scr_co=ga_ddot_patch(g_mo(ispin),'n',molo2(1),mohi2(1),
     $                        molo2(2),mohi2(2),g_citmp(ispin),'n',1,
     $                        nocc(ispin),1,nocc(ispin))
c
         coup(istate+1,jstate+1)=coup(istate+1,jstate+1)-scr_co
c
c        CI^{I}_{ia}(t),i=>occ, a=>vir
         cilo2(3)=istate; cihi2(3)=istate
         cilo2(1)=1; cihi2(1)=nocc(ispin)
         cilo2(2)=1; cihi2(2)=nvir(ispin)
c
c        (CI^{J}_{ja}(t))^+,j=>occ, a=>vir
         cilo1(3)=jstate; cihi1(3)=jstate
         cilo1(1)=1; cihi1(1)=nvir(ispin)
         cilo1(2)=1; cihi1(2)=nocc(ispin)
c
         molo2(1)=1; mohi2(1)=nocc(ispin)
         molo2(2)=1; mohi2(2)=nocc(ispin)
c
         call nga_matmul_patch('n','t',1.d0,0.d0,g_cix2(ispin),cilo2,
     $                         cihi2,g_cix2(ispin),cilo1,cihi1,
     $                         g_citmp(ispin),molo2,mohi2)
c
         if (.not.tda) then
          call nga_matmul_patch('n','t',1.d0,1.d0,g_ciy2(ispin),cilo2,
     $                          cihi2,g_ciy2(ispin),cilo1,cihi1,
     $                          g_citmp(ispin),molo2,mohi2)
          call ga_scale(g_citmp(ispin),5.d-1)
         end if
c
c        S_{ji}(t,t-dt), j=>occ, i=>occ
         molo1(1)=1; mohi1(1)=nocc(ispin)
         molo1(2)=1; mohi1(2)=nocc(ispin)
c
         scr_co=ga_ddot_patch(g_mo(ispin),'t',molo1(1),mohi1(1),
     $                        molo1(2),mohi1(2),g_citmp(ispin),'n',1,
     $                        nocc(ispin),1,nocc(ispin))
c
         coup(istate+1,jstate+1)=coup(istate+1,jstate+1)+scr_co
c
         coup(jstate+1,istate+1)=-coup(istate+1,jstate+1)
        end do
       end do
      end do
c
      call ga_sync()
      coup(:,:)=coup(:,:)/(2.d0*dt)
      if (oprint_ham) then
       write(*,*) 'HAM'
       do jstate=1,nstates
        do istate=1,nstates
         if (istate.ne.jstate) then
          write(*,'((es16.8),$)') coup(jstate,istate)
         else
          write(*,'((es16.8),$)') en(jstate)
         end if
        end do
        write(*,*)
       end do
      end if
c
      do ispin=1,ipol
       if(.not.ga_destroy(g_cix1(ispin)))
     +  call errquit(pname//'ga_destroy',g_cix1(ispin),GA_ERR)
       if(.not.ga_destroy(g_cix2(ispin)))
     +  call errquit(pname//'ga_destroy',g_cix2(ispin),GA_ERR)
       if (.not.tda) then
        if(.not.ga_destroy(g_ciy1(ispin)))
     +   call errquit(pname//'ga_destroy',g_cix1(ispin),GA_ERR)
        if(.not.ga_destroy(g_ciy2(ispin)))
     +   call errquit(pname//'ga_destroy',g_cix2(ispin),GA_ERR)
       end if
       if(.not.ga_destroy(g_citmp(ispin)))
     +  call errquit(pname//'ga_destroy',g_citmp(ispin),GA_ERR)
      end do
c
      return
      end
c
c     qmd_tdks_coup
c
      subroutine qmd_tdks_coup(nstates,dt,nmo,g_mo,coup,
     $                         swap_occ,swap_vir)
c
      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
#include "util.fh"
c
#include "qmd_common.fh"
c
      integer nstates           ! number of nstates
      integer g_mo             ! MO overlap handles
      integer nmo
      integer ipol
      integer swap_occ(nstates)
      integer swap_vir(nstates)
      double precision coup(nstates,nstates)  ! NAC handle
      double precision dt      ! nuclear time step
c
      integer istate, jstate
      integer imo, jmo
      integer amo, bmo
      double precision scr_mo1, scr_mo2
c
      character*30 pname
c
c     preliminaries
      pname='qmd_tdks_coup: '
c
      coup(:,:)=0.d0
c
      do jstate=2,nstates
       jmo=swap_occ(jstate)
       bmo=swap_vir(jstate)
       if ((jmo.gt.0).and.(bmo.gt.0)) then
        call ga_get(g_mo,jmo,jmo,bmo,bmo,scr_mo1,nmo)
        call ga_get(g_mo,bmo,bmo,jmo,jmo,scr_mo2,nmo)
        coup(1,jstate)=scr_mo2-scr_mo1
        coup(jstate,1)=-(scr_mo2-scr_mo1)
       end if
       do istate=2,nstates
        imo=swap_occ(istate)
        amo=swap_vir(istate)
        jmo=swap_occ(jstate)
        bmo=swap_vir(jstate)
        if ((imo.gt.0).and.(amo.gt.0).and.(jmo.gt.0).and.
     $      (bmo.gt.0)) then
         if ((imo.eq.jmo).and.(amo.ne.bmo)) then
          call ga_get(g_mo,amo,amo,bmo,bmo,scr_mo1,nmo)
          call ga_get(g_mo,bmo,bmo,amo,amo,scr_mo2,nmo)
          coup(istate,jstate)=scr_mo2-scr_mo1
         else if ((imo.ne.jmo).and.(amo.eq.bmo)) then
          call ga_get(g_mo,imo,imo,jmo,jmo,scr_mo1,nmo)
          call ga_get(g_mo,jmo,jmo,imo,imo,scr_mo2,nmo)
          coup(istate,jstate)=-(scr_mo2-scr_mo1)
         else
          coup(istate,jstate)=0.d0
         end if
        end if
       end do !istate
      end do !jstate
c
      coup(:,:)=coup(:,:)/(2.d0*dt)
C      if (ga_nodeid().eq.0) then
C       write(*,*) 'NAC'
C       do jstate=1,nstates
C        do istate=1,nstates
C         write(*,'((f14.8),$)') coup(jstate,istate)
C        end do
C        write(*,*)
C       end do
C      end if
C      call ga_sync()
c
      return
      end
