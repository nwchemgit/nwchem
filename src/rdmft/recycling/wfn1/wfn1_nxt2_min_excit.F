C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Optimize the state vector to minimize the total energy
C> of an excited state
C>
C> This routine minimizes the total energy of the state vector 
C> while satisfying the various constraints. The general approach is
C> to do linesearches along the gradient to find the minimum energy
C> in that particular direction. Then the step to that point is taken,
C> the energy and the gradient calculated and a new linesearch started.
C> This procedure continues until the gradient is small.
C>
      subroutine wfn1_nxt2_min_excit(rtdb,geom,basis,nbf,nea,neb,
     &           gstate,state,etot)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "wfn1_nxt.fh"
c
      integer rtdb  !< [Input] The runtime database handle
      integer geom  !< [Input] The geometry handle
      integer basis !< [Input] The basis set handle
      integer nbf   !< [Input] The number of basis functions
      integer nea   !< [Input] The number of alpha electrons
      integer neb   !< [Input] The number of beta electrons
c
      double precision gstate(8*nbf*nbf+2*nbf) !< [Input] The ground state vector
      double precision state(8*nbf*nbf+2*nbf) !< [In/Output] The state vector
      double precision etot                   !< [Output] The total energy
c
      double precision dd       !< Norm of the gradient
      double precision ds       !< Step length from gradient
      double precision step     !< Step size from linesearch
      double precision tol_conv !< Convergence tolerance
      double precision dmxocc   !< Maximum gradient in the occupation
                                !< function gradient
      double precision dmxnat   !< Maximum gradient in the natural
                                !< orbital gradient
      double precision dummy    !< Needed for unused function arguments
c
      integer l_h1,   k_h1   !< Memory handles for 1-electron integrals
      integer l_ov,   k_ov   !< Memory handles for overlap integrals
      integer l_ov12, k_ov12 !< Memory handles for sqrt(overlap)
      integer l_ovm12, k_ovm12 !< Memory handles for 1/sqrt(overlap)
      integer l_eri,  k_eri  !< Memory handles for Coulomb integrals
      integer l_erix, k_erix !< Memory handles for exchange integrals
      integer l_grad, k_grad !< Memory handles for the gradient
      integer l_fock, k_fock !< Memory handles for the Fock matrices
      integer l_stat, k_stat !< Memory handles for the local state
      integer l_m1,   k_m1   !< Scratch memory
      integer l_m2,   k_m2   !< Scratch memory
      integer l_ovl,  k_ovl  !< Scratch memory
c
      integer it     !< Iteration count
      integer maxit  !< Maximum number of iterations
      logical status !< Status variable
c
      double precision levelshift !< The levelshifter
c
      double precision wfn1_nxt_dot, ddot, dlange
      external         wfn1_nxt_dot, ddot, dlange
c
      maxit = 0000
      tol_conv = 1.0d-5
      status = rtdb_get(rtdb,"wfn1:maxit",MT_INT,1,maxit)
      status = rtdb_get(rtdb,"wfn1:thresh",MT_DBL,1,tol_conv)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"h1",l_h1,k_h1))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate h1",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ov",l_ov,k_ov))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate ov",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ov12",l_ov12,k_ov12))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate ov12",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ovm12",l_ovm12,k_ovm12))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate ovm12",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nbf*nbf,"eri",l_eri,k_eri))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate eri",
     &   ma_sizeof(MT_DBL,nbf*nbf*nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nbf*nbf,"erix",l_erix,k_erix))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate erix",
     &   ma_sizeof(MT_DBL,nbf*nbf*nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"grad",
     &   l_grad,k_grad))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate grad",
     &   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"fock",
     &   l_fock,k_fock))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate fock",
     &   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"stat",
     &   l_stat,k_stat))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate stat",
     &   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ovl",l_ovl,k_ovl))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate ovl",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"m1",l_m1,k_m1))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate m1",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"m2",l_m1,k_m2))
     &   call errquit("wfn1_nxt2_min_energy: could not allocate m2",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
c     Now initialize the memory segments and then compute the integrals c
      call dfill(wfn1_nxt_size1(),0.0d0,dbl_mb(k_grad),1)
      call dfill(wfn1_nxt_size1(),0.0d0,dbl_mb(k_fock),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_h1),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_ov),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_ov12),1)
      call dfill(nbf*nbf*nbf*nbf,0.0d0,dbl_mb(k_eri),1)
      call dfill(nbf*nbf*nbf*nbf,0.0d0,dbl_mb(k_erix),1)
c
      levelshift = 0.0d0
c
      call int_init(rtdb,1,basis)
      call wfn1_1e_tv(basis,nbf,dbl_mb(k_h1))
      call wfn1_1e_s(basis,nbf,dbl_mb(k_ov))
      call wfn1_nxt2_ov12(nbf,dbl_mb(k_ov),dbl_mb(k_ov12),
     &                    dbl_mb(k_ovm12))
      call wfn1_2e_rep(basis,nbf,dbl_mb(k_eri),dbl_mb(k_erix))
      call int_terminate
c
c     Calculate the initial energy and gradient
c
      call wfn1_nxt2_energy_gradient(rtdb,geom,nbf,nea,neb,dbl_mb(k_h1),
     &     dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),dbl_mb(k_ov12),
     &     dbl_mb(k_ovm12),state,dbl_mb(k_fock),dbl_mb(k_grad),etot,
     &     levelshift)
cDEBUG
      write(*,*)'*** Initial energy = ',etot
c     write(*,*)'*** Initial gradient'
c     call wfn1_nxt_print_state(dbl_mb(k_grad),nbf)
      write(*,*)'*** Initial fock'
      call wfn1_nxt_print_state(dbl_mb(k_fock),nbf)
cDEBUG
c
c     dd = ddot(wfn1_nxt_size1(),dbl_mb(k_grad),1,dbl_mb(k_grad),1)
c     dd = sqrt(dd)
      dd = dlange('m',1,wfn1_nxt_size1(),dbl_mb(k_grad),1,dummy)
      dmxocc = max(dlange('m',1,nbf*nbf,
     &             dbl_mb(k_grad-1+wfn1_nxt_aocc()),1),
     &             dlange('m',1,nbf*nbf,
     &             dbl_mb(k_grad-1+wfn1_nxt_bocc()),1))
      dmxnat = max(dlange('m',1,nbf*nbf,
     &             dbl_mb(k_grad-1+wfn1_nxt_anat()),1),
     &             dlange('m',1,nbf*nbf,
     &             dbl_mb(k_grad-1+wfn1_nxt_bnat()),1))
c     if (dmxocc.gt.dmxnat) then
        call dfill(nbf*nbf,0.0d0,dbl_mb(k_grad-1+wfn1_nxt_anat()),1)
        call dfill(nbf*nbf,0.0d0,dbl_mb(k_grad-1+wfn1_nxt_bnat()),1)
c     else
c       call dfill(nbf*nbf,0.0d0,dbl_mb(k_grad-1+wfn1_nxt_aocc()),1)
c       call dfill(nbf*nbf,0.0d0,dbl_mb(k_grad-1+wfn1_nxt_bocc()),1)
c     endif
c
c     Now optimize properly
c
      it = 0
      do while (dd.gt.tol_conv.and.it.le.maxit)
        it = it + 1
c       call wfn1_nxt2_linesearch_e(rtdb,geom,nbf,nea,neb,step,
c    &       dbl_mb(k_h1),dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),
c    &       dbl_mb(k_ov12),dbl_mb(k_ovm12),state,dbl_mb(k_fock),
c    &       dbl_mb(k_grad),etot,tol_conv)
c       call dscal(wfn1_nxt_size1(),-step,dbl_mb(k_grad),1)
        call wfn1_nxt2_step_excit(nea,neb,nbf,it,state,dbl_mb(k_fock),
     &       dbl_mb(k_grad),dbl_mb(k_stat))
        call dcopy(wfn1_nxt_size1(),dbl_mb(k_stat),1,state,1)
c
        call wfn1_ovlp_2_states(nbf,nbf,
     &       gstate(wfn1_nxt_aocc()),gstate(wfn1_nxt_anat()),
     &       state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &       dbl_mb(k_ov),dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_ovl))
        call wfn1_nxt2_sort_overlap(nbf,dbl_mb(k_ovl),
     &       state(wfn1_nxt_aocc()))
c
cDEBUG
c      write(*,*)"*** overlap alpha"
c      call hess_hssout(dbl_mb(k_ovl),nbf,nbf,nbf)
cDEBUG
        call wfn1_ovlp_2_states(nbf,nbf,
     &       gstate(wfn1_nxt_bocc()),gstate(wfn1_nxt_bnat()),
     &       state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &       dbl_mb(k_ov),dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_ovl))
        call wfn1_nxt2_sort_overlap(nbf,dbl_mb(k_ovl),
     &       state(wfn1_nxt_bocc()))
cDEBUG
c      write(*,*)"*** overlap beta"
c      call hess_hssout(dbl_mb(k_ovl),nbf,nbf,nbf)
cDEBUG
        call wfn1_nxt2_energy_gradient(rtdb,geom,nbf,nea,neb,
     &       dbl_mb(k_h1),dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),
     &       dbl_mb(k_ov12),dbl_mb(k_ovm12),state,dbl_mb(k_fock),
     &       dbl_mb(k_grad),etot,levelshift)
cDEBUG
c     write(*,*)'*** fock',it
c     call wfn1_nxt_print_state(dbl_mb(k_fock),nbf)
cDEBUG
c       dd = ddot(wfn1_nxt_size1(),dbl_mb(k_grad),1,dbl_mb(k_grad),1)
c       dd = sqrt(dd)
        dd = dlange('m',1,wfn1_nxt_size1(),dbl_mb(k_grad),1,dummy)
        dmxocc = max(dlange('m',1,nbf*nbf,
     &               dbl_mb(k_grad-1+wfn1_nxt_aocc()),1),
     &               dlange('m',1,nbf*nbf,
     &               dbl_mb(k_grad-1+wfn1_nxt_bocc()),1))
        dmxnat = max(dlange('m',1,nbf*nbf,
     &               dbl_mb(k_grad-1+wfn1_nxt_anat()),1),
     &               dlange('m',1,nbf*nbf,
     &               dbl_mb(k_grad-1+wfn1_nxt_bnat()),1))
c       if (dmxocc.gt.dmxnat) then
          call dfill(nbf*nbf,0.0d0,dbl_mb(k_grad-1+wfn1_nxt_anat()),1)
          call dfill(nbf*nbf,0.0d0,dbl_mb(k_grad-1+wfn1_nxt_bnat()),1)
c       else
c         call dfill(nbf*nbf,0.0d0,dbl_mb(k_grad-1+wfn1_nxt_aocc()),1)
c         call dfill(nbf*nbf,0.0d0,dbl_mb(k_grad-1+wfn1_nxt_bocc()),1)
c       endif
cDEBUG
c       if (mod(it,100).eq.0) then
c         write(*,'("wfn1 it: ",i4,f22.14,e16.4,3f12.6)')
c    &          it,etot,dd,step,dmxocc,dmxnat
          write(*,'("wfn1 it: ",i4,f22.14,e16.4,2f12.6)')
     &          it,etot,dd,dmxocc,dmxnat
c       endif
cDEBUG

      enddo
c
c     Canonicalize the occupation functions
c
      call wfn1_nxt2_canon_occ(nbf,nea,neb,state,dbl_mb(k_fock),
     &     dbl_mb(k_stat))
      call dcopy(wfn1_nxt_size1(),dbl_mb(k_stat),1,state,1)
cDEBUG
 10   continue 
      write(*,*)'*** last energy = ',etot
      call wfn1_nxt_print_occ(nea,neb,state,nbf,state(wfn1_nxt_aoccl()),
     &                        state(wfn1_nxt_boccl()))
      write(*,*)'*** last Fock'
      call wfn1_nxt_print_state(dbl_mb(k_fock),nbf)
      write(*,*)'*** last gradient'
      call wfn1_nxt_print_state(dbl_mb(k_grad),nbf)
      write(*,*)'*** last state'
      call wfn1_nxt_print_state(state,nbf)
      write(*,*)'*** last state overlaps'
      call wfn1_nxt_print_overlap(nbf,state,dbl_mb(k_ov))
cDEBUG
c
c     Clear up the memory
c
      if (.not.ma_chop_stack(l_h1)) then
        call errquit("wfn1_nxt2_min_energy: could free memory",0,MA_ERR)
      endif
      end
C>
C> @}
