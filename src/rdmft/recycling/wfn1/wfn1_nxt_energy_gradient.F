C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Compute the energy and gradient of the current state vector
C>
C> Computes the energy and gradient of the state vector given the
C> various sets of integrals. This routine can generate various
C> energy expressions including:
C> 1. Hartree-Fock
C> 2. Density Functional Theory
C> 3. Density Matrix Functional Theory
C> 4. Entropy model
C>
      subroutine wfn1_nxt_energy_gradient(rtdb,geom,nbf,nea,neb,h1,eri,
     &           erix,ov,ov12,state,grad,Etot)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "wfn1_nxt.fh"
#include "geom.fh"
c
      integer rtdb !< [Input] The RTDB handle (needed for the DFT terms)
      integer geom !< [Input] The geometry handle (needed for the
                   !< nuclear repulsion energy)
      integer nbf  !< [Input] The number of basis functions
      integer nea  !< [Input] The number of alpha electrons
      integer neb  !< [Input] The number of beta electrons
c
      double precision h1(nbf,nbf)           !< [Input] The 1-electron
                                             !< integrals
      double precision eri(nbf,nbf,nbf,nbf)  !< [Input] The Coulomb
                                             !< integrals
      double precision erix(nbf,nbf,nbf,nbf) !< [Input] The exchange
                                             !< integrals
      double precision ov(nbf,nbf)           !< [Input] The overlap
                                             !< integrals
      double precision ov12(nbf,nbf)         !< [Input] The 1/sqrt(overlap)
c
      double precision state(8*nbf*nbf) !< [Input] The state vector
      double precision grad(8*nbf*nbf)  !< [Output] The gradient
cDEBUG
c     double precision gral(8*nbf*nbf)  !< [Scratch] The gradient
cDEBUG
c
      double precision Etot !< [Output] The total energy
c
      integer l_occa, k_occa !< Memory for alpha electron occupation
                             !< numbers
      integer l_occb, k_occb !< Memory for beta electron occupation
                             !< numbers
c
      integer l_da, k_da !< Memory for the alpha electron density
                         !< matrix
      integer l_db, k_db !< Memory for the beta electron density
                         !< matrix
c
      integer l_fa, k_fa !< Memory for the alpha electron Fock
                         !< matrix
      integer l_fb, k_fb !< Memory for the beta electron Fock
                         !< matrix
c
      integer l_ta, k_ta !< Scratch matrix
      integer l_tb, k_tb !< Scratch matrix
c
      integer k_tv !< Scratch vector
c
      double precision Enuc   !< The nuclear repulsion energy
      double precision E1el   !< The 1-electron energy
      double precision E2elC  !< The 2-electron Coulomb energy
      double precision E2elX  !< The 2-electron eXchange energy
      double precision ELo    !< The occupation function Lagrangian
      double precision ELn    !< The natural orbital Lagrangian
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Da',l_da,k_da))
     &   call errquit("wfn1_nxt_energy_gradient: could not allocate "
     &        //"Da",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Db',l_db,k_db))
     &   call errquit("wfn1_nxt_energy_gradient: could not allocate "
     &        //"Db",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Fa',l_fa,k_fa))
     &   call errquit("wfn1_nxt_energy_gradient: could not allocate "
     &        //"Fa",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Fb',l_fb,k_fb))
     &   call errquit("wfn1_nxt_energy_gradient: could not allocate "
     &        //"Fb",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf,'occ a',l_occa,k_occa))
     &   call errquit("wfn1_nxt_energy_gradient: could not allocate "
     &        //"occa",ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,'occ b',l_occb,k_occb))
     &   call errquit("wfn1_nxt_energy_gradient: could not allocate "
     &        //"occb",ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      k_tv = k_occa
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Ta',l_ta,k_ta))
     &   call errquit("wfn1_nxt_energy_gradient: could not allocate "
     &        //"Ta",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Tb',l_tb,k_tb))
     &   call errquit("wfn1_nxt_energy_gradient: could not allocate "
     &        //"Tb",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
c     Zero the gradient
c
      call dfill(wfn1_nxt_size1(),0.0d0,grad,1)
c
c     Get the nuclear repulsion energy
c
      if (.not.geom_nuc_rep_energy(geom,Enuc))
     &  call errquit("wfn1_energy: Enuc failed",0,GEOM_ERR)
c
c     Compute the 1-electron energy and gradient
c
cDEBUG
c     call dfill(wfn1_nxt_size1(),0.0d0,gral,1)
cDEBUG
      call wfn1_nxt_1el_eg(nbf,nea,neb,h1,
     &     dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_da),dbl_mb(k_db),
     &     state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &     state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &     grad(wfn1_nxt_aocc()),grad(wfn1_nxt_anat()),
     &     grad(wfn1_nxt_bocc()),grad(wfn1_nxt_bnat()),E1el)
cDEBUG
c     write(*,*)'*** 1-electron energy'
c     call wfn1_nxt_print_state(gral,nbf)
c     call daxpy(wfn1_nxt_size1(),1.0d0,gral,1,grad,1)
c     call dfill(wfn1_nxt_size1(),0.0d0,gral,1)
cDEBUG
      call wfn1_nxt_2el_C_eg(nbf,nea,neb,eri,
     &     dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_da),dbl_mb(k_db),
     &     state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &     state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &     grad(wfn1_nxt_aocc()),grad(wfn1_nxt_anat()),
     &     grad(wfn1_nxt_bocc()),grad(wfn1_nxt_bnat()),E2elC,
     &     dbl_mb(k_fa),dbl_mb(k_fb))
cDEBUG
c     write(*,*)'*** 2-electron Coulomb energy'
c     call wfn1_nxt_print_state(gral,nbf)
c     call daxpy(wfn1_nxt_size1(),1.0d0,gral,1,grad,1)
c     call dfill(wfn1_nxt_size1(),0.0d0,gral,1)
cDEBUG
      call wfn1_nxt_2el_X_eg(nbf,nea,neb,wfn1_nxt_x_pow,erix,
     &     dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_da),dbl_mb(k_db),
     &     state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &     state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &     grad(wfn1_nxt_aocc()),grad(wfn1_nxt_anat()),
     &     grad(wfn1_nxt_bocc()),grad(wfn1_nxt_bnat()),E2elX,
     &     dbl_mb(k_fa),dbl_mb(k_fb))
cDEBUG
c     write(*,*)'*** 2-electron eXchange energy'
c     call wfn1_nxt_print_state(gral,nbf)
c     call daxpy(wfn1_nxt_size1(),1.0d0,gral,1,grad,1)
c     call dfill(wfn1_nxt_size1(),0.0d0,gral,1)
cDEBUG
c
c     call wfn1_nxt_Lo_eg(nbf,
c    &     state(wfn1_nxt_aocc()),state(wfn1_nxt_bocc()),
c    &     state(wfn1_nxt_aoccl()),state(wfn1_nxt_boccl()),
c    &     grad(wfn1_nxt_aocc()),grad(wfn1_nxt_bocc()),
c    &     grad(wfn1_nxt_aoccl()),grad(wfn1_nxt_boccl()),ELo)
c     call wfn1_nxt_Ln_eg(nbf,ov,
c    &     state(wfn1_nxt_anat()),state(wfn1_nxt_bnat()),
c    &     state(wfn1_nxt_anatl()),state(wfn1_nxt_bnatl()),
c    &     grad(wfn1_nxt_anat()),grad(wfn1_nxt_bnat()),
c    &     grad(wfn1_nxt_anatl()),grad(wfn1_nxt_bnatl()),ELn,
c    &     dbl_mb(k_tv))
c
      call wfn1_nxt_Po_eg(nbf,
     &     state(wfn1_nxt_aocc()),state(wfn1_nxt_bocc()),
     &     state(wfn1_nxt_aoccl()),state(wfn1_nxt_boccl()),
     &     grad(wfn1_nxt_aocc()),grad(wfn1_nxt_bocc()),
     &     grad(wfn1_nxt_aoccl()),grad(wfn1_nxt_boccl()),ELo,
     &     dbl_mb(k_ta),dbl_mb(k_tb))
cDEBUG
c     write(*,*)'*** occupation function penalty'
c     call wfn1_nxt_print_state(gral,nbf)
c     call daxpy(wfn1_nxt_size1(),1.0d0,gral,1,grad,1)
c     call dfill(wfn1_nxt_size1(),0.0d0,gral,1)
cDEBUG
      call wfn1_nxt_Pn_eg(nbf,ov,
     &     state(wfn1_nxt_anat()),state(wfn1_nxt_bnat()),
     &     state(wfn1_nxt_anatl()),state(wfn1_nxt_bnatl()),
     &     grad(wfn1_nxt_anat()),grad(wfn1_nxt_bnat()),
     &     grad(wfn1_nxt_anatl()),grad(wfn1_nxt_bnatl()),ELn,
     &     dbl_mb(k_tv),
     &     dbl_mb(k_ta),dbl_mb(k_tb))
cDEBUG
c     write(*,*)'*** natural orbital penalty'
c     call wfn1_nxt_print_state(gral,nbf)
c     call daxpy(wfn1_nxt_size1(),1.0d0,gral,1,grad,1)
c     call dfill(wfn1_nxt_size1(),0.0d0,gral,1)
cDEBUG
c
      if (.not.ma_chop_stack(l_da))
     &   call errquit("wfn1_nxt_energy_gradient: could not deallocate "
     &        //"Da",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      Etot = Enuc + E1el + E2elC + E2elX + ELo + ELn
cDEBUG
      write(*,*)'Enuc  = ',Enuc
      write(*,*)'E1el  = ',E1el
      write(*,*)'E2elC = ',E2elC
      write(*,*)'E2elX = ',E2elX
      write(*,*)'EPo   = ',ELo
      write(*,*)'EPn   = ',ELn
      write(*,*)'Etot  = ',Etot
cDEBUG
c
      end
C>
C> @}
