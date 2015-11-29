C> \ingroup wfn1_nxt2
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
      subroutine wfn1_nxt2_energy_gradient(rtdb,geom,nbf,nea,neb,h1,eri,
     &           erix,ov,ov12,ovm12,state,fock,grad,Etot,levelshift)
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
      double precision ov12(nbf,nbf)         !< [Input] The sqrt(overlap)
      double precision ovm12(nbf,nbf)        !< [Input] The 1/sqrt(overlap)
c
      double precision state(8*nbf*nbf+2*nbf) !< [Input] The state vector
      double precision fock(8*nbf*nbf+2*nbf)  !< [Output] The Fock matrices
      double precision grad(8*nbf*nbf+2*nbf)  !< [Output] The gradient
c
      double precision Etot !< [Output] The total energy
      double precision levelshift !< [In/Output] The levelshifter
c
      integer l_occa, k_occa !< Memory for alpha electron occupation
                             !< numbers
      integer l_occb, k_occb !< Memory for beta electron occupation
                             !< numbers
c
      integer l_dfa, k_dfa   !< Memory for alpha electron occupation
                             !< numbers
      integer l_dfb, k_dfb   !< Memory for beta electron occupation
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
      integer l_t1, k_t1 !< Scratch matrix
      integer l_t2, k_t2 !< Scratch matrix
      integer l_t4, k_t4 !< Scratch matrix
c
      integer k_tv !< Scratch vector
c
      double precision Enuc   !< The nuclear repulsion energy
      double precision E1el   !< The 1-electron energy
      double precision E2elC  !< The 2-electron Coulomb energy
      double precision E2elX  !< The 2-electron eXchange energy
      double precision Ef1el  !< The 1-electron energy
      double precision Ef2elC !< The 2-electron Coulomb energy
      double precision Ef2elX !< The 2-electron eXchange energy
      double precision ELo    !< The occupation function Lagrangian
      double precision ELn    !< The natural orbital Lagrangian
      double precision Ef     !< The Fock energy
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Da',l_da,k_da))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"Da",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Db',l_db,k_db))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"Db",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Fa',l_fa,k_fa))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"Fa",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Fb',l_fb,k_fb))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"Fb",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf,'occ a',l_occa,k_occa))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"occa",ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,'occ b',l_occb,k_occb))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"occb",ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      k_tv = k_occa
c
      if (.not.ma_push_get(MT_DBL,nbf,'dfa',l_dfa,k_dfa))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"dfa",ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,'dfb',l_dfb,k_dfb))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"dfb",ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Ta',l_ta,k_ta))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"Ta",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'Tb',l_tb,k_tb))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"Tb",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'T1',l_t1,k_t1))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"T1",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'T2',l_t2,k_t2))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"T2",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf**4,'T4',l_t4,k_t4))
     &   call errquit("wfn1_nxt2_energy_gradient: could not allocate "
     &        //"T4",ma_sizeof(MT_DBL,nbf**4,MT_BYTE),MA_ERR)
c
c     Zero the Fock matrices and gradients
c
      call dfill(wfn1_nxt_size1(),0.0d0,fock,1)
      call dfill(wfn1_nxt_size1(),0.0d0,grad,1)
c
c     Zero the Fock matrices
c
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_fa),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_fb),1)
c
c     Get the nuclear repulsion energy
c
      if (.not.geom_nuc_rep_energy(geom,Enuc))
     &  call errquit("wfn1_energy: Enuc failed",0,GEOM_ERR)
c
c     Compute the 1-electron energy and gradient
c
      if (wfn1_nxt_ad) then
c
c       Use automatic differentiation for matrix function derivatives
c       (the Fock matrix is calculated in AO basis)
c
        call wfn1_nxt2_1el_ef(nbf,nea,neb,h1,
     &       dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_da),dbl_mb(k_db),
     &       state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &       state(wfn1_nxt_aconf()),
     &       state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &       state(wfn1_nxt_bconf()),
     &       dbl_mb(k_fa),dbl_mb(k_fb),E1el)
        call wfn1_nxt2_2el_C_ef(nbf,nea,neb,eri,
     &       dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_da),dbl_mb(k_db),
     &       state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &       state(wfn1_nxt_aconf()),
     &       state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &       state(wfn1_nxt_bconf()),
     &       dbl_mb(k_fa),dbl_mb(k_fb),E2elC,
     &       dbl_mb(k_ta),dbl_mb(k_tb))
        call wfn1_nxt2_2el_X_ef(nbf,nea,neb,wfn1_nxt_x_pow,erix,ov,ov12,
     &       ovm12,dbl_mb(k_occa),dbl_mb(k_occb),
     &       dbl_mb(k_da),dbl_mb(k_db),
     &       state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &       state(wfn1_nxt_aconf()),
     &       state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &       state(wfn1_nxt_bconf()),
     &       dbl_mb(k_fa),dbl_mb(k_fb),E2elX,
     &       dbl_mb(k_ta),dbl_mb(k_tb),dbl_mb(k_t4))
      else
c
c       Use analytic differentiation for matrix function derivatives
c       (the Fock matrix is calculated in natural orbital basis)
c
        call wfn1_nxt2_1el_n_ef(nbf,nea,neb,h1,
     &       dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_da),dbl_mb(k_db),
     &       state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &       state(wfn1_nxt_aconf()),
     &       state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &       state(wfn1_nxt_bconf()),
     &       dbl_mb(k_fa),dbl_mb(k_fb),E1el,Ef1el,
     &       dbl_mb(k_ta),dbl_mb(k_tb))
        call wfn1_nxt2_2el_n_C_ef(nbf,nea,neb,eri,
     &       dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_da),dbl_mb(k_db),
     &       state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &       state(wfn1_nxt_aconf()),
     &       state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &       state(wfn1_nxt_bconf()),
     &       dbl_mb(k_fa),dbl_mb(k_fb),E2elC,Ef2elC,
     &       dbl_mb(k_ta),dbl_mb(k_tb),dbl_mb(k_t1),dbl_mb(k_t2))
        call wfn1_nxt2_2el_n_X_ef(nbf,nea,neb,wfn1_nxt_x_pow,erix,
     &       ov,ov12,ovm12,dbl_mb(k_occa),dbl_mb(k_occb),
     &       dbl_mb(k_dfa),dbl_mb(k_dfb),
     &       dbl_mb(k_da),dbl_mb(k_db),
     &       state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &       state(wfn1_nxt_aconf()),
     &       state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &       state(wfn1_nxt_bconf()),
     &       dbl_mb(k_fa),dbl_mb(k_fb),E2elX,Ef2elX,
     &       dbl_mb(k_ta),dbl_mb(k_tb),dbl_mb(k_t1),dbl_mb(k_t2),
     &       dbl_mb(k_t4))
      endif
c
      if (.not.ma_chop_stack(l_da))
     &   call errquit("wfn1_nxt_energy_gradient: could not deallocate "
     &        //"Da",ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      Etot   = Enuc + E1el + E2elC + E2elX
      Ef2elC = Ef2elC*0.5d0
      Ef2elX = Ef2elX*0.5d0/wfn1_nxt_x_pow
cDEBUG
c     write(*,*)
c     write(*,*)'Enuc   = ',Enuc
c     write(*,*)'E1el   = ',E1el
c     write(*,*)'E2elC  = ',E2elC
c     write(*,*)'E2elX  = ',E2elX
c     write(*,*)'Etot   = ',Etot
c     write(*,*)
c     write(*,*)'Enuc   = ',Enuc
c     write(*,*)'Ef1el  = ',Ef1el
c     write(*,*)'Ef2elC = ',Ef2elC
c     write(*,*)'Ef2elX = ',Ef2elX
c     write(*,*)'Etot   = ',Enuc + Ef1el + Ef2elC + Ef2elX
c     write(*,*)'*** after wfn1_nxt2_2el_X_ef'
c     call hess_hssout(h1,nbf,nbf,nbf)
c     call hess_hssout(dbl_mb(k_fa),nbf,nbf,nbf)
c     call hess_hssout(dbl_mb(k_fb),nbf,nbf,nbf)
cDEBUG
c
c     Construct the gradients
c
      if (wfn1_nxt_ad) then
c
c       The input Fock matrix is given in AO basis
c
        call wfn1_nxt2_fg(nbf,nea,
     &       state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &       state(wfn1_nxt_aconf()),dbl_mb(k_fa),
     &       fock(wfn1_nxt_aocc()),fock(wfn1_nxt_anat()),
     &       grad(wfn1_nxt_aocc()),grad(wfn1_nxt_anat()),
     &       dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_ta))
        call wfn1_nxt2_fg(nbf,neb,
     &       state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &       state(wfn1_nxt_bconf()),dbl_mb(k_fb),
     &       fock(wfn1_nxt_bocc()),fock(wfn1_nxt_bnat()),
     &       grad(wfn1_nxt_bocc()),grad(wfn1_nxt_bnat()),
     &       dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_tb))
      else 
c
c       The input Fock matrix is given in natural orbital basis
c
        Ef = 0.0d0
cDEBUG
c       write(*,*)'*** en_grad: alpha Fock'
cDEBUG
        call wfn1_nxt2_n_fg(nbf,nea,
     &       state(wfn1_nxt_aocc()),state(wfn1_nxt_anat()),
     &       state(wfn1_nxt_aconf()),dbl_mb(k_fa),
     &       fock(wfn1_nxt_aocc()),fock(wfn1_nxt_anat()),
     &       grad(wfn1_nxt_aocc()),grad(wfn1_nxt_anat()),
     &       dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_ta),Ef,levelshift)
cDEBUG
c       write(*,*)'*** en_grad: beta Fock'
cDEBUG
        call wfn1_nxt2_n_fg(nbf,neb,
     &       state(wfn1_nxt_bocc()),state(wfn1_nxt_bnat()),
     &       state(wfn1_nxt_bconf()),dbl_mb(k_fb),
     &       fock(wfn1_nxt_bocc()),fock(wfn1_nxt_bnat()),
     &       grad(wfn1_nxt_bocc()),grad(wfn1_nxt_bnat()),
     &       dbl_mb(k_occa),dbl_mb(k_occb),dbl_mb(k_tb),Ef,levelshift)
cDEBUG
c     write(*,*)
c     write(*,*)'Enuc  = ',Enuc
c     write(*,*)'Ef    = ',Ef
c     write(*,*)'Etot  = ',Ef+Enuc+E1el/2.0d0
c     write(*,*)
cDEBUG
      endif
c
      call wfn1_nxt2_sym_2_skw(nbf,grad(wfn1_nxt_aocc()))
      call wfn1_nxt2_sym_2_skw(nbf,grad(wfn1_nxt_bocc()))
      call wfn1_nxt2_sym_2_skw(nbf,grad(wfn1_nxt_anat()))
      call wfn1_nxt2_sym_2_skw(nbf,grad(wfn1_nxt_bnat()))
c
      end
C>
C> @}
