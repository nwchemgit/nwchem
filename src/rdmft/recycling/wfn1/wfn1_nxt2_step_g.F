C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Create a new state by moving one step from the old one
C>
C> The new state is generated from the old one by applying appropriate
C> rotations. In this case (unlike the wfn1_nxt2_step case) the 
C> rotations are determined entirely by the gradient, instead of the
C> Fock matrices.
C>
      subroutine wfn1_nxt2_step_g(nea,neb,nbf,it,statein,fock,grad,
     &                            stateout)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "wfn1_nxt.fh"
c
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
      integer nbf !< [Input] The number of basis functions
      integer it  !< [Input] The iteration number
c
      double precision statein(8*nbf*nbf) !< [Input] The initial state
      double precision fock(8*nbf*nbf)    !< [Input] The fock matrices
      double precision grad(8*nbf*nbf)    !< [Input] The gradient
c
      double precision stateout(8*nbf*nbf) !< [Input] The final state
c
      integer info  !< The status flag for DSYEV
      integer lwork !< The work space size for DSYEV
c
      integer l_m1,  k_m1
      integer l_m2,  k_m2
      integer l_exp, k_exp
      integer ir
c
      double precision levelshift1
      double precision levelshift2
      double precision levelshift3
      double precision levelshift4
      parameter(levelshift1 = 2.0d0)
      parameter(levelshift2 = 1.0d0)
      parameter(levelshift3 = 0.5d0)
      parameter(levelshift4 = 0.5d0)
      double precision levelshift
c
c     Allocate memory
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"tmp1",l_m1,k_m1))
     &  call errquit("wfn1_nxt2_step: could not allocate tmp1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"tmp2",l_m2,k_m2))
     &  call errquit("wfn1_nxt2_step: could not allocate tmp2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"exp",l_exp,k_exp))
     &  call errquit("wfn1_nxt2_step: could not allocate exp",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
c     Select the occupation function level shift
c
c     if (it.le.20) then
c       levelshift = levelshift1
c     else if (it.le.40) then
c       levelshift = levelshift2
c     else if (it.le.60) then
c       levelshift = levelshift3
c     else
c       levelshift = levelshift4
c     endif
c
c     Initialize the output state (needed for the electron
c     configuration)
c
      call dcopy(wfn1_nxt_size1(),statein,1,stateout,1)
c
c     Alpha occupation functions
c
      call wfn1_exp(nbf,grad(wfn1_nxt_aocc()),dbl_mb(k_m1),dbl_mb(k_m2),
     &              dbl_mb(k_exp))
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_aocc()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Beta occupation functions
c
      call wfn1_exp(nbf,grad(wfn1_nxt_bocc()),dbl_mb(k_m1),dbl_mb(k_m2),
     &              dbl_mb(k_exp))
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bocc()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c
c     Alpha natural orbitals
c
      call wfn1_exp(nbf,grad(wfn1_nxt_anat()),dbl_mb(k_m1),dbl_mb(k_m2),
     &              dbl_mb(k_exp))
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_anat()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_anat()),nbf)
      call dcopy(nbf*nbf,stateout(wfn1_nxt_aocc()),1,dbl_mb(k_m1),1)
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
     &           dbl_mb(k_m1),nbf,
     &           0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Beta natural orbitals
c
      call wfn1_exp(nbf,grad(wfn1_nxt_bnat()),dbl_mb(k_m1),dbl_mb(k_m2),
     &              dbl_mb(k_exp))
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bnat()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bnat()),nbf)
      call dcopy(nbf*nbf,stateout(wfn1_nxt_bocc()),1,dbl_mb(k_m1),1)
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
     &           dbl_mb(k_m1),nbf,
     &           0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c
c     Cleanup memory
c
      if (.not.ma_pop_stack(l_exp))
     &  call errquit("wfn1_nxt2_step: could not deallocate exp",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_m2))
     &  call errquit("wfn1_nxt2_step: could not deallocate tmp2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_m1))
     &  call errquit("wfn1_nxt2_step: could not deallocate tmp1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      end
C> @}
