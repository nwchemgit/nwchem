C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Create a new state by moving one step from the old one
C>
C> The new state is generated from the old one by applying appropriate
C> rotations.
C>
      subroutine wfn1_nxt2_step_excit(nea,neb,nbf,it,statein,fock,grad,
     &                                stateout)
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
      double precision statein(8*nbf*nbf+2*nbf) !< [Input] The initial state
      double precision fock(8*nbf*nbf+2*nbf)    !< [Input] The fock matrices
      double precision grad(8*nbf*nbf+2*nbf)    !< [Input] The gradient
c
      double precision stateout(8*nbf*nbf+2*nbf) !< [Input] The final state
c
      integer info  !< The status flag for DSYEV
      integer lwork !< The work space size for DSYEV
c
      integer l_m1,  k_m1
      integer l_m2,  k_m2
      integer l_exp, k_exp
      integer ir
c
      double precision elow
      double precision ehigh
c
      double precision levelshift1
      double precision levelshift2
      double precision levelshift3
      double precision levelshift4
      parameter(levelshift1 = 2.0d0)
      parameter(levelshift2 = 1.0d0)
      parameter(levelshift3 = 0.5d0)
      parameter(levelshift4 = 1.0d-5)
      double precision levelshift
c
c     Allocate memory
c
      lwork = max(nbf*nbf,3*nbf)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"tmp1",l_m1,k_m1))
     &  call errquit("wfn1_nxt2_step: could not allocate tmp1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,lwork,"tmp2",l_m2,k_m2))
     &  call errquit("wfn1_nxt2_step: could not allocate tmp2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"exp",l_exp,k_exp))
     &  call errquit("wfn1_nxt2_step: could not allocate exp",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
c     Select the occupation function level shift
c
      if (it.le.20) then
        levelshift = levelshift1
      else if (it.le.40) then
        levelshift = levelshift2
c     else if (it.le.60) then
        levelshift = levelshift3
c     else
c       levelshift = levelshift4
      endif
c     levelshift = levelshift3
c
c     Alpha occupation functions
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_aocc()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,fock(wfn1_nxt_aocc()),1,dbl_mb(k_exp),1)
      do ir = 2, nbf
c
c       Using a tiny levelshifter to fix the order of the degenerate
c       occupation functions
c
        elow = dbl_mb(k_exp+(ir-2)+(ir-2)*nbf) + levelshift4
        if (statein(wfn1_nxt_aconf()+(ir-1)).le.0.0d0) then
          elow = elow + levelshift
        endif
        ehigh = dbl_mb(k_exp+(ir-1)+(ir-1)*nbf)
        dbl_mb(k_exp+(ir-1)+(ir-1)*nbf) = max(elow,ehigh)
      enddo
      info = 0 
c     lwork = nbf*nbf
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_step: dsyev failed a: ",info,UERR)
      endif
cDEBUG
c     write(*,*)"wfn1_nxt2_step: alpha occupation functions"
c     do ir = 1, nbf
c       write(*,*)'eig = ',ir,dbl_mb(k_m1-1+ir)
c     enddo
cDEBUG
c     call wfn1_nxt2_sort_columns(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      if (mod(it,40).ne.39) then
        call dcopy(nbf*nbf,statein(wfn1_nxt_aocc()),1,
     &                     stateout(wfn1_nxt_aocc()),1)
      else
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,
     &             statein(wfn1_nxt_aocc()),nbf,
     &             dbl_mb(k_exp),nbf,0.0d0,
     &             stateout(wfn1_nxt_aocc()),nbf)
      endif
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_aocc()),nbf,
c    &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Beta occupation functions
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_bocc()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,fock(wfn1_nxt_bocc()),1,dbl_mb(k_exp),1)
      do ir = 2, nbf
c
c       Using a tiny levelshifter to fix the order of the degenerate
c       occupation functions
c
        elow = dbl_mb(k_exp+(ir-2)+(ir-2)*nbf) + levelshift4
        if (statein(wfn1_nxt_bconf()+(ir-1)).le.0.0d0) then
          elow = elow + levelshift
        endif
        ehigh = dbl_mb(k_exp+(ir-1)+(ir-1)*nbf)
        dbl_mb(k_exp+(ir-1)+(ir-1)*nbf) = max(elow,ehigh)
      enddo
      info = 0 
c     lwork = nbf*nbf
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_step: dsyev failed b: ",info,UERR)
      endif
cDEBUG
c     write(*,*)"wfn1_nxt2_step: beta occupation functions"
c     do ir = 1, nbf
c       write(*,*)'eig = ',ir,dbl_mb(k_m1-1+ir)
c     enddo
cDEBUG
c     call wfn1_nxt2_sort_columns(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      if (mod(it,40).ne.39) then
        call dcopy(nbf*nbf,statein(wfn1_nxt_bocc()),1,
     &                     stateout(wfn1_nxt_bocc()),1)
      else
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,
     &             statein(wfn1_nxt_bocc()),nbf,
     &             dbl_mb(k_exp),nbf,0.0d0,
     &             stateout(wfn1_nxt_bocc()),nbf)
      endif
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bocc()),nbf,
c    &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c
c     Alpha natural orbitals
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_anat()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,fock(wfn1_nxt_anat()),1,dbl_mb(k_exp),1)
      call wfn1_nxt2_sort_diag(nbf,dbl_mb(k_exp),
     &     statein(wfn1_nxt_anat()),stateout(wfn1_nxt_aocc()))
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call wfn1_nxt2_shift(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      info = 0 
c     lwork = nbf*nbf
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_step: dsyev failed a: ",info,UERR)
      endif
cDEBUG
c     write(*,*)"wfn1_nxt2_step: alpha natural orbitals"
cDEBUG
c     call wfn1_nxt2_sort_columns(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      if (mod(it,40).eq.39) then
        call dcopy(nbf*nbf,statein(wfn1_nxt_anat()),1,
     &                     stateout(wfn1_nxt_anat()),1)
      else
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,
     &             statein(wfn1_nxt_anat()),nbf,
     &             dbl_mb(k_exp),nbf,0.0d0,
     &             stateout(wfn1_nxt_anat()),nbf)
      endif
c     call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
c    &           statein(wfn1_nxt_aocc()),nbf,
c    &           0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Beta natural orbitals
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_bnat()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,fock(wfn1_nxt_bnat()),1,dbl_mb(k_exp),1)
      call wfn1_nxt2_sort_diag(nbf,dbl_mb(k_exp),
     &     statein(wfn1_nxt_bnat()),stateout(wfn1_nxt_bocc()))
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call wfn1_nxt2_shift(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      info = 0 
c     lwork = nbf*nbf
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_step: dsyev failed b: ",info,UERR)
      endif
cDEBUG
c     write(*,*)"wfn1_nxt2_step: beta natural orbitals"
cDEBUG
c     call wfn1_nxt2_sort_columns(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      if (mod(it,40).eq.39) then
        call dcopy(nbf*nbf,statein(wfn1_nxt_bnat()),1,
     &                     stateout(wfn1_nxt_bnat()),1)
      else
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,
     &             statein(wfn1_nxt_bnat()),nbf,
     &             dbl_mb(k_exp),nbf,0.0d0,
     &             stateout(wfn1_nxt_bnat()),nbf)
      endif
c     call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
c    &           statein(wfn1_nxt_bocc()),nbf,
c    &           0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c
c     Copy the electron configurations
c
      call dcopy(nbf,statein(wfn1_nxt_aconf()),1,
     &               stateout(wfn1_nxt_aconf()),1)
      call dcopy(nbf,statein(wfn1_nxt_bconf()),1,
     &               stateout(wfn1_nxt_bconf()),1)
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
