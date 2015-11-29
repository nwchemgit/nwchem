C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Create a new state by moving one step from the old one
C> using matrix diagonalization
C>
C> The new state is generated from the old one by applying appropriate
C> rotations. The rotations are obtained from matrix diagonalization.
C>
      subroutine wfn1_nxt2_step_dg(nea,neb,nbf,it,statein,fock,grad,
     &                             stateout,levelshift)
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
      integer l_occa, k_occa
      integer l_occb, k_occb
      integer ir, ii
c
      double precision elow
      double precision ehigh
c
      double precision levelshift1
      double precision levelshift2
      double precision levelshift3
      double precision levelshift4
      double precision levelshift5
      parameter(levelshift1 = 4.0d0)
      parameter(levelshift2 = 2.0d0)
      parameter(levelshift3 = 1.0d0)
      parameter(levelshift4 = 0.5d0)
      parameter(levelshift5 = 1.0d-5)
c     parameter(levelshift5 = 1.0d-2)
      double precision levelshift !< [Input] The levelshifter
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
      if (.not.ma_push_get(MT_DBL,nbf,"occa",l_occa,k_occa))
     &  call errquit("wfn1_nxt2_step: could not allocate occa",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,"occb",l_occb,k_occb))
     &  call errquit("wfn1_nxt2_step: could not allocate occb",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
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
c     Compute the occupation numbers
c
      call dfill(nbf,0.0d0,dbl_mb(k_occa),1)
      call dfill(nbf,0.0d0,dbl_mb(k_occb),1)
      do ir = 1, nbf
        do ii = 1, nbf
          dbl_mb(k_occa+ii-1) = dbl_mb(k_occa+ii-1)
     &       + statein(wfn1_nxt_aocc()+(ii-1)+nbf*(ir-1))**2
     &       * statein(wfn1_nxt_aconf()+ir-1)
          dbl_mb(k_occb+ii-1) = dbl_mb(k_occb+ii-1)
     &       + statein(wfn1_nxt_bocc()+(ii-1)+nbf*(ir-1))**2
     &       * statein(wfn1_nxt_bconf()+ir-1)
        enddo
      enddo

c
c     Alpha occupation functions
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_aocc()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,fock(wfn1_nxt_aocc()),1,dbl_mb(k_exp),1)
cDEBUG
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_m1),levelshift,
     &                      statein(wfn1_nxt_aconf()))
cDEBUG
c     do ir = 2, nbf
c
c       Using a tiny levelshifter to fix the order of the degenerate
c       occupation functions
c
c       elow = dbl_mb(k_exp+(ir-2)+(ir-2)*nbf) + levelshift5
c       if (statein(wfn1_nxt_aconf()+(ir-1)).le.0.0d0) then
c         elow = elow + levelshift
c       endif
c       ehigh = dbl_mb(k_exp+(ir-1)+(ir-1)*nbf)
c       dbl_mb(k_exp+(ir-1)+(ir-1)*nbf) = max(elow,ehigh)
c     enddo
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
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_aocc()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_aocc()),nbf,
c    &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Beta occupation functions
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_bocc()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,fock(wfn1_nxt_bocc()),1,dbl_mb(k_exp),1)
cDEBUG
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_m1),levelshift,
     &                      statein(wfn1_nxt_bconf()))
cDEBUG
c     do ir = 2, nbf
c
c       Using a tiny levelshifter to fix the order of the degenerate
c       occupation functions
c
c       elow = dbl_mb(k_exp+(ir-2)+(ir-2)*nbf) + levelshift5
c       if (statein(wfn1_nxt_bconf()+(ir-1)).le.0.0d0) then
c         elow = elow + levelshift
c       endif
c       ehigh = dbl_mb(k_exp+(ir-1)+(ir-1)*nbf)
c       dbl_mb(k_exp+(ir-1)+(ir-1)*nbf) = max(elow,ehigh)
c     enddo
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
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bocc()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bocc()),nbf)
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
cDEBUG
c     call dscal(nbf,max(1.0d0,levelshift),dbl_mb(k_m1),1)
cDEBUG
      call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_m1),levelshift,
     &                      dbl_mb(k_occa))
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
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_anat()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_anat()),nbf)
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
cDEBUG
c     call dscal(nbf,max(1.0d0,levelshift),dbl_mb(k_m1),1)
cDEBUG
      call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_m1),levelshift,
     &                      dbl_mb(k_occb))
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
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bnat()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bnat()),nbf)
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
      if (.not.ma_pop_stack(l_occb))
     &  call errquit("wfn1_nxt2_step: could not deallocate occb",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_occa))
     &  call errquit("wfn1_nxt2_step: could not deallocate occa",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
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
C>
C> \brief Create a new state by moving one step from the old one
C> using matrix exponents 
C>
C> The new state is generated from the old one by applying appropriate
C> rotations. The rotations are obtained from calculating the
C> exponent of the gradient matrix. This approach is used with 
C> linesearches (ls) hence the name.
C>
      subroutine wfn1_nxt2_step_ls(nea,neb,nbf,it,statein,fock,grad,
     &                             stateout,levelshift)
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
      integer l_occa, k_occa
      integer l_occb, k_occb
      integer ir, ii
c
      double precision elow
      double precision ehigh
c
      double precision levelshift !< [Input] The levelshifter
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
      if (.not.ma_push_get(MT_DBL,nbf,"occa",l_occa,k_occa))
     &  call errquit("wfn1_nxt2_step: could not allocate occa",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,"occb",l_occb,k_occb))
     &  call errquit("wfn1_nxt2_step: could not allocate occb",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
c
c     Compute the occupation numbers
c
c     call dfill(nbf,0.0d0,dbl_mb(k_occa),1)
c     call dfill(nbf,0.0d0,dbl_mb(k_occb),1)
c     do ir = 1, nbf
c       do ii = 1, nbf
c         dbl_mb(k_occa+ii-1) = dbl_mb(k_occa+ii-1)
c    &       + statein(wfn1_nxt_aocc()+(ii-1)+nbf*(ir-1))**2
c    &       * statein(wfn1_nxt_aconf()+ir-1)
c         dbl_mb(k_occb+ii-1) = dbl_mb(k_occb+ii-1)
c    &       + statein(wfn1_nxt_bocc()+(ii-1)+nbf*(ir-1))**2
c    &       * statein(wfn1_nxt_bconf()+ir-1)
c       enddo
c     enddo
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
cDEBUG
c     also rotate the rows in the occupation functions
      call dcopy(nbf*nbf,stateout(wfn1_nxt_aocc()),1,dbl_mb(k_m2),1)
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
     &           dbl_mb(k_m2),nbf,0.0d0,stateout(wfn1_nxt_aocc()),nbf)
cDEBUG
c
c     Beta natural orbitals
c
      call wfn1_exp(nbf,grad(wfn1_nxt_bnat()),dbl_mb(k_m1),dbl_mb(k_m2),
     &              dbl_mb(k_exp))
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bnat()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bnat()),nbf)
cDEBUG
c     also rotate the rows in the occupation functions
      call dcopy(nbf*nbf,stateout(wfn1_nxt_bocc()),1,dbl_mb(k_m2),1)
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
     &           dbl_mb(k_m2),nbf,0.0d0,stateout(wfn1_nxt_bocc()),nbf)
cDEBUG
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
      if (.not.ma_pop_stack(l_occb))
     &  call errquit("wfn1_nxt2_step: could not deallocate occb",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_occa))
     &  call errquit("wfn1_nxt2_step: could not deallocate occa",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
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
