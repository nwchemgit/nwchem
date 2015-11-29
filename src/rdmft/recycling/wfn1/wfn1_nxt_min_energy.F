C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Optimize the state vector to minimize the total energy
C>
C> This routine minimizes the total energy of the state vector 
C> while satisfying the various constraints. The general approach is
C> to do linesearches along the gradient to find the minimum energy
C> in that particular direction. Then the step to that point is taken,
C> the energy and the gradient calculated and a new linesearch started.
C> This procedure continues until the gradient is small.
C>
      subroutine wfn1_nxt_min_energy(rtdb,geom,basis,nbf,nea,neb,
     &           state,etot)
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
      double precision state(8*nbf*nbf) !< [In/Output] The state vector
      double precision etot             !< [Output] The total energy
c
      double precision dd       !< Norm of the gradient
      double precision ds       !< Step length from gradient
      double precision step     !< Step size from linesearch
      double precision tol_conv !< Convergence tolerance
c
      integer l_h1,   k_h1   !< Memory handles for 1-electron integrals
      integer l_ov,   k_ov   !< Memory handles for overlap integrals
      integer l_eri,  k_eri  !< Memory handles for Coulomb integrals
      integer l_erix, k_erix !< Memory handles for exchange integrals
      integer l_grad, k_grad !< Memory handles for the gradient
c
      integer it     !< Iteration count
      integer maxit  !< Maximum number of iterations
      logical status !< Status variable
c
      double precision wfn1_nxt_dot
      external         wfn1_nxt_dot
c
      maxit = 5000
      tol_conv = 1.0d-5
      status = rtdb_get(rtdb,"wfn1:maxit",MT_INT,1,maxit)
      status = rtdb_get(rtdb,"wfn1:thresh",MT_DBL,1,tol_conv)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"h1",l_h1,k_h1))
     &   call errquit("wfn1_nxt_min_energy: could not allocate h1",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"ov",l_ov,k_ov))
     &   call errquit("wfn1_nxt_min_energy: could not allocate ov",
     &   ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nbf*nbf,"eri",l_eri,k_eri))
     &   call errquit("wfn1_nxt_min_energy: could not allocate eri",
     &   ma_sizeof(MT_DBL,nbf*nbf*nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nbf*nbf,"erix",l_erix,k_erix))
     &   call errquit("wfn1_nxt_min_energy: could not allocate erix",
     &   ma_sizeof(MT_DBL,nbf*nbf*nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"grad",
     &   l_grad,k_grad))
     &   call errquit("wfn1_nxt_min_energy: could not allocate grad",
     &   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
c
c     Now initialize the memory segments and then compute the integrals
c
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_h1),1)
      call dfill(nbf*nbf,0.0d0,dbl_mb(k_ov),1)
      call dfill(nbf*nbf*nbf*nbf,0.0d0,dbl_mb(k_eri),1)
      call dfill(nbf*nbf*nbf*nbf,0.0d0,dbl_mb(k_erix),1)
c
      call int_init(rtdb,1,basis)
      call wfn1_1e_tv(basis,nbf,dbl_mb(k_h1))
      call wfn1_1e_s(basis,nbf,dbl_mb(k_ov))
      call wfn1_2e_rep(basis,nbf,dbl_mb(k_eri),dbl_mb(k_erix))
      call int_terminate
c
c     Calculate the initial energy and gradient
c
      call wfn1_nxt_energy_gradient(rtdb,geom,nbf,nea,neb,dbl_mb(k_h1),
     &     dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),state,
     &     dbl_mb(k_grad),etot)
cDEBUG
      write(*,*)'*** Initial energy = ',etot
      write(*,*)'*** Initial gradient'
      call wfn1_nxt_print_state(dbl_mb(k_grad),nbf)
cDEBUG
      call wfn1_nxt_gradient_to_search(nbf,dbl_mb(k_grad))
c
c     At this point the Lagrangian multipliers are zero, also the
c     gradient wrt. the Lagrangian multipliers is zero. So we can only
c     take a small step down hill. Doing this twice should move the
c     Lagrangian multipliers off zero. Then we can start the
c     optimization proper.
c
c     dd = wfn1_nxt_dot(wfn1_nxt_size1(),dbl_mb(k_grad),1,
c    &                  dbl_mb(k_grad),1)
      dd = wfn1_nxt_dot(nbf,dbl_mb(k_ov),dbl_mb(k_grad),
     &                  dbl_mb(k_grad))
      dd = sqrt(dd)
cDEBUG
      write(*,*)'*** Norm gradient = ',dd
cDEBUG
      ds = 0.01d0/dd
      call daxpy(wfn1_nxt_size1(),-ds,dbl_mb(k_grad),1,state,1)
c
c     Do step 2
c
      call wfn1_nxt_energy_gradient(rtdb,geom,nbf,nea,neb,dbl_mb(k_h1),
     &     dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),state,
     &     dbl_mb(k_grad),etot)
cDEBUG
      write(*,*)'*** 2nd energy = ',etot
      write(*,*)'*** 2nd gradient'
      call wfn1_nxt_print_state(dbl_mb(k_grad),nbf)
cDEBUG
      call wfn1_nxt_gradient_to_search(nbf,dbl_mb(k_grad))
c     dd = wfn1_nxt_dot(wfn1_nxt_size1(),dbl_mb(k_grad),1,
c    &                  dbl_mb(k_grad),1)
      dd = wfn1_nxt_dot(nbf,dbl_mb(k_ov),dbl_mb(k_grad),
     &                  dbl_mb(k_grad))
      dd = sqrt(dd)
cDEBUG
      write(*,*)'*** 2nd Norm gradient = ',dd
cDEBUG
      ds = 0.01d0/dd
      call daxpy(wfn1_nxt_size1(),-ds,dbl_mb(k_grad),1,state,1)
c
c     Now optimize properly
c
      it = 0
      do while (dd.gt.tol_conv.and.it.le.maxit)
        it = it + 1
        call wfn1_nxt_energy_gradient(rtdb,geom,nbf,nea,neb,
     &       dbl_mb(k_h1),dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),
     &       state,dbl_mb(k_grad),etot)
        call wfn1_nxt_gradient_to_search(nbf,dbl_mb(k_grad))
c       dd = wfn1_nxt_dot(wfn1_nxt_size1(),dbl_mb(k_grad),1,
c    &                    dbl_mb(k_grad),1)
        dd = wfn1_nxt_dot(nbf,dbl_mb(k_ov),
     &       dbl_mb(k_grad),dbl_mb(k_grad))
        dd = sqrt(dd)
cDEBUG
c       write(*,*)'wfn1 it: ',it,etot,dd
cDEBUG
        call wfn1_nxt_linesearch(rtdb,geom,nbf,nea,neb,step,
     &       dbl_mb(k_h1),dbl_mb(k_eri),dbl_mb(k_erix),dbl_mb(k_ov),
     &       state,dbl_mb(k_grad),etot,tol_conv)
        if (step.lt.tol_conv.and.dd.gt.tol_conv) then
          step = 1.0d0/100d0
        endif
        call daxpy(wfn1_nxt_size1(),-step,dbl_mb(k_grad),1,state,1)
cDEBUG
        if (mod(it,100).eq.0) then
          write(*,*)'wfn1 it: ',it,etot,dd,step
        endif
cDEBUG

      enddo
cDEBUG
      write(*,*)'*** last energy = ',etot
      write(*,*)'*** last gradient'
      call wfn1_nxt_print_state(dbl_mb(k_grad),nbf)
      write(*,*)'*** last state'
      call wfn1_nxt_print_state(state,nbf)
cDEBUG
c
      end
C>
C> @}
