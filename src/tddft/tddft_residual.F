      Subroutine tddft_residual(nroots,nbf_ao,ipol,nocc,nmo,ntrials,
     1  tda,g_trials,g_prod_p,g_prod_m,g_apbvec,apbval,g_ambvec,
     2  residual_1,residual_2,g_residuals_1,g_residuals_2,thresh,
     3  valsav,maxresidual,maxomega,converged,maxvecs,evl,
     4  nconverged,nnewtrials)
c
c $Id: tddft_residual.F,v 1.1 2002-03-22 02:10:37 sohirata Exp $
c
c Calculate the residual vector and create a new trial
c vector, if convergence is not achieved.
c
c Tamm-Dancoff approximation:
c
c R = AX - wX
c
c Full linear response theory:
c
c (X+Y) and (X-Y) are bi-orthonormalized
c
c R1 = (A+B)(X+Y) - w(X-Y)
c R2 = (A+B)(X-Y) - w(X+Y)
c
c Written by So Hirata, Mar 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
c
      integer nroots           ! Number of roots sought
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! = 1 (RDFT); =2 (UDFT)
      integer nocc(2)          ! Number of occupied orbitals
      integer nmo(2)           ! Number of orbitals
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_p(2)      ! GA handle for (A+B) product vectors
      integer g_prod_m(2)      ! GA handle for (A-B) product vectors
      integer g_residuals_1(2) ! GA handle for the residual vectors 1
      integer g_residuals_2(2) ! GA handle for the residual vectors 2
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      integer g_apbvec         ! GA handle for (X+Y)
      integer g_ambvec         ! GA handle for (X-Y)
      double precision apbval(ntrials)  ! (A-B)^(1/2)(A+B)(A-B)^(1/2) eigenvalues
      double precision residual_1(nroots) ! Norm of the residual vectors 1
      double precision residual_2(nroots) ! Norm of the residual vectors 2
      double precision thresh  ! Convergence threshold for Davidson iteration
      double precision valsav(nroots)   ! Excitation energies from previous iteration
      double precision maxresidual      ! Max Norm of the residual vectors
      double precision maxomega         ! Max change in excitation energies
      logical converged        ! True if convergence criterion is met
      integer maxvecs          ! Max number of trial vectors
      double precision evl(nbf_ao,2)    ! Orbital energies
      integer nconverged       ! Number of converged roots
      integer nnewtrials       ! Number of trial vectors in the next iteration
c
      integer i,j,k,l,n,m
      integer nov(2)
      double precision r
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c --------------------------------
c Number of occupied virtual pairs
c --------------------------------
c
      do i=1,ipol
        nov(i)=(nmo(i)-nocc(i))*nocc(i)
      enddo
      if (ipol.eq.1) nov(2)=0
c
c ------------------
c Tamm-Dancoff AX=wX
c ------------------
c
      if (tda) then
c
c Form X
c
        do i=1,ipol
          call ga_dgemm('N','N',nov(i),nroots,ntrials,
     1      1.0d0,g_trials(i),g_apbvec,0.0d0,g_residuals_1(i))
c
c Multiply (-w) to X
c
          do n=1,nroots
            call ga_scale_patch
     1        (g_residuals_1(i),1,nov(i),n,n,-apbval(n))
          enddo
c
c Add AX to form AX-wX
c
          call ga_dgemm('N','N',nov(i),nroots,ntrials,
     1      1.0d0,g_prod_p(i),g_apbvec,1.0d0,g_residuals_1(i))
        enddo  
c
c ------------------------------------------------------------
c Full linear response (A+B)(X+Y)=w(X-Y) and (A-B)(X-Y)=w(X+Y)
c ------------------------------------------------------------
c
      else
c
c Form X+Y and X-Y
c
        do i=1,ipol
          call ga_dgemm('N','N',nov(i),nroots,ntrials,
     1      1.0d0,g_trials(i),g_apbvec,0.0d0,g_residuals_1(i))
          call ga_dgemm('N','N',nov(i),nroots,ntrials,
     1      1.0d0,g_trials(i),g_ambvec,0.0d0,g_residuals_2(i))
c
c Multiply (-w) to X+Y and X-Y
c
          do n=1,nroots
            call ga_scale_patch(g_residuals_1(i),1,nov(i),n,n,
     1        -apbval(n))
            call ga_scale_patch(g_residuals_2(i),1,nov(i),n,n,
     1      -apbval(n))
          enddo
c
c Add (A-B)(X-Y) and (A+B)(X+Y)
c to form (A-B)(X-Y)-w(X+Y) and (A+B)(X+Y)-w(X-Y)
c
          call ga_dgemm('N','N',nov(i),nroots,ntrials,
     1      1.0d0,g_prod_m(i),g_ambvec,1.0d0,g_residuals_1(i))
          call ga_dgemm('N','N',nov(i),nroots,ntrials,
     1      1.0d0,g_prod_p(i),g_apbvec,1.0d0,g_residuals_2(i))
        enddo  
      endif
c
c ---------------------------------
c Compute norms of residual vectors
c ---------------------------------
c
      do n=1,nroots
        residual_1(n)=0.0d0
        do i=1,ipol
          residual_1(n)=residual_1(n)+ga_ddot_patch
     1      (g_residuals_1(i),'N',1,nov(i),n,n,
     2       g_residuals_1(i),'N',1,nov(i),n,n)
        enddo
        residual_1(n)=dsqrt(residual_1(n))
        if (.not.tda) then
          residual_2(n)=0.0d0
          do i=1,ipol
            residual_2(n)=residual_2(n)+ga_ddot_patch
     1        (g_residuals_2(i),'N',1,nov(i),n,n,
     2         g_residuals_2(i),'N',1,nov(i),n,n)
          enddo
          residual_2(n)=dsqrt(residual_2(n))
        endif
      enddo
      if (nodezero.and.util_print('convergence',print_high)) then
        write(LuOut,9000)
        write(LuOut,9010) (residual_1(n),n=1,min(nroots,20))
        if (.not.tda)
     1    write(LuOut,9010) (residual_2(n),n=1,min(nroots,20))
      endif
      if (util_print('convergence',print_debug)) then
        do i=1,ipol
          call ga_print(g_residuals_1(i))
        enddo
        if (.not.tda) then
          do i=1,ipol
            call ga_print(g_residuals_2(i))
          enddo
        endif
      endif
 9000 format(2x,'Norm of residuals')
 9010 format(2x,20f12.8)
c
c -----------------
c Check convergence
c -----------------
c
      maxresidual=0.0d0
      nconverged=0
      do n=1,nroots
        if (tda) then
          if (residual_1(n).gt.maxresidual) maxresidual=residual_1(n)
          if (residual_1(n).lt.thresh) nconverged=nconverged+1
        else
          if (residual_1(n).gt.maxresidual) maxresidual=residual_1(n)
          if (residual_2(n).gt.maxresidual) maxresidual=residual_2(n)
          if (max(residual_1(n),residual_2(n)).lt.thresh) 
     1      nconverged=nconverged+1
        endif
      enddo
      maxomega=0.0d0
      do n=1,nroots
        if (dabs(apbval(n)-valsav(n)).gt.maxomega)
     1    maxomega=dabs(apbval(n)-valsav(n))
        valsav(n)=apbval(n)
      enddo
      if (maxresidual.lt.thresh) then
        converged=.true.
      else
        converged=.false.
      endif
c
c ------------------------
c Create new trial vectors
c ------------------------
c CAUTION: this part is performed only by node zero; perhaps I should
c parallelize this part as well in the future.
c
      nnewtrials=ntrials
      do n=1,nroots
        if (residual_1(n).gt.thresh) then
          if (ntrials.gt.maxvecs) call errquit
     1      ('Out of memory',ntrials)
          if (nodezero) then
            do i=1,ipol
              l=0
              do j=1,nocc(i)
                do k=nocc(i)+1,nmo(i)
                  l=l+1
                  call ga_get(g_residuals_1(i),l,l,n,n,r,1)
                  r=r/(evl(k,i)-evl(j,i)-apbval(n))
                  call ga_put(g_residuals_1(i),l,l,n,n,r,1)
                enddo
              enddo
            enddo
          endif
c
c Schmitt orthogonalize against other trial vectors
c
          do m=1,nnewtrials
            r=0.0d0
            do i=1,ipol
              r=r+ga_ddot_patch(g_residuals_1(i),'N',1,nov(i),n,n,
     1                          g_trials(i),     'N',1,nov(i),m,m)
            enddo
            do i=1,ipol
              call ga_add_patch(1.0d0,g_residuals_1(i),1,nov(i),n,n,
     1                          -r,   g_trials(i),     1,nov(i),m,m,
     2                                g_residuals_1(i),1,nov(i),n,n)
            enddo
          enddo
c
c Normalize the new trial vector
c
          r=0.0d0
          do i=1,ipol
            r=r+ga_ddot_patch(g_residuals_1(i),'N',1,nov(i),n,n,
     1                        g_residuals_1(i),'N',1,nov(i),n,n)
          enddo
          if (r.gt.1.0d-15) then
            r=dsqrt(r)
            nnewtrials=nnewtrials+1
            do i=1,ipol
              call ga_scale_patch(g_residuals_1(i),1,nov(i),n,n,1.0d0/r)
              call ga_copy_patch('N',g_residuals_1(i),1,nov(i),n,n,
     1          g_trials(i),1,nov(i),nnewtrials,nnewtrials)
            enddo
          endif
        endif
      enddo
c
c Full linear response may generate another set of new trial vectors
c
      if (.not.tda) then
        do n=1,nroots
          if (residual_2(n).gt.thresh) then
            if (ntrials.gt.maxvecs) call errquit
     1        ('Out of memory',ntrials)
c
c Scale by 1/{w(epsilon_a-epsilon_i)^(1/2)}
c
            if (nodezero) then
              do i=1,ipol
                l=0
                do j=1,nocc(i)
                  do k=nocc(i)+1,nmo(i)
                    l=l+1
                    call ga_get(g_residuals_2(i),l,l,n,n,r,1)
                    r=r/(evl(k,i)-evl(j,i)-apbval(n))
                    call ga_put(g_residuals_2(i),l,l,n,n,r,1)
                  enddo
                enddo
              enddo
            endif
c
c Schmitt orthogonalize against other trial vectors
c
            do m=1,nnewtrials
              r=0.0d0
              do i=1,ipol
                r=r+ga_ddot_patch(g_residuals_2(i),'N',1,nov(i),n,n,
     1                            g_trials(i),     'N',1,nov(i),m,m)
              enddo
              do i=1,ipol
                call ga_add_patch(1.0d0,g_residuals_2(i),1,nov(i),n,n,
     1                            -r,   g_trials(i),     1,nov(i),m,m,
     2                                  g_residuals_2(i),1,nov(i),n,n)
              enddo
            enddo
c
c Normalize the new trial vector
c
            r=0.0d0
            do i=1,ipol
              r=r+ga_ddot_patch(g_residuals_2(i),'N',1,nov(i),n,n,
     1                          g_residuals_2(i),'N',1,nov(i),n,n)
            enddo
            if (r.gt.1.0d-15) then
              r=dsqrt(r)
              nnewtrials=nnewtrials+1
              do i=1,ipol
                call ga_scale_patch
     1            (g_residuals_2(i),1,nov(i),n,n,1.0d0/r)
                call ga_copy_patch('N',g_residuals_2(i),1,nov(i),n,n,
     1            g_trials(i),1,nov(i),nnewtrials,nnewtrials)
              enddo
            endif
          endif
        enddo
      endif
c 
c Have we expand the subspace size?
c
      if ((.not.converged).and.(nnewtrials.eq.ntrials)) call errquit
     1  ('tddft_residual: failed to expand subspace size',0)
c
c ------
c Return
c ------
c
      return
      end
