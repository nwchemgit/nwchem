      Subroutine tddft_residual(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,
     1  ntrials,tda,g_trials,g_prod_p,g_prod_m,g_buffer1,g_buffer2,
     2  d_trials,d_prod_p,d_prod_m,apbvec,apbval,
     3  ambvec,residual_1,residual_2,g_residuals_1,g_residuals_2,
     4  thresh,valsav,maxresidual,maxomega,converged,maxvecs,evl,
     5  nconverged,nnewtrials,nresiduals,algorithm,lecut,ecut,
     6  diff_max,lowin,owstart,owend,lewin,ewinl,ewinh,vshift)
c
c $Id$
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
c R2 = (A-B)(X-Y) - w(X+Y)
c
c Some parts are explicitly parallelized, taking into account 
c the specific distribution (see the ga_create calls for g_residuals_1
c and g_residuals_2 in tddft_davidson).
c
c Written by So Hirata, Mar 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "errquit.fh"
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
#include "dra.fh"
c
      integer nroots           ! Number of roots sought
      integer nresiduals       ! Number of residual vectors
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! = 1 (RDFT); =2 (UDFT)
      integer nocc(2)          ! Number of occupied orbitals
      integer nmo(2)           ! Number of orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_p(2)      ! GA handle for (A+B) product vectors
      integer g_prod_m(2)      ! GA handle for (A-B) product vectors
      integer g_buffer1(2)     ! GA handle for buffer vectors
      integer g_buffer2(2)     ! GA handle for buffer vectors
      integer d_trials(2)      ! DRA handle for trial vectors
      integer d_prod_p(2)      ! DRA handle for (A+B) product vectors
      integer d_prod_m(2)      ! DRA handle for (A-B) product vectors
      integer g_residuals_1(2) ! GA handle for the residual vectors 1
      integer g_residuals_2(2) ! GA handle for the residual vectors 2
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      double precision apbvec(ntrials,ntrials)   ! (X+Y)
      double precision ambvec(ntrials,ntrials)   ! (X-Y)
      double precision apbval(ntrials)  ! (A-B)^(1/2)(A+B)(A-B)^(1/2) eigenvalues
      double precision residual_1(nresiduals) ! Norm of the residual vectors 1
      double precision residual_2(nresiduals) ! Norm of the residual vectors 2
      double precision thresh  ! Convergence threshold for Davidson iteration
      double precision valsav(nroots)   ! Excitation energies from previous iteration
      double precision maxresidual      ! Max Norm of the residual vectors
      double precision maxomega         ! Max change in excitation energies
      logical converged        ! True if convergence criterion is met
      integer maxvecs          ! Max number of trial vectors
      double precision evl(nbf_ao,2)    ! Orbital energies
      integer nconverged       ! Number of converged roots
      integer nnewtrials       ! Number of trial vectors in the next iteration
      integer algorithm        ! Algorithm
c
      integer i,j,k,l,n,m,k1,l1
      integer nov(2)
      double precision r
      logical nodezero
      integer iproc
      integer ilo,ihi,jlo,jhi
      integer dummy,request
      double precision ecut,evl_diff
      integer vshift
      logical lecut
      logical lowin
      logical lewin
      double precision ewinl,ewinh
      integer owstart(2)
      integer owend(2)
      double precision diff_max
      double precision tga_ddot
c
      iproc=ga_nodeid()
      nodezero=(ga_nodeid().eq.0)
c
c --------------------------------
c Number of occupied virtual pairs
c --------------------------------
c
      call tddft_lentrialvec(ipol,nocc,nmo,nfc,nfv,nov)
c
c ------------------
c Tamm-Dancoff AX=wX
c ------------------
c
      if (tda) then
        call tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,ntrials,
     1                  tda,0.0d0,g_trials,g_trials,g_buffer1,d_trials,
     2                  d_trials,apbvec,ambvec,g_residuals_1,
     3                  g_residuals_2,maxvecs,nresiduals,algorithm)
c
c Multiply (-w) to X
c
        do i=1,ipol
          do n=1,min(nresiduals,ntrials)
            call ga_scale_patch
     1        (g_residuals_1(i),1,nov(i),n,n,-apbval(n))
          enddo
        enddo
c
c Add AX to form AX-wX
c
        call tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,ntrials,
     1                  tda,1.0d0,g_prod_p,g_prod_m,g_buffer1,d_prod_p,
     2                  d_prod_m,apbvec,ambvec,g_residuals_1,
     3                  g_residuals_2,maxvecs,nresiduals,algorithm)
c
c ------------------------------------------------------------
c Full linear response (A+B)(X+Y)=w(X-Y) and (A-B)(X-Y)=w(X+Y)
c ------------------------------------------------------------
c
      else
        call tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,ntrials,
     1                  tda,0.0d0,g_trials,g_trials,g_buffer1,d_trials,
     2                  d_trials,apbvec,ambvec,g_residuals_1,
     3                  g_residuals_2,maxvecs,nresiduals,algorithm)
c
c Multiply (-w) to X+Y and X-Y
c
        do i=1,ipol
          do n=1,min(nresiduals,ntrials)
            call ga_scale_patch(g_residuals_1(i),1,nov(i),n,n,
     1        -apbval(n))
            call ga_scale_patch(g_residuals_2(i),1,nov(i),n,n,
     1        -apbval(n))
          enddo
        enddo
c
c Add (A-B)(X-Y) and (A+B)(X+Y)
c to form (A-B)(X-Y)-w(X+Y) and (A+B)(X+Y)-w(X-Y)
c
        call tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,ntrials,
     1                  tda,1.0d0,g_prod_m,g_prod_p,g_buffer1,d_prod_m,
     2                  d_prod_p,ambvec,apbvec,g_residuals_1,
     3                  g_residuals_2,maxvecs,nresiduals,algorithm)
      endif
c
c ---------------------------------
c Compute norms of residual vectors
c ---------------------------------
c
      do n=1,min(nresiduals,ntrials)
        residual_1(n)=0.0d0
        do i=1,ipol
          residual_1(n)=residual_1(n)+tga_ddot
     &      (g_residuals_1(i),g_residuals_1(i),n,n)
        enddo
        residual_1(n)=dsqrt(residual_1(n))
        if (.not.tda) then
          residual_2(n)=0.0d0
          do i=1,ipol
            residual_2(n)=residual_2(n)+tga_ddot
     &        (g_residuals_2(i),g_residuals_2(i),n,n)
          enddo
          residual_2(n)=dsqrt(residual_2(n))
        endif
      enddo
      if (nodezero.and.util_print('convergence',print_debug)) then
        write(LuOut,9000)
        write(LuOut,9010) (residual_1(n),n=1,min(nroots,20))
        if (.not.tda)
     1    write(LuOut,9010) (residual_2(n),n=1,min(nroots,20))
      endif
      if (util_print('residuals',print_debug)) then
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
c ------------------------------------------------------------------
c Check whether we have enough space available for new trial vectors
c if not form nroots eigenvectors and their matrix vector products
c and start Davidson from them.
c ------------------------------------------------------------------
c
      if (tda) then
        if (ntrials+min(nresiduals,ntrials)-nconverged.gt.maxvecs) then
          if (ga_nodeid().eq.0) then
            write(luout,6010)maxvecs,ntrials,
     &                       min(nresiduals,ntrials)-nconverged
c6000       format("maximum number of vectors: ",i4,/,
c    &             "current number of vectors: ",i4,/,
c    &             "number of new vectors    : ",i4,/,
c    &             "total number of vectors  : ",i4,/,
c    &             "maximum exceeded restarting Davidson solver")
          endif
c
c         Due to space constraints the algorithm will require a certain
c         amount of rework. Ideally we would want to use the residuals
c         calculated already but we do not have enough space to store
c         them yet. So we are forced to use the residual storage to 
c         construct the eigenvectors X and AX that we will replace the
c         trial vectors with. Subsequently resetting apbvec, ambvec, and
c         apbval we can then recalculate the residuals and continue as
c         usual.
c
c         Compute matrix vector products AX of eigenvectors X first
c
          call tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,ntrials,
     1                    tda,0.0d0,g_prod_p,g_prod_m,g_buffer1,
     2                    d_prod_p,d_prod_m,apbvec,ambvec,g_residuals_1,
     3                    g_residuals_2,maxvecs,nresiduals,algorithm)
          if ((algorithm.eq.1).or.(algorithm.eq.2)) then
            do i=1,ipol
              call ga_copy_patch('N',g_residuals_1(i),1,nov(i),1,nroots,
     +                               g_prod_p(i),1,nov(i),1,nroots)
            enddo
          else
            do i=1,ipol
              dummy=dra_write_section(.false.,
     1              g_residuals_1(i),1,nov(i),1,nroots,
     2              d_prod_p(i),1,nov(i),1,nroots,request)
              dummy=dra_wait(request)
            enddo
          endif
c
c         Compute eigenvectors X now
c
          call tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,ntrials,
     1                  tda,0.0d0,g_trials,g_trials,g_buffer1,d_trials,
     2                  d_trials,apbvec,ambvec,g_residuals_1,
     3                  g_residuals_2,maxvecs,nresiduals,algorithm)
          if ((algorithm.eq.1).or.(algorithm.eq.2)) then
            do i=1,ipol
              call ga_copy_patch('N',g_residuals_1(i),1,nov(i),1,nroots,
     +                               g_trials(i),1,nov(i),1,nroots)
            enddo
          else
            do i=1,ipol
              dummy=dra_write_section(.false.,
     1              g_residuals_1(i),1,nov(i),1,nroots,
     2              d_trials(i),1,nov(i),1,nroots,request)
              dummy=dra_wait(request)
            enddo
          endif
c
c         Scale the eigenvectors to obtain -wX
c
          do i=1,ipol
            do n=1,nroots
              call ga_scale_patch(g_residuals_1(i),1,nov(i),n,n,
     +                            -apbval(n))
            enddo
          enddo
c
c         Add the product vectors AX to reconstruct the residual vectors
c
          do i=1,ipol
            if ((algorithm.eq.1).or.(algorithm.eq.2)) then
              call tga_add(1.0d0,g_prod_p(i),g_residuals_1(i),
     &               1,1,nroots)
            else
              do n=1,nroots
                dummy=dra_read_section(.false.,
     1                g_buffer1(i),1,nov(i),1,1,
     2                d_prod_p(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
                call tga_add(1.0d0,g_buffer1(i),g_residuals_1(i),
     &               1,n,1)
              enddo
            endif
          enddo
c
c         Reset apbvec to correspond to the nroots*nroots unit matrix
c         (note this thing is declared as an ntrials*ntrials matrix)
c
          l1=1
          k1=1
          do l=1,nroots
            do k=1,nroots
              apbvec(k1,l1)=0.0d0
              if (k.eq.l) apbvec(k1,l1)=1.0d0
              k1=k1+1
              if (k1.gt.ntrials) then
                k1=1
                l1=l1+1
              endif
            enddo
          enddo
          ntrials=nroots
        endif
c      
      else ! (.not.tda i.e. full tddft)
c
        if (ntrials+2*(min(nresiduals,ntrials)-nconverged).gt.maxvecs)
     &  then
          if (ga_nodeid().eq.0) then
            write(luout,6010)maxvecs,ntrials,
     &                       2*(min(nresiduals,ntrials)-nconverged)
c
c           6010 currently the same as 6000 but this is expected to
c           change
c
 6010       format("maxvecs: ",i4," ntrials: ",i4," nnew: ",i4,/,
     &             "maximum exceeded restarting Davidson solver")
          endif
c
c         Due to space constraints the algorithm will require a certain
c         amount of rework. Ideally we would want to use the residuals
c         calculated already but we do not have enough space to store
c         them yet. So we are forced to use the residual storage to 
c         construct the eigenvectors (X+Y) and (X-Y), and the matrix-
c         vector products (A+B)(X+Y) and (A-B)(X-Y) that we will
c         replace the trial vectors with. Subsequently we can then
c         recalculate the residuals and continue as usual.
c
c         Compute matrix vector products AX of eigenvectors X first
c         1. Transform the (X+Y) vectors
c
          call tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,ntrials,
     1                    tda,0.0d0,g_prod_p,g_prod_p,g_buffer1,
     2                    d_prod_p,d_prod_p,apbvec,ambvec,g_residuals_1,
     3                    g_residuals_2,maxvecs,nresiduals,algorithm)
          if ((algorithm.eq.1).or.(algorithm.eq.2)) then
            do i=1,ipol
              call ga_copy_patch('N',g_residuals_1(i),1,nov(i),1,nroots,
     +                               g_prod_p(i),1,nov(i),1,nroots)
              call ga_copy_patch('N',g_residuals_2(i),1,nov(i),1,nroots,
     +                               g_prod_p(i),1,nov(i),
     +                                           nroots+1,2*nroots)
            enddo
          else
            do i=1,ipol
              dummy=dra_write_section(.false.,
     1              g_residuals_1(i),1,nov(i),1,nroots,
     2              d_prod_p(i),1,nov(i),1,nroots,request)
              dummy=dra_wait(request)
              dummy=dra_write_section(.false.,
     1              g_residuals_2(i),1,nov(i),1,nroots,
     2              d_prod_p(i),1,nov(i),nroots+1,2*nroots,request)
              dummy=dra_wait(request)
            enddo
          endif
c
c         2. Transform the (X-Y) vectors
c
          call tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,ntrials,
     1                    tda,0.0d0,g_prod_m,g_prod_m,g_buffer1,
     2                    d_prod_m,d_prod_m,apbvec,ambvec,g_residuals_1,
     3                    g_residuals_2,maxvecs,nresiduals,algorithm)
          if ((algorithm.eq.1).or.(algorithm.eq.2)) then
            do i=1,ipol
              call ga_copy_patch('N',g_residuals_1(i),1,nov(i),1,nroots,
     +                               g_prod_m(i),1,nov(i),1,nroots)
              call ga_copy_patch('N',g_residuals_2(i),1,nov(i),1,nroots,
     +                               g_prod_m(i),1,nov(i),
     +                                           nroots+1,2*nroots)
            enddo
          else
            do i=1,ipol
              dummy=dra_write_section(.false.,
     1              g_residuals_1(i),1,nov(i),1,nroots,
     2              d_prod_m(i),1,nov(i),1,nroots,request)
              dummy=dra_wait(request)
              dummy=dra_write_section(.false.,
     1              g_residuals_2(i),1,nov(i),1,nroots,
     2              d_prod_m(i),1,nov(i),nroots+1,2*nroots,request)
              dummy=dra_wait(request)
            enddo
          endif
c
c         Compute eigenvectors X now
c
          call tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,ntrials,
     1                  tda,0.0d0,g_trials,g_trials,g_buffer1,d_trials,
     2                  d_trials,apbvec,ambvec,g_residuals_1,
     3                  g_residuals_2,maxvecs,nresiduals,algorithm)
          if ((algorithm.eq.1).or.(algorithm.eq.2)) then
            do i=1,ipol
              call ga_copy_patch('N',g_residuals_1(i),1,nov(i),1,nroots,
     +                               g_trials(i),1,nov(i),1,nroots)
              call ga_copy_patch('N',g_residuals_2(i),1,nov(i),1,nroots,
     +                               g_trials(i),1,nov(i),
     +                                           nroots+1,2*nroots)
            enddo
          else
            do i=1,ipol
              dummy=dra_write_section(.false.,
     1              g_residuals_1(i),1,nov(i),1,nroots,
     2              d_trials(i),1,nov(i),1,nroots,request)
              dummy=dra_wait(request)
              dummy=dra_write_section(.false.,
     1              g_residuals_2(i),1,nov(i),1,nroots,
     2              d_trials(i),1,nov(i),nroots+1,2*nroots,request)
              dummy=dra_wait(request)
            enddo
          endif
c
c         Scale the eigenvectors to obtain -w(X+Y) and -w(X-Y)
c
          do i=1,ipol
            do n=1,nroots
              call ga_scale_patch(g_residuals_1(i),1,nov(i),n,n,
     +                            -apbval(n))
              call ga_scale_patch(g_residuals_2(i),1,nov(i),n,n,
     +                            -apbval(n))
            enddo
          enddo
c
c         Add the product vectors to reconstruct the residual vectors
c         (A-B)(X-Y)-w(X+Y) and (A+B)(X+Y)-w(X-Y)
c
          do i=1,ipol
            if ((algorithm.eq.1).or.(algorithm.eq.2)) then
              call tga_add(1.0d0,g_prod_m(i),g_residuals_1(i),
     &             nroots+1,1,nroots)
              call tga_add(1.0d0,g_prod_p(i),g_residuals_2(i),
     &             1,1,nroots)
            else
              do n=1,nroots
                dummy=dra_read_section(.false.,
     1                g_buffer1(i),1,nov(i),1,1,
     2                d_prod_m(i),1,nov(i),nroots+n,nroots+n,request)
                dummy=dra_wait(request)
                call tga_add(1.0d0,g_buffer1(i),g_residuals_1(i),
     &               1,n,1)
              enddo
              do n=1,nroots
                dummy=dra_read_section(.false.,
     1                g_buffer1(i),1,nov(i),1,1,
     2                d_prod_p(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
                call tga_add(1.0d0,g_buffer1(i),g_residuals_2(i),
     &               1,n,1)
              enddo
            endif
          enddo
c
c         Reset apbvec and ambvec to correspond to the 
c         (2*nroots)*(2*nroots) unit matrix (note these matrices are
c         declared as ntrials*ntrials matrices)
c
c         The intended approach above will not work. Due to the fact
c         that the Davidson subspace is not orthonormal apbvec and
c         ambvec should not correspond to unit matrices but something
c         else. Rather than trying to work the details of this out
c         we reset the matrices to zero so that attempts to use this
c         data will result in obvious errors.
c
          l1=1
          k1=1
          do l=1,2*nroots
            do k=1,2*nroots
              apbvec(k1,l1)=0.0d0
              ambvec(k1,l1)=0.0d0
c             if (k.eq.l) apbvec(k1,l1)=1.0d0
c             if (k+nroots.eq.l) ambvec(k1,l1)=1.0d0
c             if (k.eq.l+nroots) ambvec(k1,l1)=1.0d0
              k1=k1+1
              if (k1.gt.ntrials) then
                k1=1
                l1=l1+1
              endif
            enddo
          enddo
          ntrials=2*nroots
c
c         Re-orthonormalize the whole subspace
c
          do l=1,ntrials
            if ((algorithm.eq.1).or.(algorithm.eq.2)) then
              r = 0.0d0
              do i=1,ipol
                r = r + tga_ddot(g_trials(i),g_trials(i),l,l)
              enddo ! i=1,ipol
              r = 1.0d0/dsqrt(r)
              do i=1,ipol
                call ga_scale_patch(g_trials(i),1,nov(i),l,l,r)
                call ga_scale_patch(g_prod_p(i),1,nov(i),l,l,r)
                call ga_scale_patch(g_prod_m(i),1,nov(i),l,l,r)
              enddo ! i=1,ipol
            else
              r = 0.0d0
              do i=1,ipol
                dummy = dra_read_section(.false.,
     &                  g_buffer2(i),1,nov(i),1,1,
     &                  d_trials(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
                r = r + tga_ddot(g_buffer2(i),g_buffer2(i),1,1)
              enddo ! i=1,ipol
              r = 1.0d0/dsqrt(r)
              do i=1,ipol
                call ga_scale_patch(g_buffer2(i),1,nov(i),1,1,r)
                dummy = dra_write_section(.false.,
     &                  g_buffer2(i),1,nov(i),1,1,
     &                  d_trials(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
c
                dummy = dra_read_section(.false.,
     &                  g_buffer1(i),1,nov(i),1,1,
     &                  d_prod_p(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
                call ga_scale_patch(g_buffer1(i),1,nov(i),1,1,r)
                dummy = dra_write_section(.false.,
     &                  g_buffer1(i),1,nov(i),1,1,
     &                  d_prod_p(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
c
                dummy = dra_read_section(.false.,
     &                  g_buffer1(i),1,nov(i),1,1,
     &                  d_prod_m(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
                call ga_scale_patch(g_buffer1(i),1,nov(i),1,1,r)
                dummy = dra_write_section(.false.,
     &                  g_buffer1(i),1,nov(i),1,1,
     &                  d_prod_m(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
              enddo ! i=1,ipol
            endif ! ((algorithm.eq.1).or.(algorithm.eq.2))
            do k=1,l-1
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                r = 0.0d0
                do i=1,ipol
                  r = r + tga_ddot(g_trials(i),g_trials(i),k,l)
                enddo ! i=1,ipol
                do i=1,ipol
                  call tga_add(-r,g_trials(i),g_trials(i),k,l,1)
                  call tga_add(-r,g_prod_p(i),g_prod_p(i),k,l,1)
                  call tga_add(-r,g_prod_m(i),g_prod_m(i),k,l,1)
                enddo ! i=1,ipol
              else
                r = 0.0d0
                do i=1,ipol
                  dummy = dra_read_section(.false.,
     &                    g_buffer2(i),1,nov(i),1,1,
     &                    d_trials(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
                  dummy = dra_read_section(.false.,
     &                    g_buffer1(i),1,nov(i),1,1,
     &                    d_trials(i),1,nov(i),k,k,request)
                  dummy = dra_wait(request)
                  r = r + tga_ddot(g_buffer1(i),g_buffer2(i),1,1)
                enddo ! i=1,ipol
                do i=1,ipol
                  call tga_add(-r,g_buffer1(i),g_buffer2(i),1,1,1)
                  dummy = dra_write_section(.false.,
     &                    g_buffer2(i),1,nov(i),1,1,
     &                    d_trials(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
c
                  dummy = dra_read_section(.false.,
     &                    g_buffer2(i),1,nov(i),1,1,
     &                    d_prod_p(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
                  dummy = dra_read_section(.false.,
     &                    g_buffer1(i),1,nov(i),1,1,
     &                    d_prod_p(i),1,nov(i),k,k,request)
                  dummy = dra_wait(request)
                  call tga_add(-r,g_buffer1(i),g_buffer2(i),1,1,1)
                  dummy = dra_write_section(.false.,
     &                    g_buffer2(i),1,nov(i),1,1,
     &                    d_prod_p(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
c
                  dummy = dra_read_section(.false.,
     &                    g_buffer2(i),1,nov(i),1,1,
     &                    d_prod_m(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
                  dummy = dra_read_section(.false.,
     &                    g_buffer1(i),1,nov(i),1,1,
     &                    d_prod_m(i),1,nov(i),k,k,request)
                  dummy = dra_wait(request)
                  call tga_add(-r,g_buffer1(i),g_buffer2(i),1,1,1)
                  dummy = dra_write_section(.false.,
     &                    g_buffer2(i),1,nov(i),1,1,
     &                    d_prod_m(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
                enddo ! i=1,ipol
              endif
c             write(*,'(" overlap (",i2,",",i2,") = ",f22.14)')l,k,r
            enddo ! k=1,l-1
            if ((algorithm.eq.1).or.(algorithm.eq.2)) then
              r = 0.0d0
              do i=1,ipol
                r = r + tga_ddot(g_trials(i),g_trials(i),l,l)
              enddo ! i=1,ipol
              r = 1.0d0/dsqrt(r)
              do i=1,ipol
                call ga_scale_patch(g_trials(i),1,nov(i),l,l,r)
                call ga_scale_patch(g_prod_p(i),1,nov(i),l,l,r)
                call ga_scale_patch(g_prod_m(i),1,nov(i),l,l,r)
              enddo ! i=1,ipol
            else ! ((algorithm.eq.1).or.(algorithm.eq.2))
              r = 0.0d0
              do i=1,ipol
                dummy = dra_read_section(.false.,
     &                  g_buffer2(i),1,nov(i),1,1,
     &                  d_trials(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
                r = r + tga_ddot(g_buffer2(i),g_buffer2(i),1,1)
              enddo ! i=1,ipol
              r = 1.0d0/dsqrt(r)
              do i=1,ipol
                call ga_scale_patch(g_buffer2(i),1,nov(i),1,1,r)
                dummy = dra_write_section(.false.,
     &                  g_buffer2(i),1,nov(i),1,1,
     &                  d_trials(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
c
                dummy = dra_read_section(.false.,
     &                  g_buffer1(i),1,nov(i),1,1,
     &                  d_prod_p(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
                call ga_scale_patch(g_buffer1(i),1,nov(i),1,1,r)
                dummy = dra_write_section(.false.,
     &                  g_buffer1(i),1,nov(i),1,1,
     &                  d_prod_p(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
c
                dummy = dra_read_section(.false.,
     &                  g_buffer1(i),1,nov(i),1,1,
     &                  d_prod_m(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
                call ga_scale_patch(g_buffer1(i),1,nov(i),1,1,r)
                dummy = dra_write_section(.false.,
     &                  g_buffer1(i),1,nov(i),1,1,
     &                  d_prod_m(i),1,nov(i),l,l,request)
                dummy = dra_wait(request)
              enddo ! i=1,ipol
            endif ! ((algorithm.eq.1).or.(algorithm.eq.2))
            if (r.gt.1.0d5) then
              call errquit(
     &        "tddft_residual: linear dependency in restart",0,UERR)
            endif
            if (r.gt.dsqrt(2.0d0)) then
c
c             Repeat Gramm-Schmidt to ensure accurate ortogonality
c
              do k=1,l-1
                if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                  r = 0.0d0
                  do i=1,ipol
                    r = r + tga_ddot(g_trials(i),g_trials(i),k,l)
                  enddo ! i=1,ipol
                  do i=1,ipol
                    call tga_add(-r,g_trials(i),g_trials(i),k,l,1)
                    call tga_add(-r,g_prod_p(i),g_prod_p(i),k,l,1)
                    call tga_add(-r,g_prod_m(i),g_prod_m(i),k,l,1)
                  enddo ! i=1,ipol
                else ! ((algorithm.eq.1).or.(algorithm.eq.2))
                  r = 0.0d0
                  do i=1,ipol
                    dummy = dra_read_section(.false.,
     &                      g_buffer2(i),1,nov(i),1,1,
     &                      d_trials(i),1,nov(i),l,l,request)
                    dummy = dra_wait(request)
                    dummy = dra_read_section(.false.,
     &                      g_buffer1(i),1,nov(i),1,1,
     &                      d_trials(i),1,nov(i),k,k,request)
                    dummy = dra_wait(request)
                    r = r + tga_ddot(g_buffer1(i),g_buffer2(i),1,1)
                  enddo ! i=1,ipol
                  do i=1,ipol
                    call tga_add(-r,g_buffer1(i),g_buffer2(i),1,1,1)
                    dummy = dra_write_section(.false.,
     &                      g_buffer2(i),1,nov(i),1,1,
     &                      d_trials(i),1,nov(i),l,l,request)
                    dummy = dra_wait(request)
c
                    dummy = dra_read_section(.false.,
     &                      g_buffer2(i),1,nov(i),1,1,
     &                      d_prod_p(i),1,nov(i),l,l,request)
                    dummy = dra_wait(request)
                    dummy = dra_read_section(.false.,
     &                      g_buffer1(i),1,nov(i),1,1,
     &                      d_prod_p(i),1,nov(i),k,k,request)
                    dummy = dra_wait(request)
                    call tga_add(-r,g_buffer1(i),g_buffer2(i),1,1,1)
                    dummy = dra_write_section(.false.,
     &                      g_buffer2(i),1,nov(i),1,1,
     &                      d_prod_p(i),1,nov(i),l,l,request)
                    dummy = dra_wait(request)
c
                    dummy = dra_read_section(.false.,
     &                      g_buffer2(i),1,nov(i),1,1,
     &                      d_prod_m(i),1,nov(i),l,l,request)
                    dummy = dra_wait(request)
                    dummy = dra_read_section(.false.,
     &                      g_buffer1(i),1,nov(i),1,1,
     &                      d_prod_m(i),1,nov(i),k,k,request)
                    dummy = dra_wait(request)
                    call tga_add(-r,g_buffer1(i),g_buffer2(i),1,1,1)
                    dummy = dra_write_section(.false.,
     &                      g_buffer2(i),1,nov(i),1,1,
     &                      d_prod_m(i),1,nov(i),l,l,request)
                    dummy = dra_wait(request)
                  enddo ! i=1,ipol
                endif ! ((algorithm.eq.1).or.(algorithm.eq.2))
c               write(*,'(" overlap (",i2,",",i2,") = ",f22.14)')l,k,r
              enddo ! k=1,l-1
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                r = 0.0d0
                do i=1,ipol
                  r = r + tga_ddot(g_trials(i),g_trials(i),l,l)
                enddo ! i=1,ipol
                r = 1.0d0/dsqrt(r)
                do i=1,ipol
                  call ga_scale_patch(g_trials(i),1,nov(i),l,l,r)
                  call ga_scale_patch(g_prod_p(i),1,nov(i),l,l,r)
                  call ga_scale_patch(g_prod_m(i),1,nov(i),l,l,r)
                enddo ! i=1,ipol
              else ! ((algorithm.eq.1).or.(algorithm.eq.2))
                r = 0.0d0
                do i=1,ipol
                  dummy = dra_read_section(.false.,
     &                    g_buffer2(i),1,nov(i),1,1,
     &                    d_trials(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
                  r = r + tga_ddot(g_buffer2(i),g_buffer2(i),1,1)
                enddo ! i=1,ipol
                r = 1.0d0/dsqrt(r)
                do i=1,ipol
                  call ga_scale_patch(g_buffer2(i),1,nov(i),1,1,r)
                  dummy = dra_write_section(.false.,
     &                    g_buffer2(i),1,nov(i),1,1,
     &                    d_trials(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
c
                  dummy = dra_read_section(.false.,
     &                    g_buffer1(i),1,nov(i),1,1,
     &                    d_prod_p(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
                  call ga_scale_patch(g_buffer1(i),1,nov(i),1,1,r)
                  dummy = dra_write_section(.false.,
     &                    g_buffer1(i),1,nov(i),1,1,
     &                    d_prod_p(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
c
                  dummy = dra_read_section(.false.,
     &                    g_buffer1(i),1,nov(i),1,1,
     &                    d_prod_m(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
                  call ga_scale_patch(g_buffer1(i),1,nov(i),1,1,r)
                  dummy = dra_write_section(.false.,
     &                    g_buffer1(i),1,nov(i),1,1,
     &                    d_prod_m(i),1,nov(i),l,l,request)
                  dummy = dra_wait(request)
                enddo ! i=1,ipol
              endif ! ((algorithm.eq.1).or.(algorithm.eq.2))
            endif ! (r.gt.dsqrt(2.0d0))
          enddo
c      
        endif
      endif ! tda
c
c ------------------------
c Create new trial vectors
c ------------------------
c
      nnewtrials=ntrials
      do n=1,min(nresiduals,ntrials)
        if ((n.gt.nroots).and.(apbval(n)-apbval(nroots).gt.1.0d-4))
     1    continue
        if (residual_1(n).gt.thresh) then
          if (nnewtrials.gt.maxvecs) call errquit
     1      ('tddft_residual: out of memory',nnewtrials, MEM_ERR)
          call ga_sync()
c
          do i=1,ipol
           call ga_distribution
     1         (g_residuals_1(i),iproc,ilo,ihi,jlo,jhi)
           if ((jlo.le.n).and.(jhi.ge.n)) then
            do l=ilo,ihi
              if (abs(vshift).gt.0) then
                 k=nocc(i)+1+vshift
              else
                 k=mod(l-1,nmo(i)-nfv(i)-nocc(i))+nocc(i)+1
              end if
              j=(l-k+nocc(i))/(nmo(i)-nfv(i)-nocc(i))+nfc(i)+1
              call ga_get(g_residuals_1(i),l,l,n,n,r,1)
              evl_diff = evl(nmo(i),i)-evl(1,i) ! span the eigenvalues
              if (lecut) then
                if (evl(j,i).le.ecut) then
                  evl_diff = evl(k,i)-evl(j,i)
                end if
              else if (lowin) then
                if ((j.ge.owstart(i)).and.(j.le.owend(i))) then
                  evl_diff = evl(k,i)-evl(j,i)
                endif
              else if (lewin) then
                if (evl(j,i).le.ewinh .and. evl(j,i).ge.ewinl) then
                  evl_diff = evl(k,i)-evl(j,i)
                end if
              else
                evl_diff = evl(k,i)-evl(j,i)
              endif
              r=r/(evl_diff-apbval(n))
              call ga_put(g_residuals_1(i),l,l,n,n,r,1)
            enddo
           endif ! jlo, jhi check
          enddo
          call ga_sync()
c
c Schmidt orthogonalize against other trial vectors
c
          do m=1,nnewtrials
            r=0.0d0
            do i=1,ipol
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                r=r+tga_ddot(g_residuals_1(i),g_trials(i),n,m)
              else
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_trials(i),1,nov(i),m,m,request)
                dummy=dra_wait(request)
                r=r+tga_ddot(g_residuals_1(i),g_buffer1(i),n,1)
              endif
            enddo
            do i=1,ipol
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                call tga_add(-r,g_trials(i),g_residuals_1(i),m,n,1)
              else
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_trials(i),1,nov(i),m,m,request)
                dummy=dra_wait(request)
                call tga_add(-r,g_buffer1(i),g_residuals_1(i),1,n,1)
              endif
            enddo
          enddo
c
c Normalize the new trial vector
c
          r=0.0d0
          do i=1,ipol
            r=r+tga_ddot(g_residuals_1(i),g_residuals_1(i),n,n)
          enddo
          if (r.gt.1.0d-15) then
            r=dsqrt(r)
            nnewtrials=nnewtrials+1
            do i=1,ipol
              call ga_scale_patch(g_residuals_1(i),1,nov(i),n,n,1.0d0/r)
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                call ga_copy_patch('N',g_residuals_1(i),1,nov(i),n,n,
     1            g_trials(i),1,nov(i),nnewtrials,nnewtrials)
              else
                dummy=dra_write_section(.false.,
     1            g_residuals_1(i),1,nov(i),n,n,
     2            d_trials(i),1,nov(i),nnewtrials,nnewtrials,request)
                dummy=dra_wait(request)
              endif
            enddo
          endif
        endif
      enddo
c
c Full linear response may generate another set of new trial vectors
c
      if (.not.tda) then
        do n=1,min(nresiduals,ntrials)
        if ((n.gt.nroots).and.(apbval(n)-apbval(nroots).gt.1.0d-4))
     1    continue
          if (residual_2(n).gt.thresh) then
            if (nnewtrials.gt.maxvecs) call errquit
     1        ('tddft_residual: out of memory',nnewtrials, MEM_ERR)
c
c Scale by 1/{w(epsilon_a-epsilon_i)^(1/2)}
c
            call ga_sync()
            do i=1,ipol
             call ga_distribution
     1           (g_residuals_2(i),iproc,ilo,ihi,jlo,jhi)
             if ((jlo.le.n).and.(jhi.ge.n)) then
              do l=ilo,ihi
                if (abs(vshift).gt.0) then
                   k=nocc(i)+1+vshift
                else
                   k=mod(l-1,nmo(i)-nfv(i)-nocc(i))+nocc(i)+1
                end if
                j=(l-k+nocc(i))/(nmo(i)-nfv(i)-nocc(i))+nfc(i)+1
                call ga_get(g_residuals_2(i),l,l,n,n,r,1)
                evl_diff = evl(nmo(i),i)-evl(1,i) ! span the eigenvalues
                if (lecut) then
                  if (evl(j,i).le.ecut) then
                    evl_diff = evl(k,i)-evl(j,i)
                  end if
                else if (lowin) then
                  if((j.ge.owstart(i)).and.(j.le.owend(i))) then
                    evl_diff = evl(k,i)-evl(j,i)
                  endif
                else if (lewin) then
                  if (evl(j,i).le.ewinh .and. evl(j,i).ge.ewinl) then
                    evl_diff = evl(k,i)-evl(j,i)
                  end if
                else
                  evl_diff = evl(k,i)-evl(j,i)
                endif
                r=r/(evl_diff-apbval(n))
                call ga_put(g_residuals_2(i),l,l,n,n,r,1)
              enddo
             endif ! jlo, jhi check
            enddo
            call ga_sync()
c
c Schmidt orthogonalize against other trial vectors
c
            do m=1,nnewtrials
              r=0.0d0
              do i=1,ipol
                if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                  r=r+tga_ddot(g_residuals_2(i),g_trials(i),n,m)
                else
                  dummy=dra_read_section(.false.,
     1              g_buffer1(i),1,nov(i),1,1,
     2              d_trials(i),1,nov(i),m,m,request)
                  dummy=dra_wait(request)
                  r=r+tga_ddot(g_residuals_2(i),g_buffer1(i),n,1)
                endif
              enddo
              do i=1,ipol
                if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                  call tga_add(-r,g_trials(i),g_residuals_2(i),m,n,1)
                else
                  dummy=dra_read_section(.false.,
     1              g_buffer1(i),1,nov(i),1,1,
     2              d_trials(i),1,nov(i),m,m,request)
                  dummy=dra_wait(request)
                  call tga_add(-r,g_buffer1(i),g_residuals_2(i),1,n,1)
                endif
              enddo
            enddo
c
c Normalize the new trial vector
c
            r=0.0d0
            do i=1,ipol
              r=r+tga_ddot(g_residuals_2(i),g_residuals_2(i),n,n)
            enddo
            if (r.gt.1.0d-15) then
              r=dsqrt(r)
              nnewtrials=nnewtrials+1
              do i=1,ipol
                call ga_scale_patch
     1            (g_residuals_2(i),1,nov(i),n,n,1.0d0/r)
                if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                  call ga_copy_patch('N',g_residuals_2(i),1,nov(i),n,n,
     1              g_trials(i),1,nov(i),nnewtrials,nnewtrials)
                else
                  dummy=dra_write_section(.false.,
     1              g_residuals_2(i),1,nov(i),n,n,
     2              d_trials(i),1,nov(i),nnewtrials,nnewtrials,request)
                  dummy=dra_wait(request)
                endif
              enddo
            endif
          endif
        enddo
      endif
c 
c Have we expand the subspace size?
c
      if ((.not.converged).and.(nnewtrials.eq.ntrials)) call errquit
     1  ('tddft_residual: failed to expand subspace size',0, MEM_ERR)
c
c ------
c Return
c ------
c
      return
      end
