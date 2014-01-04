      Subroutine tddft_fold(nroots,nbf_ao,ipol,nocc,nmo,nfc,nfv,
     1  ntrials,tda,alpha,
     &  g_prod_p,g_prod_m,g_buffer1,
     2  d_prod_p,d_prod_m,apbvec,
     3  ambvec,g_residuals_1,g_residuals_2,
     4  maxvecs,
     5  nresiduals,algorithm)
c
c $Id$
c
c Calculate linear combinations of (trial or product) vectors.
c The vectors are returned in g_residuals_1 and g_residuals_2.
c
c Tamm-Dancoff approximation:
c
c R(j) = alpha*R(j) + (i=1:ntrials,apbvec(i,j)*g_prod_p(i))
c
c Full linear response theory:
c
c (X+Y) and (X-Y) are bi-orthonormalized
c
c R1(j) = alpha*R1(j) + (i=1:ntrials,apbvec(i,j)*g_prod_p(i))
c R2(j) = alpha*R2(j) + (i=1:ntrials,ambvec(i,j)*g_prod_m(i))
c
c Some parts are explicitly parallelized, taking into account 
c the specific distribution (see the ga_create calls for g_residuals_1
c and g_residuals_2 in tddft_davidson).
c
c Created by Huub van Dam, Jan 2011. (c) Battelle, PNNL, 2011.
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
      integer g_prod_p(2)      ! GA handle for (A+B) product vectors
      integer g_prod_m(2)      ! GA handle for (A-B) product vectors
      integer g_buffer1(2)     ! GA handle for buffer vectors
      integer d_prod_p(2)      ! DRA handle for (A+B) product vectors
      integer d_prod_m(2)      ! DRA handle for (A-B) product vectors
      integer g_residuals_1(2) ! GA handle for the residual vectors 1
      integer g_residuals_2(2) ! GA handle for the residual vectors 2
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      double precision apbvec(ntrials,ntrials)   ! (X+Y)
      double precision ambvec(ntrials,ntrials)   ! (X-Y)
      double precision alpha  ! scale factor for data currently in g_residual
      integer maxvecs          ! Max number of trial vectors
      integer algorithm        ! Algorithm
c
      integer i,j,k,l,n,m
      integer nov(2)
      double precision r
      logical nodezero
      integer iproc
      integer ilo,ihi,jlo,jhi
      integer dummy,request
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
        do i=1,ipol
          call ga_scale(g_residuals_1(i),alpha)
c
c Form X
c
          do n=1,min(nresiduals,ntrials)
            do m=1,ntrials
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                call tga_add(apbvec(m,n),g_prod_p(i),g_residuals_1(i),
     &               m,n,1)
              else
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_prod_p(i),1,nov(i),m,m,request)
                dummy=dra_wait(request)
                call tga_add(apbvec(m,n),g_buffer1(i),g_residuals_1(i),
     &               1,n,1)
              endif
            enddo ! m=1,ntrials
          enddo ! n=1,min(nresiduals,ntrials)
        enddo ! i=1,ipol
c
c ------------------------------------------------------------
c Full linear response (A+B)(X+Y)=w(X-Y) and (A-B)(X-Y)=w(X+Y)
c ------------------------------------------------------------
c
      else
        do i=1,ipol
          call ga_scale(g_residuals_1(i),alpha)
          call ga_scale(g_residuals_2(i),alpha)
c
c Form (X+Y) and (X-Y)
c
          do n=1,min(nresiduals,ntrials)
            do m=1,ntrials
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                call tga_add(apbvec(m,n),g_prod_p(i),g_residuals_1(i),
     &               m,n,1)
                call tga_add(ambvec(m,n),g_prod_m(i),g_residuals_2(i),
     &               m,n,1)
              else
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_prod_p(i),1,nov(i),m,m,request)
                dummy=dra_wait(request)
                call tga_add(apbvec(m,n),g_buffer1(i),g_residuals_1(i),
     &               1,n,1)
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,
     2            d_prod_m(i),1,nov(i),m,m,request)
                dummy=dra_wait(request)
                call tga_add(ambvec(m,n),g_buffer1(i),g_residuals_2(i),
     &               1,n,1)
              endif
            enddo ! m=1,ntrials
          enddo ! n=1,min(nresiduals,ntrials)
        enddo ! i=1,ipol
      endif
c
c ------
c Return
c ------
c
      return
      end
