      subroutine tddft_subspace(ipol,nocc,nmo,nfc,nfv,g_trials,
     1  g_prod_p,g_prod_m,ntrials,tda,apb,amb)
c
c $Id: tddft_subspace.F,v 1.4 2002-04-16 01:46:28 sohirata Exp $
c
c Form subspace representation of (A+B) and (A-B) matrices
c using trial and product vectors.  The elements of (A+B)
c and (A-B) matrices are simply:
c
c (A+B)_pq = (A+B) product vector (p) * trial vector (q)
c (A-B)_pq = (A-B) product vector (p) * trial vector (q)
c
c This part is explicitly parallelized.  Note that the g_trials
c and g_prod_p and g_prod_m are chunked in the first argument,
c but each processor has contributions from all vectors. 
c This subroutine will break, if someone changes the definition
c of these vectors, especially how they are chunked.
c
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
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
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_p(2)      ! GA handle for A+B product vectors
      integer g_prod_m(2)      ! GA handle for A-B product vectors
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      double precision apb(ntrials,ntrials) ! Subspace (A+B) matrix
      double precision amb(ntrials,ntrials) ! Subspace (A-B) matrix
c
      integer nov(2)           ! Number of occupied-virtual pairs
      integer iproc            ! Process number
      integer ilo,ihi,jlo,jhi  ! Range of distribution
      integer maxnov           ! Leading dimension of work array
      integer l_trial,k_trial  ! MA for trial vector strips
      integer l_prod,k_prod    ! MA for product vector strips
      integer g_apb_tmp        ! Temporary GA array to accumulate apb
      integer g_amb_tmp        ! Temporary GA array to accumulate amb
      integer i,l,m,n
c
c --------------------
c Allocate work arrays
c --------------------
c
      do i=1,ipol
        nov(i)=(nmo(i)-nfv(i)-nocc(i))*(nocc(i)-nfc(i))
      enddo
      if (ipol.eq.1) nov(2)=0
      maxnov=max(nov(1),nov(2))
      if (.not.ma_push_get(mt_dbl,maxnov*ntrials,'trial',
     1  l_trial,k_trial)) call errquit
     2  ('tddft_subspace: failed to allocate trial',0)
      if (.not.ma_push_get(mt_dbl,maxnov*ntrials,'product',
     1  l_prod,k_prod)) call errquit
     2  ('tddft_subspace: failed to allocate prod',0)
c
c ---------------------------------------------
c Form a local contribution of the dot products 
c ---------------------------------------------
c
      do m=1,ntrials
        do n=1,ntrials
          apb(m,n)=0.0d0
          if (.not.tda) amb(m,n)=0.0d0
        enddo
      enddo
      iproc=ga_nodeid()
      call ga_sync()
      do i=1,ipol
        call ga_distribution(g_trials(i),iproc,ilo,ihi,jlo,jhi)
        if ((jlo.gt.1).or.(jhi.lt.ntrials))
     1    call errquit('tddft_subspace: an internal error detected',0)
        call ga_get(g_trials(i),ilo,ihi,1,ntrials,
     1    dbl_mb(k_trial),maxnov)
        call ga_get(g_prod_p(i),ilo,ihi,1,ntrials,
     1    dbl_mb(k_prod),maxnov)
        do m=1,ntrials
          do n=1,ntrials
            do l=ilo,ihi
              apb(m,n)=apb(m,n)
     1          +dbl_mb(k_trial+(m-1)*maxnov+(l-ilo))
     2          *dbl_mb(k_prod +(n-1)*maxnov+(l-ilo))
            enddo
          enddo
        enddo
        if (.not.tda) then
          call ga_get(g_prod_m(i),ilo,ihi,1,ntrials,
     1      dbl_mb(k_prod),maxnov)
          do m=1,ntrials
            do n=1,ntrials
              do l=ilo,ihi
                amb(m,n)=amb(m,n)
     1            +dbl_mb(k_trial+(m-1)*maxnov+(l-ilo))
     2            *dbl_mb(k_prod +(n-1)*maxnov+(l-ilo))
              enddo
            enddo
          enddo
        endif
      enddo
      call ga_sync()
c
c -----------------------------------------------------------------
c Now we accumulate the local contributions to form the full matrix
c -----------------------------------------------------------------
c
      if (.not.ga_create(mt_dbl,ntrials,ntrials,'Work A+B',
     1  ntrials,ntrials,g_apb_tmp)) call errquit
     2  ('tddft_subspace: failed to create g_apb_tmp',0)
      call ga_zero(g_apb_tmp)
      call ga_acc(g_apb_tmp,1,ntrials,1,ntrials,apb,ntrials,1.0d0)
      if (.not.tda) then
        if (.not.ga_create(mt_dbl,ntrials,ntrials,'Work A-B',
     1    ntrials,ntrials,g_amb_tmp)) call errquit
     2    ('tddft_subspace: failed to create g_amb_tmp',0)
        call ga_zero(g_amb_tmp)
        call ga_acc(g_amb_tmp,1,ntrials,1,ntrials,amb,ntrials,1.0d0)
      endif
      call ga_sync()
      call ga_get(g_apb_tmp,1,ntrials,1,ntrials,apb,ntrials)
      if (.not.tda)
     1  call ga_get(g_amb_tmp,1,ntrials,1,ntrials,amb,ntrials)
      call ga_sync()
      if (.not.ga_destroy(g_apb_tmp)) call errquit
     1  ('tddft_subspace: failed to destroy g_apb_tmp',0)
      if (.not.tda) then
        if (.not.ga_destroy(g_amb_tmp)) call errquit
     1    ('tddft_subspace: failed to destroy g_amb_tmp',0)
      endif
c
c --------------
c Print matrices
c --------------
c
      if (util_print('subspace',print_debug)) then
        call ma_print(apb,ntrials,ntrials,'A+B')
        if (.not.tda) call ma_print(amb,ntrials,ntrials,'A-B')
      endif
c           
c ------
c Return
c ------
c
      if (.not.ma_pop_stack(l_prod)) call errquit
     1  ('tddft_subspace: failed to deallocate prod',0)
      if (.not.ma_pop_stack(l_trial)) call errquit
     1  ('tddft_subspace: failed to deallocate trial',0)
      return
      end
