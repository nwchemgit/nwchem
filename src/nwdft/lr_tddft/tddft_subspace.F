      subroutine tddft_subspace(ipol,nocc,nmo,nfc,nfv,g_trials,
     1  g_prod_p,g_prod_m,g_buffer1,g_buffer2,d_trials,d_prod_p,
     2  d_prod_m,ntrials,tda,apb,amb,algorithm)
c
c $Id$
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
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_p(2)      ! GA handle for A+B product vectors
      integer g_prod_m(2)      ! GA handle for A-B product vectors
      integer g_buffer1(2)     ! GA handle for buffer vectors
      integer g_buffer2(2)     ! GA handle for buffer vectors
      integer d_trials(2)      ! DRA handle for trial vectors
      integer d_prod_p(2)      ! DRA handle for A+B product vectors
      integer d_prod_m(2)      ! DRA handle for A-B product vectors
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      integer algorithm        ! Algorithm
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
      integer dummy,request
      integer g_ab
      double precision mynorm
c
c     Determine the length of the trial vector
      call tddft_lentrialvec(ipol,nocc,nmo,nfc,nfv,nov)
c
c =================
c Incore algorithms
c =================
c
      if ((algorithm.eq.1).or.(algorithm.eq.2)) then
        iproc=ga_nodeid()
        maxnov=0
        do i=1,ipol
          call ga_distribution(g_trials(i),iproc,ilo,ihi,jlo,jhi)
          maxnov=max(maxnov,ihi-ilo+1)
        enddo
        if (.not.ma_push_get(mt_dbl,maxnov*ntrials,'trial',
     1    l_trial,k_trial)) call errquit
     2    ('tddft_subspace: failed to allocate trial',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,maxnov*ntrials,'product',
     1    l_prod,k_prod)) call errquit
     2    ('tddft_subspace: failed to allocate prod',0, MA_ERR)
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
#if 0
c
        call ga_sync()
        do i=1,ipol
         call ga_distribution(g_trials(i),iproc,ilo,ihi,jlo,jhi)
         if ((jlo.le.1).and.(jhi.ge.ntrials)) then
          call ga_get(g_trials(i),ilo,ihi,1,ntrials,
     1      dbl_mb(k_trial),maxnov)
          call ga_get(g_prod_p(i),ilo,ihi,1,ntrials,
     1      dbl_mb(k_prod),maxnov)
          do m=1,ntrials
            do n=1,ntrials
              do l=ilo,ihi
                apb(m,n)=apb(m,n)
     1            +dbl_mb(k_trial+(m-1)*maxnov+(l-ilo))
     2            *dbl_mb(k_prod +(n-1)*maxnov+(l-ilo))
              enddo
            enddo
          enddo
          if (.not.tda) then
            call ga_get(g_prod_m(i),ilo,ihi,1,ntrials,
     1        dbl_mb(k_prod),maxnov)
            do m=1,ntrials
              do n=1,ntrials
                do l=ilo,ihi
                  amb(m,n)=amb(m,n)
     1              +dbl_mb(k_trial+(m-1)*maxnov+(l-ilo))
     2              *dbl_mb(k_prod +(n-1)*maxnov+(l-ilo))
                enddo
              enddo
            enddo
          endif ! .not. tda
         endif ! jlo,jhi check
        enddo
c
c -----------------------------------------------------------------
c Now we accumulate the local contributions to form the full matrix
c -----------------------------------------------------------------
c
c use gdops for amb and apb instead of costly accs        
        call ga_dgop(913,apb,ntrials*ntrials, '+')
        if (.not.tda)
     C       call ga_dgop(914,amb,ntrials*ntrials, '+')
        if(ga_nodeid().eq.0) write(6,*) apb
#else
        if (.not.ga_create(mt_dbl,ntrials,ntrials,
     1       'g_ab',0,0,g_ab))
     2       call errquit('tddft_subspace: failed to create g_amb',0,
     &       GA_ERR)
        call ga_zero(g_ab)
        do i=1,ipol
           call ga_dgemm('t', 'n', ntrials, ntrials, nov(i), 1.d0,
     D          g_trials(i),g_prod_p(i),1d0,g_ab)
        enddo
c        call ga_print(g_ab)
        call util_mygabcast(g_ab,ntrials,ntrials,apb,ntrials)
        if (.not.tda) then
           call ga_zero(g_ab)
           do i=1,ipol
              call ga_dgemm('t', 'n', ntrials, ntrials, nov(i), 1.d0,
     D             g_trials(i),g_prod_m(i),1d0,g_ab)
           enddo
c           call ga_print(g_ab)
           call util_mygabcast(g_ab,ntrials,ntrials,amb,ntrials)
        endif
        if (.not.ga_destroy(g_ab))
     &    call errquit('ga_dest g_ab fail',0, GA_ERR)
#endif
        if (.not.ma_pop_stack(l_prod)) call errquit
     1    ('tddft_subspace: failed to deallocate prod',0, MA_ERR)
        if (.not.ma_pop_stack(l_trial)) call errquit
     1    ('tddft_subspace: failed to deallocate trial',0, MA_ERR)
c
c ===============
c Disk algorithms
c ===============
c This part could be highly inefficient -- consider saving previously
c computed blocks of A+B and A-B matrices.
c
      else
        do m=1,ntrials
          do n=1,ntrials
            apb(m,n)=0.0d0
            if (.not.tda) amb(m,n)=0.0d0
          enddo
        enddo
        do i=1,ipol
          do m=1,ntrials
            dummy=dra_read_section(.false.,g_buffer1(i),1,nov(i),1,1,
     1        d_trials(i),1,nov(i),m,m,request)
            dummy=dra_wait(request)
            do n=1,ntrials
              dummy=dra_read_section(.false.,g_buffer2(i),1,nov(i),1,1,
     1          d_prod_p(i),1,nov(i),n,n,request)
              dummy=dra_wait(request)
              apb(m,n)=apb(m,n)+ga_ddot(g_buffer1(i),g_buffer2(i))
            enddo
            if (.not.tda) then
              do n=1,ntrials
                dummy=dra_read_section(.false.,
     1            g_buffer2(i),1,nov(i),1,1,
     2            d_prod_m(i),1,nov(i),n,n,request)
                dummy=dra_wait(request)
                amb(m,n)=amb(m,n)+ga_ddot(g_buffer1(i),g_buffer2(i))
              enddo
            endif
          enddo
        enddo
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
      return
      end
