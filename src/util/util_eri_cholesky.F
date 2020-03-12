C> \ingroup cholesky
C> @{
C>
C> \brief Calculate the Cholesky representation of the electron
C> repulsion integrals
C>
C> This subroutine calculates the Cholesky decomposition of the 
C> 2-electron repulsion integral (ERI) tensor. The ERI are provided
C> in the atomic orbital basis as
C> \f{eqnarray*}{
C>   \mathrm{eri}(\mu,\rho,\nu,\lambda) &=& (\mu,\rho|\nu,\lambda) \\\\
C>   &=& \int_{-\infty}^{\infty}
C>    g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}
C>    g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
C> \f}
C> where \f$\mu\f$, \f$\nu\f$, \f$\rho\f$, and \f$\lambda\f$ label
C> spatial Gaussian atomic orbitals \f$g\f$. The quantities \f$X\f$ 
C> label the atomic coordinates on which the obitals are centered.
C>
C> The Cholesky decomposition approximates these integrals as
C> \f{eqnarray*}{
C>   (\mu,\rho|\nu,\lambda) &\approx &
C>        \sum_{L=1}^M (\mu,\rho|L)(L|\nu,\lambda) \\\\
C>   &=& \sum_{L=1}^M B_{\mu\rho}^L B_{\nu\lambda}^L
C> \f}
C> This subroutine calculates \f$B\f$ in an integral direct way. I.e.
C> the ERI are recalculated whenever needed rather than stored (with
C> the exception of the diagonal of the ERI tensor).
C>
C> The algorithm to compute the decomposition was taken from [1].
C> It can be summarized as:
C>
C> - A. Set \f$k = 0\f$
C> - B. Compute the ERI diagonal
C>   \f$D^0_{\nu\lambda,\nu\lambda} = (\nu\lambda|\nu\lambda)\f$
C> - C. Choose compound index \f$\nu\lambda_k\f$ such that
C>   \f$D^k_{\nu\lambda_k,\nu\lambda_k}\f$ is maximal
C> - D. Do While (\f$\max\left(D^{k}_{\nu\lambda,\nu\lambda}\right) 
C>   > \delta\f$)
C>   - E. Compute all integrals \f$R^k_{\mu\rho,\nu\lambda_k} = 
C>     (\mu,\rho|\nu\lambda_k)\f$
C>   - F. Compute \f$R^k_{\mu\rho,\nu\lambda_k}
C>     = R^k_{\mu\rho,\nu\lambda_k}
C>     - \sum_{i=1}^{k-1}B_{\mu\rho}^i B_{\nu\lambda_k}^i\f$
C>   - G. Compute Cholesky vector \f$B_{\mu\rho}^{k+1} = 
C>     R^k_{\mu\rho,\nu\lambda_k} /
C>     \sqrt{D^k_{\nu\lambda_k,\nu\lambda_k}}\f$
C>   - H. \f$k=k+1\f$
C>   - I. Compute \f$D^k_{\nu\lambda,\nu\lambda} = 
C>     D^{k-1}_{\nu\lambda,\nu\lambda} -
C>     B^k_{\nu\lambda}B^k_{\nu\lambda}\f$
C>   - J. Choose compound index \f$\nu\lambda_k\f$ such that
C>     \f$D^k_{\nu\lambda_k,\nu\lambda_k}\f$ is maximal
C> - K. End Do
C>
C> This algorithm requires \f$D\f$ and \f$R\f$ as temporary variables.
C> The most appropriate data distribution is to distribute the vectors.
C> I.e. if we indicate the indeces being partitioned by square brackets
C> then the distributions are \f$D_{[\nu][\lambda],[\nu][\lambda]}\f$,
C> \f$R_{[\mu][\rho],\nu\lambda_k}\f$ and we assume 
C> \f$B^k_{[\mu][\rho]}\f$.
C>
C> Clearly the algorithm involves 3 distributed data structures \f$B\f$,
C> \f$D\f$, and \f$R\f$. All 3 data structures are stored in global
C> arrays. These arrays have the following properties:
C> - \f$B^k_{[\mu][\rho]}\f$:
C>   - A 3-dimensional global array
C>   - The 1st dimension runs over basis functions and is distributed
C>     over processors
C>   - The 2nd dimension runs over basis functions and is distributed
C>     over processors
C>   - The 3rd dimension runs over Cholesky vectors and is not 
C>     distributed
C> - \f$D_{[\mu][\rho],[\mu][\rho]}\f$:
C>   - A 2-dimensional global array representing the diagonal of a
C>     4-dimensional tensor
C>   - The 1st dimension runs over basis functions and is distributed
C>     over processors, and the distribution has to be the same as
C>     the distribution of the 1st dimension of \f$B^k_{[\mu][\rho]}\f$.
C>   - The 2nd dimension runs over basis functions and is distributed
C>     over processors, and the distribution has to be the same as
C>     the distribution of the 2nd dimension of \f$B^k_{[\mu][\rho]}\f$.
C> - \f$R_{[\mu][\rho],\nu\lambda_k}\f$:
C>   - A 2-dimensional global array representing \f$R\f$ of all
C>     \f$\mu\f$ and \f$\rho\f$ for a single fixed value of the 
C>     compound index \f$\nu\lambda_k\f$.
C>   - The 1st dimension runs over basis functions and is distributed
C>     over processors, and the distribution has to be the same as
C>     the distribution of the 1st dimension of \f$B^k_{[\mu][\rho]}\f$.
C>   - The 2nd dimension runs over basis functions and is distributed
C>     over processors, and the distribution has to be the same as
C>     the distribution of the 2nd dimension of \f$B^k_{[\mu][\rho]}\f$.
C>
C> Currently the algorithm does not explicitly ensure that the data
C> distribution of \f$D\f$ and \f$R\f$ match the one of \f$B\f$. The
C> algorithm only tests whether the distributions match and if they do
C> not it aborts. Of course if needed the algorithm can be enhanced to
C> take the data distribution of \f$B\f$ and use that in the creation of
C> \f$D\f$ and \f$R\f$.
C>
C> Comment: NWChem evaluates the integrals in shell blocks. So depending
C> on which \f$\nu\lambda_k\f$ is picked a multiple of the required
C> number of integrals is calculated. E.g. if \f$\nu\f$ and 
C> \f$\lambda\f$ are in p-shells then to evaluate the \f$\nu\lambda_k\f$
C> integrals all 9 shell-pair components have to be evaluated.
C> Initially we will not exploit this but for efficiency reasons we
C> will have to consider this at some point.
C> 
C> ### References ###
C>
C> [1] E. Epifanovsky, D. Zuev, X. Feng, K. Khistyaev, Y. Shao,
C>     A. I. Krylov, "General implementation of the 
C>     resolution-of-the-identity and Cholesky representations of 
C>     electron repulsion integral within coupled-cluster and
C>     equation-of-motion methods: Theory and benchmarks",
C>     J. Chem. Phys. <b>139</b> (2013) 134105, DOI:
C>     <a href="https://doi.org/10.1063/1.4820484">
C>     10.1063/1.4820484</a>.
C>
      subroutine util_eri_cholesky(rtdb,geom,basis,tol,g_chol,k)
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
c
c     Arguments
c
      integer rtdb   !< [Input] The RTDB handle
      integer geom   !< [Input] The geometry handle
      integer basis  !< [Input] The basis set handle
c
      double precision tol !< [Input] The tolerance \f$\delta\f$ in the
                           !< algorithm outline
c
      integer g_chol !< [Output] The Cholesky vectors
                     !< \f$B_{[\nu\lambda]}^{1:M}\f$
      integer k      !< [Output] The number of Cholesky vectors
c
c     Local
c
      integer iproc  !< The rank of this process
c
      integer kk     !< Counter
      integer g_d    !< GA for the diagonal of the ERI tensor
      integer g_r    !< GA for the residual
      integer m      !< The maximum number of Cholesky vectors
c
      integer natoms !< The number of atoms
      integer mshlen !< Maximum number of basis function in a shell
      integer nshell !< The number of shells
      integer ish    !< Index of I-th shell
      integer jsh    !< Index of J-th shell
      integer ksh    !< Index of K-th shell
      integer lsh    !< Index of L-th shell
      integer mxbl   !< The maximum size of a shell pair
c
      integer max2e        !< The buffer size for ERI
      integer mem2         !< The scratch memory size
      integer l_eri, k_eri !< Buffer for ERI
      integer l_mem, k_mem !< Buffer for scratch memory
      integer l_row, k_row !< Buffer for scratch memory, row of g_chol
c
      integer nbf    !< The number of basis functions
c
      integer ndim             !< The number of dimensions of a GA
      integer dims(ga_max_dim) !< The dimensions of a GA
      integer chnk(ga_max_dim) !< The chunking of a GA
c
      integer lo_b(ga_max_dim) !< The lower limits of blocks of B
      integer lo_r(ga_max_dim) !< The lower limits of blocks of R
      integer lo_d(ga_max_dim) !< The lower limits of blocks of D
      integer hi_b(ga_max_dim) !< The upper limits of blocks of B
      integer hi_r(ga_max_dim) !< The upper limits of blocks of R
      integer hi_d(ga_max_dim) !< The upper limits of blocks of D
      integer ld_b(ga_max_dim) !< The leading dims of blocks of B
      integer ld_r(ga_max_dim) !< The leading dims of blocks of R
      integer ld_d(ga_max_dim) !< The leading dims of blocks of D
c
      integer lo_x(ga_max_dim) !< The lower limits of blocks
      integer hi_x(ga_max_dim) !< The upper limits of blocks
      integer ld_x(ga_max_dim) !< The leading dims of blocks
c
      integer itype               !< MA data type
      integer ioff                !< Local memory offset
      integer ijoff               !< Local memory offset
      integer ijijoff             !< Local memory offset
      integer indx_b              !< Memory offset for Cholesky vectors
      integer indx_r              !< Memory offset for residual
      integer indx_d              !< Memory offset for diag
      integer indx_d0(ga_max_dim) !< The position max diag
      double precision val_d0     !< Max diag value
c
      integer ibflo !< The lowest basis function of a shell
      integer ibfhi !< The highest basis function of a shell
      integer jbflo !< The lowest basis function of a shell
      integer jbfhi !< The highest basis function of a shell
      integer kbflo !< The lowest basis function of a shell
      integer kbfhi !< The highest basis function of a shell
      integer lbflo !< The lowest basis function of a shell
      integer lbfhi !< The highest basis function of a shell
      integer ibfnm !< The number of basis functions in i
      integer jbfnm !< The number of basis functions in j
      integer kbfnm !< The number of basis functions in k
      integer lbfnm !< The number of basis functions in l
      integer ijnm  !< The number of i j basis function pairs
      integer ijknm !< The number of i j k basis function triples
c
      integer i, j  !< Counters
c
      character*(17) pname
      parameter (pname="util_eri_cholesky")
c
      iproc = ga_nodeid()
      call ga_zero(g_chol)
      call int_init(rtdb,1,basis)
c
      if (.not.geom_ncent(geom, natoms)) 
     &  call errquit(pname//": geom_ncent failed",0,GEOM_ERR)
      if (.not.bas_numcont(basis, nshell))
     &  call errquit(pname//": bas_numcont failed",0,BASIS_ERR)
      if (.not.bas_nbf_cn_max(basis,mshlen))
     &  call errquit(pname//": bas_nbf_cn_max failed",0,BASIS_ERR)
      if (.not.bas_numbf(basis,nbf))
     &  call errquit(pname//": bas_numbf failed",0,BASIS_ERR)
c
      call int_mem_2e4c(max2e, mem2)
c
      mxbl  = mshlen*mshlen
c
      if (.not.ga_create(MT_DBL,nbf,nbf,"ERI diag",-1,-1,g_d))
     &  call errquit(pname//": ga_create ERI diag failed",
     &       ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),GA_ERR)
      ndim = 2
      dims(1) = nbf
      dims(2) = nbf
c     dims(3) = mxbl
      chnk(1) = -1
      chnk(2) = -1
c     chnk(3) = mxbl
      if (.not.nga_create(MT_DBL,ndim,dims,"ERI R",chnk,g_r))
     &  call errquit(pname//": ga_create ERI R failed",
     &       ma_sizeof(MT_DBL,nbf*nbf*1,MT_BYTE),GA_ERR)
c
c     Check the data distribution
c
      call nga_inquire(g_chol,itype,ndim,dims)
      m = dims(3)
      call nga_distribution(g_chol,iproc,lo_b,hi_b)
      call nga_distribution(g_d,iproc,lo_d,hi_d)
      call nga_distribution(g_r,iproc,lo_r,hi_r)
      if ((lo_b(1).eq.lo_d(1).and.lo_b(1).eq.lo_r(1)).and.
     &    (hi_b(1).eq.hi_d(1).and.hi_b(1).eq.hi_r(1)).and.
     &    (lo_b(2).eq.lo_d(2).and.lo_b(2).eq.lo_r(2)).and.
     &    (hi_b(2).eq.hi_d(2).and.hi_b(2).eq.hi_r(2))) then
c       data distribution is as expected
      else
        call errquit(pname//": unexpected data distribution",0,GA_ERR)
      endif
c
      if (.not.ma_push_get(MT_DBL,max2e,pname//" ERI",l_eri,k_eri))
     &  call errquit(pname//": allocate ERI failed",
     &       ma_sizeof(MT_DBL,max2e,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,mem2,pname//" Scratch",l_mem,k_mem))
     &  call errquit(pname//": allocate Scratch failed",
     &       ma_sizeof(MT_DBL,mem2,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,m,pname//" Scratch row",l_row,k_row))
     &  call errquit(pname//": allocate Scratch Row failed",
     &       ma_sizeof(MT_DBL,m,MT_BYTE),MA_ERR)
c
c     Step A. Initialize k
c
      k = 0
c
c     Step B. Compute the diagonal
c
      do ish = 1, nshell
        if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &     call errquit(pname//": bas_cn2bfr ish failed",
     &                  ish,BASIS_ERR)
        if ((lo_d(1).le.ibflo).and.(ibflo.le.hi_d(1))) then
          do jsh = 1, nshell
            if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &         call errquit(pname//": bas_cn2bfr jsh failed",
     &                      jsh,BASIS_ERR)
            if ((lo_d(2).le.jbflo).and.(jbflo.le.hi_d(2))) then
              call int_2e4c(basis,jsh,ish,basis,jsh,ish,
     &                      mem2,dbl_mb(k_mem),max2e,dbl_mb(k_eri))
c
c             We only need the (i,j|i,j) integrals
c             so compress the buffer
c
              ibfnm = ibfhi-ibflo+1
              jbfnm = jbfhi-jbflo+1
              kbfnm = ibfhi-ibflo+1
              lbfnm = jbfhi-jbflo+1
              ijnm  = ibfnm*jbfnm
              ijknm = ijnm*kbfnm
              do i = 0, ibfnm-1
                do j = 0, jbfnm-1
                  ijoff = j*ibfnm+i
                  ijijoff = ijoff+i*ijnm+j*ijknm
                  dbl_mb(k_eri+ijoff) = dbl_mb(k_eri+ijijoff)
                enddo
              enddo
              call ga_put(g_d,ibflo,ibfhi,jbflo,jbfhi,dbl_mb(k_eri),
     &                    ibfnm)
            endif
          enddo ! jsh
        endif
      enddo ! ish
c
c     Step C. Find the coordinates of the maximum element of the
c             diagonal.
c
      call nga_idamax(g_d,indx_d0,val_d0)
c
c     Step D. Start the while loop
c
      do while(val_d0.gt.tol.and.k.lt.m)
c
c       Find the corresponding shells (k,l)
c
        if (.not.bas_bf2cn(basis,indx_d0(1),ksh)) 
     &     call errquit(pname//": bas_bf2cn failed on ksh",indx_d0(1),
     &                  BASIS_ERR)
        if (.not.bas_bf2cn(basis,indx_d0(2),lsh)) 
     &     call errquit(pname//": bas_bf2cn failed on lsh",indx_d0(2),
     &                  BASIS_ERR)
        if (.not.bas_cn2bfr(basis,ksh,kbflo,kbfhi))
     &     call errquit(pname//": bas_cn2bfr ksh failed",
     &                  ksh,BASIS_ERR)
        if (.not.bas_cn2bfr(basis,lsh,lbflo,lbfhi))
     &     call errquit(pname//": bas_cn2bfr lsh failed",
     &                  lsh,BASIS_ERR)
c
c       Step E. Compute the integrals (i,j|k,l)
c
        call ga_zero(g_r)
        do ish = 1, nshell
          if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &       call errquit(pname//": bas_cn2bfr ish failed",
     &                    ish,BASIS_ERR)
          if ((lo_r(1).le.ibflo).and.(ibflo.le.hi_r(1))) then
            do jsh = 1, nshell
              if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &           call errquit(pname//": bas_cn2bfr jsh failed",
     &                        jsh,BASIS_ERR)
              if ((lo_r(2).le.jbflo).and.(jbflo.le.hi_r(2))) then
                call int_2e4c(basis,lsh,ksh,basis,jsh,ish,
     &                        mem2,dbl_mb(k_mem),max2e,dbl_mb(k_eri))
                ibfnm = ibfhi-ibflo+1
                jbfnm = jbfhi-jbflo+1
                kbfnm = kbfhi-kbflo+1
                lbfnm = lbfhi-lbflo+1
                ijnm  = ibfnm*jbfnm
                ijknm = ijnm*kbfnm
                ioff  = (indx_d0(2)-lbflo)*ijknm+(indx_d0(1)-kbflo)*ijnm
                lo_x(1) = ibflo
                lo_x(2) = jbflo
                hi_x(1) = ibfhi
                hi_x(2) = jbfhi
                ld_x(1) = ibfhi-ibflo+1
                call nga_put(g_r,lo_x,hi_x,dbl_mb(k_eri+ioff),ld_x)
              endif
            enddo ! jsh
          endif
        enddo ! ish
        call ga_sync()
c
c       Step F. Update the residual
c
        lo_x(1) = indx_d0(1)
        lo_x(2) = indx_d0(2)
        lo_x(3) = 1
        hi_x(1) = indx_d0(1)
        hi_x(2) = indx_d0(2)
        hi_x(3) = max(k,1)
        ld_x(1) = 1
        ld_x(2) = 1
        call nga_get(g_chol,lo_x,hi_x,dbl_mb(k_row),ld_x)
        call nga_access(g_r,lo_r,hi_r,indx_r,ld_r)
        call nga_access(g_chol,lo_b,hi_b,indx_b,ld_b)
        do kk = 0, k-1
          do i = 0, hi_r(1)-lo_r(1)
            do j = 0, hi_r(2)-lo_r(2)
              dbl_mb(indx_r+i+j*ld_r(1)) = dbl_mb(indx_r+i+j*ld_r(1))
     &              - dbl_mb(indx_b+i+j*ld_b(1)+kk*ld_b(2)*ld_b(1))
     &              * dbl_mb(k_row+kk)
            enddo
          enddo
        enddo
        call nga_release(g_chol,lo_b,hi_b)
        call nga_release_update(g_r,lo_r,hi_r)
c
c       Step G. Compute the new Cholesky vector
c
        call nga_access(g_r,lo_r,hi_r,indx_r,ld_r)
        call nga_access(g_chol,lo_b,hi_b,indx_b,ld_b)
        do i = 0, hi_r(1)-lo_r(1)
          do j = 0, hi_r(2)-lo_r(2)
            dbl_mb(indx_b+i+j*ld_b(1)+k*ld_b(2)*ld_b(1))
     &      = dbl_mb(indx_r+i+j*ld_r(1))/sqrt(val_d0)
          enddo
        enddo
        call nga_release_update(g_chol,lo_b,hi_b)
        call nga_release(g_r,lo_r,hi_r)
c
c       Step H. Increment k
c
        k = k + 1
c
c       Step I. Update the diagonal
c
        call nga_access(g_d,lo_d,hi_d,indx_d,ld_d)
        call nga_access(g_chol,lo_b,hi_b,indx_b,ld_b)
        do i = 0, hi_d(1)-lo_d(1)
          do j = 0, hi_d(2)-lo_d(2)
            dbl_mb(indx_d+i+j*ld_d(1)) = dbl_mb(indx_d+i+j*ld_d(1))
     &            - dbl_mb(indx_b+i+j*ld_b(1)+(k-1)*ld_b(2)*ld_b(1))**2
          enddo
        enddo
        call nga_release(g_chol,lo_b,hi_b)
        call nga_release_update(g_d,lo_d,hi_d)
c
c       Step J. Find the coordinates of the maximum element of the
c               diagonal.
c
        call nga_idamax(g_d,indx_d0,val_d0)
c
      enddo
c
      call int_terminate
c
      if (.not.ma_pop_stack(l_row))
     &  call errquit(pname//": deallocate Scratch Row failed",0,MA_ERR)
      if (.not.ma_pop_stack(l_mem))
     &  call errquit(pname//": deallocate Scratch failed",0,MA_ERR)
      if (.not.ma_pop_stack(l_eri))
     &  call errquit(pname//": deallocate ERI failed",0,MA_ERR)
c
      if (.not.ga_destroy(g_r))
     &  call errquit(pname//": ga_destroy ERI R failed",0,GA_ERR)
      if (.not.ga_destroy(g_d))
     &  call errquit(pname//": ga_destroy ERI diag failed",0,GA_ERR)
c
      end
C>
C> @}
c $Id$
