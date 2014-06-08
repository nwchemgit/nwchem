C> \ingroup cholesky
C> @{
C>
C> \brief Test the Cholesky representation of the electron
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
C> This subroutine calculates \f$B\f$ by invoking the Cholesky routine.
C> The Cholesky decomposition is tested against the 4-dimension ERI
C> tensor.
C>
      subroutine util_test_cholesky(rtdb,geom,basis,tol,g_chol,k)
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
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
      integer g_eri  !< GA for the ERI tensor
c
      integer natoms !< The number of atoms
      integer mshlen !< Maximum number of basis function in a shell
      integer nshell !< The number of shells
      integer ish    !< Index of I-th shell
      integer jsh    !< Index of J-th shell
      integer ksh    !< Index of K-th shell
      integer lsh    !< Index of L-th shell
      integer ibf    !< Index of I-th basis function
      integer jbf    !< Index of J-th basis function
      integer kbf    !< Index of K-th basis function
      integer lbf    !< Index of L-th basis function
      integer mxbl   !< The maximum size of a shell pair
      integer m      !< The maximum number of Cholesky vectors
c
      integer max2e        !< The buffer size for ERI
      integer mem2         !< The scratch memory size
      integer l_eri, k_eri !< Buffer for ERI
      integer l_mem, k_mem !< Buffer for scratch memory
      integer l_row1, k_row1 !< Buffer for scratch memory, row of g_chol
      integer l_row2, k_row2 !< Buffer for scratch memory, row of g_chol
c
      integer nbf    !< The number of basis functions
c
      integer ndim             !< The number of dimensions of a GA
      integer dims(ga_max_dim) !< The dimensions of a GA
      integer chnk(ga_max_dim) !< The chunking of a GA
c
      integer lo_eri(ga_max_dim) !< The lower limits of blocks of eri
      integer hi_eri(ga_max_dim) !< The upper limits of blocks of eri
      integer ld_eri(ga_max_dim) !< The leading dims of blocks of eri
c
      integer lo_x(ga_max_dim) !< The lower limits of blocks
      integer hi_x(ga_max_dim) !< The upper limits of blocks
      integer ld_x(ga_max_dim) !< The leading dims of blocks
c
      integer itype               !< MA data type
      integer ioff                !< Local memory offset
      integer indx_b              !< Memory offset for Cholesky vectors
      integer indx_r              !< Memory offset for residual
      integer indx_d              !< Memory offset for diag
      integer indx_d0(ga_max_dim) !< The position max diag
      double precision val_d0     !< Max diag value
      double precision eri        !< Value of integral
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
      integer i, j !< Counters
c
      character*(17) pname
      parameter (pname="util_eri_cholesky")
c
      double precision ddot
      external         ddot
c
      iproc = ga_nodeid()
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
      call nga_inquire(g_chol,itype,ndim,dims)
      m = dims(3)
c
      ndim = 4
      dims(1) = nbf
      dims(2) = nbf
      dims(3) = nbf
      dims(4) = nbf
      chnk(1) = -1
      chnk(2) = -1
      chnk(3) = -1
      chnk(4) = -1
      if (.not.nga_create(MT_DBL,ndim,dims,"ERI R",chnk,g_eri))
     &  call errquit(pname//": ga_create ERI failed",
     &       ma_sizeof(MT_DBL,nbf*nbf*nbf*nbf,MT_BYTE),GA_ERR)
c
c     Check the data distribution
c
      call nga_distribution(g_eri,iproc,lo_eri,hi_eri)
c
      if (.not.ma_push_get(MT_DBL,max2e,pname//" ERI",l_eri,k_eri))
     &  call errquit(pname//": allocate ERI failed",
     &       ma_sizeof(MT_DBL,max2e,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,mem2,pname//" Scratch",l_mem,k_mem))
     &  call errquit(pname//": allocate Scratch failed",
     &       ma_sizeof(MT_DBL,mem2,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,m,pname//" Scratch row1",
     &                     l_row1,k_row1))
     &  call errquit(pname//": allocate Scratch Row1 failed",
     &       ma_sizeof(MT_DBL,m,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,m,pname//" Scratch row2",
     &                     l_row2,k_row2))
     &  call errquit(pname//": allocate Scratch Row2 failed",
     &       ma_sizeof(MT_DBL,m,MT_BYTE),MA_ERR)
c
c     Step B. Compute the ERI
c
      do ish = 1, nshell
        if (.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &     call errquit(pname//": bas_cn2bfr ish failed",
     &                  ish,BASIS_ERR)
        if ((lo_eri(1).le.ibflo).and.(ibflo.le.hi_eri(1))) then
          do jsh = 1, nshell
            if (.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &         call errquit(pname//": bas_cn2bfr jsh failed",
     &                      jsh,BASIS_ERR)
            if ((lo_eri(2).le.jbflo).and.(jbflo.le.hi_eri(2))) then
              do ksh = 1, nshell
                if (.not.bas_cn2bfr(basis,ksh,kbflo,kbfhi))
     &             call errquit(pname//": bas_cn2bfr ksh failed",
     &                          ksh,BASIS_ERR)
                if ((lo_eri(3).le.kbflo).and.(kbflo.le.hi_eri(3))) then
                  do lsh = 1, nshell
                    if (.not.bas_cn2bfr(basis,lsh,lbflo,lbfhi))
     &                 call errquit(pname//": bas_cn2bfr lsh failed",
     &                              lsh,BASIS_ERR)
                    if ((lo_eri(4).le.lbflo).and.(lbflo.le.hi_eri(4)))
     &              then
                      call int_2e4c(basis,lsh,ksh,basis,jsh,ish,
     &                              mem2,dbl_mb(k_mem),
     &                              max2e,dbl_mb(k_eri))
                      lo_x(1) = ibflo
                      lo_x(2) = jbflo
                      lo_x(3) = kbflo
                      lo_x(4) = lbflo
                      hi_x(1) = ibfhi
                      hi_x(2) = jbfhi
                      hi_x(3) = kbfhi
                      hi_x(4) = lbfhi
                      ld_x(1) = ibfhi-ibflo+1
                      ld_x(2) = jbfhi-jbflo+1
                      ld_x(3) = kbfhi-kbflo+1
                      ld_x(4) = lbfhi-lbflo+1
                      call nga_put(g_eri,lo_x,hi_x,dbl_mb(k_eri),ld_x)
                    endif
                  enddo ! lsh
                endif
              enddo ! ksh
            endif
          enddo ! jsh
        endif
      enddo ! ish
c
      call int_terminate
c
c     Calculate the Cholesky decomposition
c
      call util_eri_cholesky(rtdb,geom,basis,tol,g_chol,k)
      if (ga_nodeid().eq.0) then
        write(luout,*)'Number of Cholesky vectors = ',k
      endif
c
c     Subtract the reconstituted ERI from the originals
c
      kk = 0
      do ibf = 1, nbf
        do jbf = 1, nbf
          do kbf = 1, nbf
            do lbf = 1, nbf
              if (mod(kk,ga_nnodes()).eq.ga_nodeid()) then
                lo_x(1) = ibf
                lo_x(2) = jbf
                lo_x(3) = 1
                hi_x(1) = ibf
                hi_x(2) = jbf
                hi_x(3) = k
                ld_x(1) = 1
                ld_x(2) = 1
                call nga_get(g_chol,lo_x,hi_x,dbl_mb(k_row1),ld_x)
                lo_x(3) = kbf
                lo_x(4) = lbf
                lo_x(5) = 1
                hi_x(3) = kbf
                hi_x(4) = lbf
                hi_x(5) = k
                ld_x(1) = 1
                ld_x(2) = 1
                ld_x(3) = 1
                ld_x(4) = 1
                call nga_get(g_chol,lo_x(3),hi_x(3),dbl_mb(k_row2),
     &                       ld_x(3))
                eri = ddot(k,dbl_mb(k_row1),1,dbl_mb(k_row2),1)
                call nga_acc(g_eri,lo_x,hi_x,eri,ld_x,-1.0d0)
              endif
              kk = kk + 1
            enddo ! lbf
          enddo ! kbf
        enddo ! jbf
      enddo ! ibf
c
c     Step C. Find the coordinates of the maximum element of the
c             diagonal.
c
      call nga_idamax(g_eri,indx_d0,val_d0)
c
c     Check result
c
      if (val_d0.gt.tol) then
        write(*,*)'*** ERROR: indx_d0 = ',indx_d0(1:4)
        write(*,*)'*** ERROR: val_d0  = ',val_d0
      endif
c
      if (.not.ma_pop_stack(l_row2))
     &  call errquit(pname//": deallocate Scratch Row2 failed",0,MA_ERR)
      if (.not.ma_pop_stack(l_row1))
     &  call errquit(pname//": deallocate Scratch Row1 failed",0,MA_ERR)
      if (.not.ma_pop_stack(l_mem))
     &  call errquit(pname//": deallocate Scratch failed",0,MA_ERR)
      if (.not.ma_pop_stack(l_eri))
     &  call errquit(pname//": deallocate ERI failed",0,MA_ERR)
c
      if (.not.ga_destroy(g_eri))
     &  call errquit(pname//": ga_destroy ERI failed",0,GA_ERR)
c
      end
C>
C> @}
c $Id$
