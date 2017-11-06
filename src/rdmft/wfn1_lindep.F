c-----------------------------------------------------------------------
c
      subroutine wfn1_lindep(g_s,g_s12,g_v,nindep,tol)
      implicit none
C>
C> \brief Diagonalize the overlap matrix and find from the eigenvalues
C> how many linearly independent vectors there are
C>
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer g_s   ! the GA handle for the overlap matrix
      integer g_s12 ! the GA handle for the S^{-1/2} matrix
      integer g_v   ! the GA handle for the eigenvectors
c
      integer, intent(out)      :: nindep ! the number of linearly 
                                          ! independent vectors
      double precision, intent(in) :: tol ! the eigenvalue cutoff
c
      integer itype
      integer ndim
      integer idims(ga_max_dim)
      integer nbf       ! the number of basis functions
      integer ii        ! counter
      integer jj        ! counter
      integer kk        ! counter
      integer nproc     ! the number of processor ranks
      integer iproc     ! the current processor rank
      integer ilo,ihi   ! row limits on the local block
      integer jlo,jhi   ! column limits on the local block
      integer islo,ishi ! row limits on the local block
      integer jslo,jshi ! column limits on the local block
      integer k_v       ! memory location for the local block of v
      integer ldv       ! leading dimension of the local block of v
c
      double precision              :: dnrm
      double precision, allocatable :: eigval(:)
      double precision, allocatable :: vi(:,:)
      double precision, allocatable :: vj(:,:)
      double precision, allocatable :: s12(:,:)
c
      character*(12) pname
      pname = "wfn1_lindep:"
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      call nga_inquire(g_s,itype,ndim,idims)
      nbf = idims(1)
      if (ndim.ne.2) then
        call errquit(pname//" wrong number of dimensions",ndim,UERR)
      endif
      if (idims(1).ne.idims(2)) then
        call errquit(pname//" non-square overlap matrix",
     &               idims(2)-idims(1),UERR)
      endif
c
      allocate(eigval(1:nbf))
c
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
      call ga_pdsyev(g_s, g_v, eigval, 0)
#else
      call ga_diag_std(g_s, g_v, eigval)
#endif
#else
      call ga_diag_std_seq(g_s, g_v, eigval)
#endif
c
      nindep = 0
      do ii = 1, nbf
        if (eigval(ii).ge.tol) nindep = nindep + 1
      enddo
      if (nbf-nindep.gt.0.and.iproc.eq.0) then
        write(luout,'(" !! The overlap matrix has ",i5,
     &                " vectors deemed linearly dependent with"/
     &                "    eigenvalues:")')nbf-nindep
        write(luout,'(1p,8d9.2)')(eigval(ii),ii=0,nbf-nindep-1)
        write(luout,*)
      endif
c
c     compute S^{1/2} matrix
c
      call ga_distribution(g_s12,iproc,islo,ishi,jslo,jshi)
      allocate(vi(islo:ishi,nbf))
      allocate(vj(jslo:jshi,nbf))
      allocate(s12(islo:ishi,jslo:jshi))
      call ga_get(g_v,islo,ishi,1,nbf,vi,ishi-islo+1)
      call ga_get(g_v,jslo,jshi,1,nbf,vj,jshi-jslo+1)
      s12 = 0.0d0
      do kk = 1, nbf
        do jj = jslo, jshi
          do ii = islo, ishi
            s12(ii,jj) = s12(ii,jj)
     &                 + vi(ii,kk)*vj(jj,kk)*sqrt(eigval(kk))
          enddo ! ii
        enddo ! jj
      enddo ! kk
      call ga_put(g_s12,islo,ishi,jslo,jshi,s12,ishi-islo+1)
      deallocate(vi,vj,s12)
c
c     scale eigenvectors to obtain S^{-1/2}
c
      call ga_distribution(g_v,iproc,ilo,ihi,jlo,jhi)
      call ga_access(g_v,ilo,ihi,jlo,jhi,k_v,ldv)
      do jj = 0, jhi-jlo
        dnrm = 1.0d0/sqrt(eigval(jlo+jj))
        do ii = 0, ihi-ilo
          dbl_mb(k_v+jj*ldv+ii) = dbl_mb(k_v+jj*ldv+ii)*dnrm
        enddo
      enddo
c
      deallocate(eigval)
c
      end
c
c-----------------------------------------------------------------------
