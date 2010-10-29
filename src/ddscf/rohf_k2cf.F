      subroutine rohf_k2cf(basis, nbf, nmo, nclosed, nopen,
     $                     g_kvec, g_vecs, g_result)
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "cscfps.fh"
c     
c     Transform mo-coeffs from k-vector using exp(K)
c     
c     Arguments
c
      integer nbf, basis, nmo, nclosed, nopen
      integer g_kvec, g_vecs, g_result
c     
      integer g_temp
      double precision one, zero
      parameter (one=1.0d0, zero=0.0d0)
c     
      if (oscfps) call pstat_on(ps_k2cf)
c
      if (.not. ga_create(MT_DBL, nmo, nmo, 'k2cf: temp',
     $     32, 32, g_temp)) call errquit('rohf_k2cf: GA temp',0, GA_ERR)
      call rohf_k2u( basis, nbf, nmo, nclosed, nopen,
     $               g_kvec, g_temp )
c
c     Transform MO's
c     
      call ga_dgemm('n', 'n', nbf, nmo, nmo, one, g_vecs,
     $               g_temp, zero, g_result)
c
      if (.not.ga_destroy(g_temp))
     $     call errquit('rohf_k2cf: cannot destroy temp',0, GA_ERR)
c
*      call ga_orthog_mos(basis, g_result)
c
      if (oscfps) call pstat_off(ps_k2cf)
c
      end
      subroutine rohf_k2u( basis, nbf, nmo, nclosed, nopen,
     $                     g_kvec, g_u )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
c  
c  Generate U = exp(K)
c  K is the antisymmetric matrix form of k
c
c     
c     Arguments
c
      integer nbf, basis, nmo, nclosed, nopen
      integer g_kvec, g_u
c     
      integer nvir
      integer i,l_tmp,k_tmp,iioff,ivoff,nnnodes
      integer tmp_siz
      double precision one, mone, zero
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0)
      double precision kmax
c     
      call ga_sync()
      nnnodes = ga_nnodes()
c     
c     Form K-matrix from vector
c     
      nvir = nmo - nclosed - nopen
      ivoff = nclosed + nopen + 1
      tmp_siz =  max(nvir,nopen)
c
c
      if (.not. ma_push_get(MT_DBL,tmp_siz,'temp k',l_tmp,k_tmp))
     $     call errquit('k2cf: ma failed on tmp',tmp_siz, MA_ERR)
c
      call ga_zero(g_u)
      if (nvir .gt. 0) then
         do i=ga_nodeid()+1,nclosed+nopen,nnnodes
            iioff = (i-1)*nvir
            call ga_get(g_kvec,iioff+1,iioff+nvir,1,1,
     $           dbl_mb(k_tmp),nvir)
            call ga_put(g_u,ivoff,nmo,i,i,dbl_mb(k_tmp),nvir)
            call dscal(nvir,mone,dbl_mb(k_tmp),1)
            call ga_put(g_u,i,i,ivoff,nmo,dbl_mb(k_tmp),1)
         enddo
      endif
      if (nopen .gt. 0) then
         do i=ga_nodeid()+1,nclosed,nnnodes
            iioff = (nclosed+nopen)*nvir + (i-1)*nopen
            call ga_get(g_kvec,iioff+1,iioff+nopen,1,1,
     $           dbl_mb(k_tmp),nopen)
            call ga_put(g_u,(nclosed+1),(nclosed+nopen),i,i,
     $           dbl_mb(k_tmp),nopen)
            call dscal(nopen,mone,dbl_mb(k_tmp),1)
            call ga_put(g_u,i,i,(nclosed+1),(nclosed+nopen),
     $           dbl_mb(k_tmp),1)
         enddo
      endif
      if (.not. ma_pop_stack(l_tmp))
     $     call errquit('k2cf: pop failed', 0, MA_ERR)
c
c     Make near zeoes exactly zero so that sparsity tests in the
c     matrix multiply work well ... we only need to maintain
c     quadratic convergence.
c     Cannot screen with a large number (greater than accuracy in integrals) 
c     since in a non-abelian group we could break symmetry ... either
c     screen with a small threshold or screen with a large one and then
c     enforce symmetry (which is expensive in high order groups).
c     Also, screening must tend to zero for high convergence
c
      call ga_maxelt(g_u, kmax)
      kmax = min(kmax*1d-2,kmax*kmax,1d-12)
      call ga_screen(g_u, kmax)
c     
c     Get a unitary 2nd order approximation to exp(K)
c     
      call matrix_exp(g_u )
      end
