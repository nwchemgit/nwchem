      subroutine rohf_k2cf(basis, nbf, nclosed, nopen,
     $     g_kvec, g_vecs, g_result)
C     $Id: rohf_k2cf.F,v 1.1 1995-03-31 01:36:10 d3g681 Exp $
      implicit none
#include "mafdecls.h"
#include "global.fh"
c     
c     Transform mo-coeffs from k-vector using exp(K)
c     
c     Arguments
c
      integer nbf, basis, nclosed, nopen
      integer g_kvec, g_vecs, g_result
c     
      integer g_temp
      integer nvir
      integer i,l_tmp,k_tmp,iioff,ivoff,nnnodes
      integer tmp_siz
      double precision one, mone, zero
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0)
c     
      call ga_sync()
      nnnodes = ga_nnodes()
c     
c     Form K-matrix from vector
c     
      nvir = nbf - nclosed - nopen
      ivoff = nclosed + nopen + 1
      tmp_siz =  max(nvir,nopen)
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'k2cf: temp',
     $     1, 1, g_temp)) call errquit('rohf_k2cf: GA temp',0)
c
      if (.not. ma_push_get(MT_DBL,tmp_siz,'temp k',l_tmp,k_tmp))
     $     call errquit('k2cf: ma failed on tmp',tmp_siz)
c
      call ga_zero(g_temp)
      do i=ga_nodeid()+1,nclosed+nopen,nnnodes
         iioff = (i-1)*nvir
         call ga_get(g_kvec,iioff+1,iioff+nvir,1,1,dbl_mb(k_tmp),nvir)
         call ga_put(g_temp,ivoff,nbf,i,i,dbl_mb(k_tmp),nvir)
         call dscal(nvir,mone,dbl_mb(k_tmp),1)
         call ga_put(g_temp,i,i,ivoff,nbf,dbl_mb(k_tmp),1)
      enddo
      if (nopen .gt. 0) then
         do i=ga_nodeid()+1,nclosed,nnnodes
            iioff = (nclosed+nopen)*nvir + (i-1)*nopen
            call ga_get(g_kvec,iioff+1,iioff+nopen,1,1,
     $           dbl_mb(k_tmp),nopen)
            call ga_put(g_temp,(nclosed+1),(nclosed+nopen),i,i,
     $           dbl_mb(k_tmp),nopen)
            call dscal(nopen,mone,dbl_mb(k_tmp),1)
            call ga_put(g_temp,i,i,(nclosed+1),(nclosed+nopen),
     $           dbl_mb(k_tmp),1)
         enddo
      endif
      if (.not. ma_pop_stack(l_tmp))
     $     call errquit('k2cf: pop failed', 0)
c     
c     Get exp(K)
c     
      call matrix_exp( basis, g_temp )
c     
c     Transform MO's
c     
      call ga_dgemm('n', 'n', nbf, nbf, nbf, one, g_vecs,
     $     g_temp, zero, g_result)
c
      if (.not.ga_destroy(g_temp))
     $     call errquit('rohf_k2cf: cannot destroy temp',0)
c
      call ga_orthog_mos(basis, g_result)
c
      end


