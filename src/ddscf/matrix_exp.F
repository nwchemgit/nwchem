      subroutine matrix_exp(basis, g_k)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "cscfps.fh"
c     
c     Replace the input matrix with a unitary approximation
c     to its exponential.  Current algorithm is
c     
c     a) Determine n (a power of two) such that exp(K/n) is
c     well approximated by a taylor series
c     
c     b) Use third order taylor series on K/n
c     
c     c) Make the approximation to exp(K/n) unitary with
c     iterative orthogonalization
c     
c     d) Compute (exp(K/n))**n
c     
c     e) Reapply iterative orthogonalization
c     
      integer g_k, basis
c     
      integer nbf, nmo, type, g_tmp1, g_tmp2, i
      double precision knorm, scale, third, one, half, tmp, zero
      parameter (third = 1.0d0/3.0d0, one = 1.0d0, half = 0.5d0,
     $     zero = 0.0d0)
      integer power
c
      if (oscfps) call pstat_on(ps_matexp)
c
      call ga_inquire(g_k, type, nbf, nmo)
      if (nbf.ne.nmo) call errquit('matrix_exp: must be square', 0)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'mat_exp:tmp1',
     $     1, 1, g_tmp1)) call errquit('mat_exp:ga_create?',nbf)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'mat_exp:tmp2',
     $     1, 1, g_tmp2)) call errquit('mat_exp:ga_create?',nbf)
c     
c     Determine how to scale the matrix
c     
      knorm = dsqrt(ga_ddot(g_k, g_k)/(nbf*nbf))
      power = 1
 10   if (knorm .gt. 0.01d0) then
         power = power * 2
         knorm = knorm * 0.5d0
         goto 10
      endif
*      if (ga_nodeid() .eq. 0) write(6,*) ' power ', power
      scale = 1.0d0/power
      call ga_dscal(g_k, scale)
c     
c     Compute taylor series approximation thru 3rd order
c     
c     1 + K/3 -> g_tmp1
c     
      call ga_copy(g_k, g_tmp1)
      call ga_dscal(g_tmp1, third)
      do i = ga_nodeid()+1, nbf, ga_nnodes()
         call ga_get(g_tmp1, i, i, i, i, tmp, 1)
         tmp = tmp + one
         call ga_put(g_tmp1, i, i, i, i, tmp, 1)
      enddo
c     
c     1 + K/2 * g_tmp1 -> g_tmp2
c     
      call ga_dgemm('n', 'n', nbf, nbf, nbf, half, g_k, g_tmp1,
     $     zero, g_tmp2)
      do i = ga_nodeid()+1, nbf, ga_nnodes()
         call ga_get(g_tmp2, i, i, i, i, tmp, 1)
         tmp = tmp + one
         call ga_put(g_tmp2, i, i, i, i, tmp, 1)
      enddo
c     
c     1 + K*g_tmp2 -> g_k
c     
      call ga_dgemm('n', 'n', nbf, nbf, nbf, one, g_k, g_tmp2,
     $     zero, g_tmp1)
      do i = ga_nodeid()+1, nbf, ga_nnodes()
         call ga_get(g_tmp1, i, i, i, i, tmp, 1)
         tmp = tmp + one
         call ga_put(g_tmp1, i, i, i, i, tmp, 1)
      enddo
      call ga_copy(g_tmp1, g_k)
c     
      if (.not. (ga_destroy(g_tmp1) .and. ga_destroy(g_tmp2)))
     $     call errquit('mat_exp: ga_destroy failed', 0)
c     
c     Now have in g_k an approximation to exp(K/scale)
c     
c     Make this unitary and then scale as necessary
c     
      call ga_iter_orthog(basis, g_k, .false.)
c     
      if (power .ne. 1) then
         if (.not. ga_create(MT_DBL, nbf, nbf, 'mat_exp:tmp1',
     $        1, 1, g_tmp1)) call errquit('mat_exp:ga_create?',nbf)
 20      call ga_dgemm('n', 'n', nbf, nbf, nbf, one, g_k, g_k,
     $        zero, g_tmp1)
         call ga_copy(g_tmp1, g_k)
         power = power / 2
         if (power .ne. 1) goto 20
         if (.not. ga_destroy(g_tmp1))
     $        call errquit('mat_exp: ga_destroy failed', 0)
c     
c     Make sure the final result is also unitary
c     
         call ga_iter_orthog(basis, g_k, .false.)
      endif
c     
      if (oscfps) call pstat_off(ps_matexp)
c
      end
