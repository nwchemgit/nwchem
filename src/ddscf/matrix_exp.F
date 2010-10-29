      subroutine matrix_exp(g_k)
C$Id$
      implicit none
#include "errquit.fh"
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
c     c) If necessary make exp(-K/n) unitary
c
c     d) Compute (exp(K/n))**n
c     
c     e) If necessary make exp(-K) unitary
c
      integer g_k               ! [input] GA handle to matrix
c     
      integer nbf, nmo, type, g_tmp1, g_tmp2, i
      double precision knorm, scale, third, one, half, zero
*     integer j
*     double precision xx
*      logical obad
      parameter (third = 1.0d0/3.0d0, one = 1.0d0, half = 0.5d0,
     $     zero = 0.0d0)
      integer power
      integer dummy_basis
      integer ilo, ihi, jlo, jhi ! Used to speed extraction of diagonals
c
      dummy_basis = -1          ! Not used but passed to ga_iter_orthog
c
      if (oscfps) call pstat_on(ps_matexp)
c
      call ga_inquire(g_k, type, nbf, nmo)
      if (nbf.ne.nmo) call errquit('matrix_exp: must be square', 0,
     &       UNKNOWN_ERR)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'mat_exp:tmp1',
     $     32, 32, g_tmp1)) call errquit('mat_exp:ga_create?',nbf,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'mat_exp:tmp2',
     $     32, 32, g_tmp2)) call errquit('mat_exp:ga_create?',nbf,
     &       GA_ERR)
c
c     Determine range of local diagonal values ... assume same for tmp1/2
c
      call ga_distribution(g_tmp1, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         ilo = max(ilo,jlo)
         ihi = min(ihi,jhi)
      endif
c     
c     Determine how to scale the matrix
c     
      call ga_maxelt(g_k, knorm)
      power = 1
 10   if (knorm .gt. 0.01d0) then
         power = power * 2
         knorm = knorm * 0.5d0
         goto 10
      endif
c
      knorm = knorm*power       ! recover original value
*      if (ga_nodeid() .eq. 0) write(6,*) ' power, knorm ',power,knorm
      scale = 1.0d0/power
      call ga_dscal(g_k, scale)
c     
c     Compute taylor series approximation thru 3rd order
c     
c     1 + K/3 -> g_tmp1
c     
      call ga_copy(g_k, g_tmp1)
      call ga_dscal(g_tmp1, third)
      do i = ilo, ihi
         call ga_acc(g_tmp1, i, i, i, i, 1.0d0, 1, 1.0d0)
      enddo
c     
c     1 + K/2 * g_tmp1 -> g_tmp2
c     
      call ga_dgemm('n', 'n', nbf, nbf, nbf, half, g_k, g_tmp1,
     $     zero, g_tmp2)
      do i = ilo, ihi
         call ga_acc(g_tmp2, i, i, i, i, 1.0d0, 1, 1.0d0)
      enddo
c     
c     1 + K*g_tmp2 -> g_k
c     
      call ga_dgemm('n', 'n', nbf, nbf, nbf, one, g_k, g_tmp2,
     $     zero, g_tmp1)
      do i = ilo, ihi
         call ga_acc(g_tmp1, i, i, i, i, 1.0d0, 1, 1.0d0)
      enddo
      if (.not. ga_destroy(g_tmp2))
     $     call errquit('mat_exp: ga_destroy failed', 0, GA_ERR)
c
      call ga_copy(g_tmp1, g_k)
c     
c     Now have in g_k an approximation to exp(K/scale)
c     
c     Make this unitary and then scale if necessary
c     
      if (power.gt.1 .or. knorm.gt.1d-3) then
*         if (ga_nodeid().eq.0)write(6,*) ' matexp doing first orthog '
         call ga_iter_orthog(dummy_basis, g_k, .false., .true.)
      endif
c     
      if (power .ne. 1) then
 20      call ga_dgemm('n', 'n', nbf, nbf, nbf, one, g_k, g_k,
     $        zero, g_tmp1)
         call ga_copy(g_tmp1, g_k)
         power = power / 2
         if (power .ne. 1) goto 20
c     
c     Make sure the final result is also unitary
c     
         if (knorm.gt.1) then
*            if (ga_nodeid().eq.0)write(6,*) ' matexp doing 2nd orthog'
            call ga_iter_orthog(dummy_basis, g_k, .false., .true.)
         endif
      endif
c$$$c
c$$$c     test for unitary property ... paranoia at work
c$$$c
c$$$      call ga_dgemm('n', 't', nbf, nbf, nbf,
c$$$     $     one, g_k, g_k, zero, g_tmp1)
c$$$      call ga_distribution(g_tmp1, ga_nodeid(), ilo, ihi, jlo, jhi)
c$$$      obad = .false.
c$$$      do i = ilo, ihi
c$$$         do j = jlo, jhi
c$$$            call ga_get(g_tmp1, i, i, j, j, xx, 1)
c$$$            if (i .eq. j) xx = xx - one
c$$$            if (abs(xx).gt.1d-12) then 
c$$$               write(6,*) i, j, xx
c$$$               call util_flush(6)
c$$$               obad = .true.
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$      call ga_sync()
c$$$      if (obad) call errquit('matexp: bad', 0)
c     
      if (.not. ga_destroy(g_tmp1))
     $     call errquit('mat_exp: ga_destroy failed', 0, GA_ERR)
c     
      if (oscfps) call pstat_off(ps_matexp)
c
      end
