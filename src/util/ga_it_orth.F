      subroutine ga_iter_orthog(basis, g_c, ometric, osym)
C$Id$
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer basis
      integer g_c
      logical ometric
      logical osym
c     
c     Iteratively orthogonalize the vectors in g_c using quadratically
c     convergent algorithms
c
c     if (ometric) then
c       use the overlap integrals from basis as a metric
c     else
c       use a unit matrix as a metric
c     endif
c
c     if (osym) then
c       symmetrically orthogonalize the vectors using
c
c       C <-- C*(1.5*I - 0.5*S)
c
c     else
c       orthogonalize higher index vectors to lower
c       index vectors (e.g., so as to avoid mixing
c       virtual orbitals into occupied orbitals)
c
c       C <-- C*(I + U(I-S)) 
c
c       where U(M) = upper triangle with halved diagonal
c     endif
c     
c     Two n*n scratch GAs are used.
c     
c     
      integer g_s, g_tmp, g_ms
      integer iter, maxiter, nbf, nmo, type, i
      integer ilo, ihi, jlo, jhi
      double precision tol, zero, one, mone, mhalf, sii, err
      parameter (tol = 1.0d-15, maxiter = 9)
      parameter (zero=0.0d0, one=1.0d0, mhalf=-0.5d0, mone=-1.0d0)
c     
c     Make scratch arrays
c     
      call ga_inquire(g_c, type, nbf, nmo)
      if (nmo .ne. nbf) call errquit('iter_orthog:nmo.ne.nbf',nbf,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'iter_orthog:over',
     $     32, 32, g_s)) call errquit('iter_orthog:ga_create?',nbf,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, nbf, nmo, 'iter_orthog:tmp',
     $     32, 32, g_tmp)) call errquit('iter_orthog:ga_create?',nbf,
     &       GA_ERR)
      if (nmo .ne. nbf) then
         if (.not. ga_create(MT_DBL, nmo, nmo, 'iter_orthog:ms',
     $        32, 32, g_ms)) call errquit('iter_orthog:ga_create?',nbf,
     &       GA_ERR)
      else
         g_ms = g_s             ! Save space usually
      endif
c
c     Determine range of local diagonals for shifting g_ms
c
      call ga_distribution(g_ms, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         ilo = max(ilo,jlo)
         ihi = min(ihi,jhi)
      endif
c     
      do iter = 1, maxiter
c     
c     Make the overlap matrix (S) in the MO basis in g_s
c     
         if (ometric) then
            call ga_zero(g_s)
            call int_1e_ga(basis, basis, g_s, 'overlap', .false.)
            call ga_dgemm('n','n',nbf,nmo,nbf,one,g_s,g_c,
     $           zero,g_tmp)
            call ga_dgemm('t','n',nmo,nmo,nbf,one,g_c,g_tmp,
     $           zero,g_ms)
         else
            call ga_dgemm('t','n',nmo,nmo,nbf,one,g_c,g_c,
     $           zero,g_ms)
         endif
c     
c     Form S - I and determine the error norm
c     
         do i = ilo, ihi
            call ga_acc(g_ms, i, i, i, i, mone, 1, one)
         enddo
         call ga_maxelt(g_ms, err)
*         if ((ga_nodeid().eq.0))
*     $        write(6,*) ' iter_orthog: err = ', err
c     
c     Form and apply the correction using locality
c     
         if (osym) then
            call ga_dscal(g_ms, mhalf)
            do i = ilo, ihi
               call ga_acc(g_ms, i, i, i, i, one, 1, one)
            enddo
         else
            call ga_dscal(g_ms, mone)
            do i = ilo, ihi
               call ga_get(g_ms, i, i, i, i, sii, 1)
               sii = 0.5d0*sii + one
               call ga_put(g_ms, i, i, i, i, sii, 1)
            enddo
         endif
c            
         call ga_dgemm('n','n',nbf,nmo,nmo,one,g_c,g_ms,zero,g_tmp)
         call ga_copy(g_tmp, g_c)
c     
c     Finally check for convergence (noting quadratic convergence)
c     
         if (err.lt.sqrt(tol)) goto 1000
c     
      enddo
      write(6,*) ' iter_orthog: failed to converge, error = ', err
      call errquit('iter_orthog: failed to converge', 0, CALC_ERR)
c     
 1000 continue
      if (.not. (ga_destroy(g_s) .and. ga_destroy(g_tmp)))
     $     call errquit('iter_orthog: ga_destroy failed', 0, GA_ERR)
      if (nbf .ne. nmo) then
         if (.not. ga_destroy(g_ms))
     $     call errquit('iter_orthog: ga_destroy failed', 0, GA_ERR)
      endif
c
      call ga_sync()
c
      end
