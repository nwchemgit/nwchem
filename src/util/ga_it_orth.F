      subroutine ga_iter_orthog(basis, g_c, ometric, osym)
C$Id: ga_it_orth.F,v 1.2 1996-01-19 07:32:48 d3g681 Exp $
      implicit none
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer basis
      integer g_c
      logical ometric
      logical osym
c     
c     Iteratively orthogonalize the vectors in g_c
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
      integer l_buf, k_buf      ! MA handles for workspace
      integer iter, maxiter, nbf, nmo, type, i
      double precision tol, zero, one, mone, mhalf, sii, err
      parameter (tol = 1.0d-15, maxiter = 9)
      parameter (zero=0.0d0, one=1.0d0, mhalf=-0.5d0, mone=-1.0d0)
c     
c     Make scratch arrays
c     
      call ga_inquire(g_c, type, nbf, nmo)
      if (nmo .ne. nbf) call errquit('iter_orthog:nmo.ne.nbf',nbf)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'iter_orthog:over',
     $     32, 32, g_s)) call errquit('iter_orthog:ga_create?',nbf)
      if (.not. ga_create(MT_DBL, nbf, nmo, 'iter_orthog:tmp',
     $     32, 32, g_tmp)) call errquit('iter_orthog:ga_create?',nbf)
      if (nmo .ne. nbf) then
         if (.not. ga_create(MT_DBL, nmo, nmo, 'iter_orthog:ms',
     $        32, 32, g_ms)) call errquit('iter_orthog:ga_create?',nbf)
      else
         g_ms = g_s             ! Save space usually
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
         call ga_sync()
         do i = ga_nodeid()+1, nmo, ga_nnodes()
            call ga_get(g_ms, i, i, i, i, sii, 1)
            sii = sii - one
            call ga_put(g_ms, i, i, i, i, sii, 1)
         enddo
         call ga_sync()
	 call ga_maxelt(g_ms, err)
*         if ((ga_nodeid().eq.0))
*     $        write(6,*) ' iter_orthog: err = ', err
c     
c     Form and apply the correction ... SHOULD USE LOCALITY
c     
         if (osym) then
            call ga_dscal(g_ms, mhalf)
            do i = ga_nodeid()+1, nmo, ga_nnodes()
               call ga_get(g_ms, i, i, i, i, sii, 1)
               sii = sii + one
               call ga_put(g_ms, i, i, i, i, sii, 1)
            enddo
            call ga_sync()
         else
            call ga_zero(g_tmp)
            call ga_dscal(g_ms, mone)
            if (.not. ma_push_get(MT_DBL, nbf, 'ga_iter_orthog',
     $           l_buf, k_buf)) call errquit
     $           ('ga_iter_orthog: insufficient memory', nbf)
            do i = ga_nodeid()+1, nmo, ga_nnodes()
               call ga_get(g_ms, 1, i, i, i, dbl_mb(k_buf), 1)
               sii = dbl_mb(k_buf+i-1)
               sii = 0.5d0*sii + one
               dbl_mb(k_buf+i-1) = sii
               call ga_put(g_tmp, 1, i, i, i, dbl_mb(k_buf), 1)
            enddo
            if (.not. ma_pop_stack(l_buf)) call errquit
     $           ('ga_iter_orthog: pop stack?', nbf)
            call ga_sync()
            call ga_copy(g_tmp, g_ms)
         endif
c            
         call ga_zero(g_tmp)
         call ga_dgemm('n','n',nbf,nmo,nmo,one,g_c,g_ms,zero,g_tmp)
         call ga_copy(g_tmp, g_c)
c     
c     Finally check for convergence (noting quadratic convergence)
c     
         if (err.lt.sqrt(tol)) goto 1000
c     
      enddo
      write(6,*) ' iter_orthog: failed to converge, error = ', err
      call errquit('iter_orthog: failed to converge', 0)
c     
 1000 continue
      if (.not. (ga_destroy(g_s) .and. ga_destroy(g_tmp)))
     $     call errquit('iter_orthog: ga_destroy failed', 0)
      if (nbf .ne. nmo) then
         if (.not. ga_destroy(g_ms))
     $     call errquit('iter_orthog: ga_destroy failed', 0)
      endif
c
      end
