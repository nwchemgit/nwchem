C$Id$
      subroutine ga_init_project(n, nroot, maxiter, maxsub, tol,
     $     precond, product, oprint, eval0)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cproject.fh"
      integer n                 ! Matrix size
      integer nroot             ! No. of roots to project
      integer maxiter
      integer maxsub
      double precision tol
      external product, precond ! Functions for solvers
      logical oprint
      double precision eval0
c     
      double precision acc_proj(max_proj)
      integer iter, iroot
      logical ga_iter_diag
      external ga_iter_diag
      external ga_block_project ! For T3D linker
c     
c     Comment out the next line to enable roots to be improved
c     by calling this routine multiple times
c
***** if (oproject) return
      if (.not. oproject) then
         if (nroot .gt. max_proj) call errquit
     $        ('ga_init_project: max_proj exceeded', nroot, INPUT_ERR)
c     
         nroot_proj = nroot
         n_proj     = n
*ga:1:0
         if (.not. ga_create(MT_DBL, n, nroot, 'ga_project', 0, nroot,
     $        g_proj)) call errquit('ga_init_project: ga_create?', 0,
     &       GA_ERR)
         call ga_zero(g_proj)
      else
         if (nroot_proj.ne.nroot .or. n_proj.ne.n) call errquit
     $        ('ga_init_project: mismatch of dimensions', 0, GA_ERR)
      endif
c     
      if (.not. ga_iter_diag(n, nroot, maxiter, maxsub, tol,
     $     precond, product, oprint, eval0, g_proj, e_proj,
     $     acc_proj, iter)) then
         do iroot = nroot,1,-1
            if (acc_proj(iroot).gt.1.0d0) nroot_proj = nroot_proj - 1
         enddo
         if (nroot_proj .ne. nroot) then
            if (ga_nodeid() .eq. 0) then
               write(6,1) nroot_proj, nroot
 1             format('  ga_init_project: only got',i3,' roots of',i3)
               call util_flush(6)
            endif
            if (nroot_proj .eq. 0) call errquit
     $           ('ga_init_project: no roots extracted', 0, GA_ERR)
         endif
      endif
c     
      oproject = .true.
c     
      end
      subroutine ga_term_project
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "cproject.fh"
c
      if (oproject) then
         if(.not. ga_destroy(g_proj))
     &      call errquit('ga_term_project: ga_destroy failed',911,
     &       GA_ERR)
      endif
      oproject = .false.
c
      end
      subroutine ga_precond_project(g_x, shift, precond, g_work)
      implicit none
#include "global.fh"
#include "cproject.fh"
      integer g_x
      double precision shift
      external precond
      integer g_work
c     
c     The user provides the routine precond that computes an
c     approximation to the shifted inverse of the matrix
c     being solved.  This routine combines this preconditioner
c     with the eigenvectors computed by ga_init_project to
c     come up with an improved preconditioner
c     
c     Let A be the matrix being inverted, e(u) be the u'th 
c     eigenvector with eigenvalue v(u), D be the approximation
c     to A inverted by precond() and P be a projection operator
c     onto the subspace of known eigenvectors.
c     
c     Then, this routine uses as an approximation to A
c     
c     PDP + sum(u=1,nroot) v(u) e(u) e(u)T
c     
c     where
c     
c     P = 1 - sum(u=1,nroot) e(u) e(u)T
c     
c     The inverse of this matrix is
c     
c     PD^-1P + sum(u=1,nroot) v(u)^-1 e(u) e(u)T
c     
      integer iroot, n, nroot
      double precision scale
c     
      if (.not. oproject) then
         call precond(g_x, shift)
      else
         n = n_proj
         nroot = nroot_proj
c     
c     Form P D^-1 P b in work
c     
         call ga_copy_patch('n', g_x, 1, n, 1, 1, g_work, 1, n, 1, 1)
c     
         do iroot = 1, nroot
            scale = ga_ddot_patch(
     $           g_x, 'n', 1, n, 1, 1,
     $           g_proj, 'n', 1, n, iroot, iroot)
            call ga_dadd_patch(
     $           1.0d0, g_work, 1, n, 1, 1,
     $           -scale, g_proj, 1, n, iroot, iroot,
     $           g_work, 1, n, 1, 1)
         enddo
c     
         call precond(g_work, shift)
c     
         do iroot = 1, nroot
            scale = ga_ddot_patch(
     $           g_work, 'n', 1, n, 1, 1,
     $           g_proj, 'n', 1, n, iroot, iroot)
            call ga_dadd_patch(
     $           1.0d0, g_work, 1, n, 1, 1,
     $           -scale, g_proj, 1, n, iroot, iroot,
     $           g_work, 1, n, 1, 1)
         enddo
c     
c     Add on the piece due to the eigenvectors
c     
         do iroot = 1, nroot
            scale = ga_ddot_patch(
     $           g_x, 'n', 1, n, 1, 1,
     $           g_proj, 'n', 1, n, iroot, iroot)
            scale = scale / e_proj(iroot)
            call ga_dadd_patch(
     $           1.0d0, g_work, 1, n, 1, 1,
     $           scale, g_proj, 1, n, iroot, iroot,
     $           g_work, 1, n, 1, 1)
         enddo
c     
c     Copy the result back
c     
         call ga_copy_patch('n', g_work, 1, n, 1, 1, g_x, 1, n, 1, 1)
      endif
c     
      end
