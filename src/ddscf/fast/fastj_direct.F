      subroutine fastj_direct_sum_init(basis, g_dens, eps, data)
      implicit none
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "fudge.fh"
      integer basis, g_dens, data(2)
      double precision eps
c
c     Initializ the direct sum of the Coulomb potential for
c     the DFT by replicating the density (eventually lower
c     triangular and prescreening for significant pairs
c     of shells)
c     
      integer k_dens, l_dens, nbf
      logical ofirstcall
      data ofirstcall/.true./
c
      if (ofirstcall) then
         call xlm_init()
         call anl_fit_init()
         call xlm_coeff_inv(lmax2,d,dinv)
         ofirstcall = .false.
      endif
c
      if (.not. bas_numbf(basis, nbf)) call errquit('fastj:ds:nbf',0)
c
      if (.not. ma_push_get(mt_dbl, nbf**2, 'fjdsden', l_dens, k_dens))
     $     call errquit('fastj: ds: ma ', nbf*nbf, 0)
c
      if (ga_nodeid() .eq. 0)
     $     call ga_get(g_dens, 1, nbf, 1, nbf, dbl_mb(k_dens), nbf)
      call ga_brdcst(11, dbl_mb(k_dens), 8*nbf**2, 0)
c
      data(1) = l_dens
      data(2) = nint(log10(eps))
c
      nq = 0
      ncall = 0
c
      end
      subroutine fastj_direct_sum_tidy(data)
      implicit none
#include "mafdecls.fh"
#include "fudge.fh"
      integer data(2)
c
      if (.not. ma_pop_stack(data(1))) call errquit('fastj:dstidy?',0)
      write(6,*) ' nq, ncall, nq/ncall ', nq, ncall, nq/ncall
c
      end
      subroutine fastj_direct_sum(basis, dens, npts, r, potential, tol)
      implicit none
#include "bas.fh"
#include "util.fh"
#include "fudge.fh"
#include "schwarz.fh"
c
      integer basis, npts
      double precision r(3,npts), potential(npts), dens(*)
      double precision tol
c
c     Compute the electrostatic potential due to the density
c     at the given list of points.
c
c     Done by direct summation of the potential due to the
c     primitive charge distributions with screening on the primitives.
c
      integer nshell
      integer ishell, jshell
      integer ilo, ihi, jlo, jhi, nbf
      double precision fact, tols
c
      nq = nq + npts
      ncall = ncall + 1
c
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('fastj_ds: bas_numcont failed for basis', basis)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('fastj_ds: bas_numbf failed for basis', basis)
c     
      call dfill(npts, 0.0d0, potential, 1)
      tols = tol / schwarz_max()
c     
      do jshell = 1, nshell
         if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     $        call errquit('fastj_ds: bas_cn2bfr ?', basis)
         do ishell = 1, jshell
            if (schwarz_shell(jshell,ishell) .gt. tols) then
               if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) 
     $           call errquit('fastj_ds: bas_cn2bfr ?', basis)
               fact = 1.0d0
               if (ishell .ne. jshell) fact = 2.0d0
               call potential_from_shell_pair(basis, ishell, jshell,
     $              dens(ilo + (jlo-1)*nbf), nbf,
     $              tol, d, dinv, lmax2, 
     $              npts, r, potential, fact)
            endif
         enddo
      enddo
c
      end

