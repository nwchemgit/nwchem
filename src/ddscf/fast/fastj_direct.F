      subroutine fastj_direct_sum_init(basis, g_dens, eps, data)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "fudge.fh"
#include "pottimes.fh"
#include "util.fh"
      integer basis, g_dens, data(2)
      double precision eps
c
c     Initialize the direct sum of the Coulomb potential for
c     the DFT by replicating the density (eventually lower
c     triangular and prescreening for significant pairs
c     of shells)
c     
c$$$      integer g_xlm
c$$$      double precision center(3), q(-2:2,0:2), dd, mp(-2:2,0:2)
c$$$      integer i, j, l, m
c
      integer k_dens, l_dens, nbf
      logical ofirstcall
      data ofirstcall/.true./

      if (.not. bas_numbf(basis, nbf)) call errquit('fastj:ds:nbf',0,
     &       BASIS_ERR)

      otime = util_print('pottime',print_debug)
      time_init = 0.0d0
      time_make_poles = 0.0d0
      time_poles = 0.0d0
      time_direct = 0.0d0

      gg_dens = g_dens

c
c$$$      if (.not. ga_create(mt_dbl, 9*nbf, nbf, 'xlm', 0, 0, 
c$$$     $     g_xlm)) call errquit(' ga_xlm ? ',0)
c$$$c
c$$$      center(1) = 0.0d0
c$$$      center(2) = 0.0d0
c$$$      center(3) = 0.0d0
c$$$      call dfill(9, 0.0d0, mp, 1)
c$$$      call my_xlm_make_ao_poles(basis, center, 2, g_xlm)
c$$$      do i = 1, 3
c$$$         do j = 1, 3
c$$$            call ga_get(g_dens, i, i, j, j, dd, 1)
c$$$            call my_xlm_get_ij(g_xlm, i, j, 2, q)
c$$$            do l = 0, 2
c$$$               do m = -l, l
c$$$                  mp(m,l) = mp(m,l) + q(m,l)*dd
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      write(6,*) ' Poles from first shell ?'
c$$$      call xlm_print(2,mp,2)
c$$$c
c$$$      call ga_print(g_dens)
c$$$      call ga_print(g_xlm)
c$$$      if (.not. ga_destroy(g_xlm)) call errquit('fred',0)
c
      if (ofirstcall) then
         call xlm_init()
         call anl_fit_init()
         call xlm_coeff_inv(lmax2,d,dinv)
         ofirstcall = .false.
      endif
c
c
      if (.not. ma_push_get(mt_dbl, nbf**2, 'fjdsden', l_dens, k_dens))
     $     call errquit('fastj: ds: ma ', nbf*nbf, MA_ERR)
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
#include "errquit.fh"
#include "mafdecls.fh"
#include "fudge.fh"
#include "pottimes.fh"
      integer data(2)
c
      if (.not. ma_pop_stack(data(1))) call errquit('fastj:dstidy?',0,
     &       MA_ERR)
      write(6,*) ' nq, ncall, nq/ncall ', nq, ncall, nq/ncall
      if (otime) then
         write(6,1) time_init, time_make_poles, time_poles, time_direct
 1       format(2x,4f10.4)
      endif
c
      end
      subroutine fastj_direct_sum(basis, dens, npts, r, potential, tol)
      implicit none
#include "errquit.fh"
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
      call potential_list(basis, gg_dens, npts, r, potential, tol)
c
      nq = nq + npts
      ncall = ncall + 1
      return
c
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('fastj_ds: bas_numcont failed for basis', basis, BASIS_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('fastj_ds: bas_numbf failed for basis', basis, BASIS_ERR)
c     
      call dfill(npts, 0.0d0, potential, 1)
      tols = tol / schwarz_max()

      do jshell = 1, nshell
         if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     $        call errquit('fastj_ds: bas_cn2bfr ?', basis, BASIS_ERR)
         do ishell = 1, jshell
            if (schwarz_shell(jshell,ishell) .gt. tols) then
               if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) 
     $           call errquit('fastj_ds: bas_cn2bfr ?', basis,
     &       BASIS_ERR)
               fact = 1.0d0
               if (ishell .ne. jshell) fact = 2.0d0
               call potential_from_shell_pair(basis, ishell, jshell,
     $              dens(ilo + (jlo-1)*nbf), nbf,
     $              tol, d, dinv, lmax2, 
     $              npts, r, potential, fact, .false., .false., 1)
            endif
         enddo
      enddo
**            stop 99
c
      end




