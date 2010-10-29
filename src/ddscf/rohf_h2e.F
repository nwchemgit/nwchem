      subroutine rohf_hessv_2e( basis, geom, nbf, nmo, nclosed, nopen,
     $     g_movec, oskel, g_x, acc, g_ax )
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "cscfps.fh"
#include "bgj.fh"
#include "rtdb.fh"
c     
c     Return the ROHF orbital 2e-Hessian vector product, g_ax = A * g_x
c     
      integer basis, geom       ! basis & geom handle
      integer nbf, nclosed, nopen ! Basis size and occupation
      integer nmo               ! No. of linearly dependent MOs
      integer g_movec           ! MO coefficients
      logical oskel
      integer g_x               ! Argument
      double precision acc      ! Accuracy of "Fock" construction
      integer g_ax              ! Hessian product
c     
      integer nvir, voff, xoff, xend
      integer ooff, oend
c     
      integer nfock, ifock
      integer g_tmp1, g_tmp2, g_tmp3
c
      integer ivec, nvec, gtype, vlen
      integer nset, max_nvec
      parameter(nset=3, max_nvec=300)
      integer iv_dens(nset*max_nvec),iv_fock(nset*max_nvec)
      double precision jfac(nset*max_nvec),kfac(nset*max_nvec)
c
      integer g_dcv,  g_dcp,  g_dpv
      double precision tol2e
      logical odebug
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision one, zero, mone, four, half, mhalf, two
      double precision itol_floor, itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0)
c     
c     This get's cleaned up a lot when the Fock build accepts
c     a single multiple dimension GA for its input.
c
c     Lots more to do on memory management.
c
      call ga_inquire(g_x,gtype,vlen,nvec)
      if (nvec .gt. max_nvec) call errquit
     $     ('rohf_h2e: too many vectors',nvec, UNKNOWN_ERR)
      odebug = util_print('rohf_hessv', print_debug)
c     
      tol2e = min(max(acc,itol_floor),itol_ceil)
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
      ooff = nclosed + 1
      oend = nclosed + nopen
c     
c     Allocate global for work space
c     
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp2')
      if (nopen .gt. 0) then
         g_tmp3 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp3')
      endif
c
c     Allocate arrays and generate the AO densities
c     
      nfock = 0
      do ivec = 1, nvec
         g_dcv = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dcv')
         nfock          = nfock + 1
         iv_dens(nfock) = g_dcv
         jfac(nfock)    =  1.0d0   ! CV
         kfac(nfock)    = -0.5d0   ! CV
         if (nopen .gt. 0) then
            g_dcp = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dcp')
            g_dpv = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dpv')
            nfock          = nfock + 1
            iv_dens(nfock) = g_dcp
            jfac(nfock)    =  1.0d0 ! CP
            kfac(nfock)    = -0.5d0 ! CP
            nfock          = nfock + 1
            iv_dens(nfock) = g_dpv
            jfac(nfock) =  0.0d0 ! PV
            kfac(nfock) =  0.5d0 ! PV
         endif
c     
c     Compute CV, PV & CP "densities" from argument vector
c     
         call ga_zero(g_tmp1)
         call ga_zero(g_tmp2)
         call ga_zero(g_dcv)
         call ga_vec_to_mat(g_tmp1, 1, nvir, 1, nclosed,
     $                      g_x, 1, ivec)
*         call ga_copy_patch('n', g_x, 1, (nclosed*nvir), ivec, ivec,
*     $        g_tmp1, 1, nvir, 1, nclosed )
         call ga_matmul_patch('n','t',one,zero,
     $        g_tmp1, 1, nvir, 1, nclosed,
     $        g_movec, 1, nclosed, 1, nbf,
     $        g_tmp2, 1, nvir, 1, nbf )
         call ga_matmul_patch('n','n',one,zero,
     $        g_movec, 1, nbf, voff, nmo,
     $        g_tmp2, 1, nvir, 1, nbf,
     $        g_dcv, 1, nbf, 1, nbf)
c     
         if (nopen .gt. 0) then
c     
c     PV density
c     
            xoff = nclosed*nvir + 1
            xend = nclosed*nvir + nopen*nvir
            call ga_zero(g_tmp1)
            call ga_zero(g_tmp2)
            call ga_zero(g_dpv)
            call ga_copy_patch('n',
     $           g_x, xoff, xend, ivec, ivec,
     $           g_tmp1, 1, nvir, 1, nopen )
            call ga_matmul_patch('n', 't', one, zero,
     $           g_tmp1, 1, nvir, 1, nopen,
     $           g_movec, ooff, oend, 1, nbf,
     $           g_tmp2, 1, nvir, 1, nbf )
            call ga_matmul_patch('n', 'n', one, zero,
     $           g_movec, 1, nbf, voff, nmo,
     $           g_tmp2, 1, nvir, 1, nbf,
     $           g_dpv, 1, nbf, 1, nbf )
c     
c     CP density
c     
            call ga_zero(g_tmp1)
            call ga_zero(g_tmp2)
            call ga_zero(g_dcp)
            xoff = (nclosed+nopen)*nvir + 1
            xend = (nclosed+nopen)*nvir + nclosed*nopen
            call ga_copy_patch('n', g_x, xoff, xend, ivec, ivec,
     $           g_tmp1, 1, nopen, 1, nclosed )
            call ga_matmul_patch('n', 't', one, zero,
     $           g_tmp1, 1, nopen, 1, nclosed,
     $           g_movec, 1, nclosed, 1, nbf,
     $           g_tmp2, 1, nopen, 1, nbf )
            call ga_matmul_patch('n', 'n', one, zero,
     $           g_movec, 1, nbf, ooff, oend,
     $           g_tmp2, 1, nopen, 1, nbf,
     $           g_dcp, 1, nbf, 1, nbf )
         endif
c     
c     We actually use dcp+dpv and dcp-dpv to build fock matrix
c     
         if (nopen .gt. 0) then
            call ga_dadd(one, g_dcp, one, g_dpv, g_tmp1)
            call ga_dadd(one, g_dcp,mone, g_dpv, g_tmp2)
            call ga_copy(g_tmp1, g_dcp)
            call ga_copy(g_tmp2, g_dpv)
         endif
c
         call ga_symmetrize(g_dcv)
         if (nopen .gt. 0) then
            call ga_symmetrize(g_dpv)
            call ga_symmetrize(g_dcp)
         endif
      enddo
c
      if (.not. ga_destroy(g_tmp1)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_tmp3)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      endif
c
c     Fock build
c     
      do ifock = 1, nfock
         iv_fock(ifock) = ga_create_atom_blocked
     $        (geom, basis, 'rohf:h2e fock')
         call ga_zero(iv_fock(ifock))
      enddo
c
      call fock_2e( geom, basis, nfock, jfac, kfac,
     $     tol2e, oskel, iv_dens, iv_fock,.false.)
c     
      if (oskel) then
         if (oscfps) call pstat_on(ps_sym_sym)
         do ifock = 1, nfock
            call sym_symmetrize(geom, basis, .false., iv_fock(ifock))
         enddo
         if (oscfps) call pstat_off(ps_sym_sym)
      endif
c
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp2')
      if (nopen .gt. 0) then
         g_tmp3 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp3')
      endif
c
      ifock = 0
      do ivec = 1, nvec
         ifock = ifock + 1      ! History ... copy into tmp arrays
         call ga_copy(iv_fock(ifock),g_tmp1)
         g_dcv = iv_dens(ifock)
         if (nopen .gt. 0) then 
            ifock = ifock + 1
            call ga_copy(iv_fock(ifock),g_tmp2)
            g_dcp = iv_dens(ifock)
            ifock = ifock + 1
            call ga_copy(iv_fock(ifock),g_tmp3)
            g_dpv = iv_dens(ifock)
         endif
c     
c     P      =  4(ij|kl) - (ik|jl) - (il|kj)
c     ij,kl
c     
c     K      =  (ik|jl) + (il|kj)
c     ij,kl
c     
c     cv         cv          pv   cp
c     Z   =  2P.[D  ]  +  P.[D  + D  ]
c     
c     pv          cv           cp   pv
c     Z   =  0.5d0*Z   + 0.5*K.[D  - D  ]
c     
c     cp          cv           cp   pv
c     Z   =  0.5d0*Z   - 0.5*K.[D  - D  ]
c     
c     Add the Fock matrices together overwriting the density
c     matrices to form the results above
c     
         call ga_copy(g_tmp1, g_dcv)
         call ga_dscal(g_dcv, two)
         if (nopen .gt. 0) then
            call ga_dadd(one,  g_dcv,  one, g_tmp2, g_dcv)
            call ga_dadd(half, g_dcv, half, g_tmp3, g_dpv)
            call ga_dadd(half, g_dcv,mhalf, g_tmp3, g_dcp)
         endif
c     
c     Closed-Virtual bit
c     
         call ga_matmul_patch('n','n',two,zero,
     $        g_dcv, 1, nbf, 1, nbf,
     $        g_movec, 1, nbf, 1, nclosed,
     $        g_tmp1, 1, nbf, 1, nclosed )
         call ga_matmul_patch('t','n',one,zero,
     $        g_movec, voff, nmo, 1, nbf, 
     $        g_tmp1, 1, nbf, 1, nclosed,
     $        g_tmp2, 1, nvir, 1, nclosed )
         xoff = 1
         xend = nvir*nclosed
         call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed,
     $        g_ax, xoff, ivec, four, '+')
*         call ga_dadd_patch( four, g_tmp2, 1, nvir, 1, nclosed,
*     $        one, g_ax, xoff, xend, ivec, ivec,
*     $        g_ax, xoff, xend, ivec, ivec )
c     
         if (nopen .gt. 0) then
c     
c     Open-Virtual bit
c     
            call ga_matmul_patch('n','n',two,zero,
     $           g_dpv, 1, nbf, 1, nbf,
     $           g_movec, 1, nbf, ooff, oend,
     $           g_tmp1, 1, nbf, 1, nopen )
            call ga_matmul_patch('t','n',one,zero,
     $           g_movec, voff, nmo, 1, nbf,
     $           g_tmp1, 1, nbf, 1, nopen,
     $           g_tmp2, 1, nvir, 1, nopen )
            xoff = xend + 1
            xend = xend + nvir*nopen
            call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nopen,
     $          g_ax, xoff, ivec, four, '+')
*            call ga_dadd_patch( four, g_tmp2, 1, nvir, 1, nopen,
*     $           one, g_ax, xoff, xend, ivec, ivec,
*     $           g_ax, xoff, xend, ivec, ivec)
c     
c     Closed-Open bit
c     
            call ga_matmul_patch('n','n',two,zero,
     $           g_dcp, 1, nbf, 1, nbf,
     $           g_movec, 1, nbf, 1, nclosed,
     $           g_tmp1, 1, nbf, 1, nclosed )
            call ga_zero(g_tmp2)
            call ga_matmul_patch('t','n',one,zero,
     $           g_movec, ooff, oend, 1, nbf,
     $           g_tmp1, 1, nbf, 1, nclosed,
     $           g_tmp2, 1, nopen, 1, nclosed )
            xoff = xend + 1
            xend = xend + nopen*nclosed
            call ga_mat_to_vec(g_tmp2, 1, nopen, 1, nclosed,
     $          g_ax, xoff, ivec, four, '+')
*            call ga_dadd_patch( four, g_tmp2, 1, nopen, 1, nclosed,
*     $           one, g_ax, xoff, xend, ivec, ivec,
*     $           g_ax, xoff, xend, ivec, ivec)
         endif
      enddo
c
      do ifock = 1, nfock
         if (.not. ga_destroy(iv_dens(ifock)))
     $        call errquit('rohf:destroy?',0, GA_ERR)
         if (.not. ga_destroy(iv_fock(ifock)))
     $        call errquit('rohf:destroy?',0, GA_ERR)
      enddo
c
      if (.not. ga_destroy(g_tmp1)) call errquit('rohf:destroy?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit('rohf:destroy?',0,
     &       GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_tmp3)) call errquit('rohf:destroy?',0,
     &       GA_ERR)
      endif
c     
      end

