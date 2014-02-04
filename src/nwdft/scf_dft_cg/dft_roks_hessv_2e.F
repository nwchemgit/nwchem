      subroutine dft_roks_hessv_2e( basis, geom, nbf, nmo, nclosed,
     $                      nopen, g_movec, oskel, g_x, acc, g_ax )
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
c     Return the ROKS orbital 2e-Hessian vector product, g_ax = A * g_x
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
      integer nfock, ifock, jfock
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
      double precision time1_2e, time2_2e, time1_xc, time2_xc
      logical odebug
c
      integer dims(3), chnk(3)
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer g_dd, g_ff
      integer rtdb, calc_type
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision one, zero, mone, four, half, mhalf, two
      double precision itol_floor, itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0)
c
      double precision Exc(2)
      integer nExc
      logical xc_gotxc,xc_got2nd,rodft
      double precision xc_hfexch
      external xc_gotxc,xc_hfexch,xc_got2nd
c     
c     This get's cleaned up a lot when the Fock build accepts
c     a single multiple dimension GA for its input.
c
c     Lots more to do on memory management.
c
      call ga_inquire(g_x,gtype,vlen,nvec)
      if (nvec .gt. max_nvec) call errquit
     $     ('dft_roks_hessv_2e: too many vectors',nvec, UNKNOWN_ERR)
      odebug = util_print('dft_roks_hessv', print_debug)
c     
      tol2e = min(max(acc,itol_floor),itol_ceil)
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
      ooff = nclosed + 1
      oend = nclosed + nopen
c     
c     Allocate global for work space
c     
      g_tmp1 = ga_create_atom_blocked(geom, basis, 
     $         'dft_roks_hess_2e: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis,
     $         'dft_roks_hess_2e: tmp2')
      if (nopen .gt. 0) then
         g_tmp3 = ga_create_atom_blocked(geom, basis,
     $            'dft_roks_hess_2e: tmp3')
      endif
c
c     Allocate arrays and generate the AO densities
c     
      nfock = 0
      if (xc_gotxc()) then
        dims(1) = nvec
        chnk(1) = nvec
        if (nopen.gt.0) then
          dims(1) = 4*nvec
          chnk(1) = 4*nvec
        endif
        dims(2) = nbf
        dims(3) = nbf
        chnk(2) = -1
        chnk(3) = -1
        if (.not.nga_create(MT_DBL,3,dims,"densities",chnk,g_dd))
     +    call errquit("dft_roks_hessv_2e: failed to allocate densities"
     +                ,0,GA_ERR)
        call ga_zero(g_dd)
        if (nopen.gt.0) then
          dims(1) = 4*nvec
          chnk(1) = 4*nvec
        endif
        if (.not.nga_create(MT_DBL,3,dims,"focks",chnk,g_ff))
     +    call errquit("dft_roks_hessv_2e: failed to allocate focks",
     +                 0,GA_ERR)
      endif
      do ivec = 1, nvec
         g_dcv = ga_create_atom_blocked(geom, basis, 
     +           'dft_roks_hessv_2e: dcv')
         nfock          = nfock + 1
         iv_dens(nfock) = g_dcv
         jfac(nfock)    =  1.0d0   ! CV
         kfac(nfock)    = -0.5d0*xc_hfexch()   ! CV
         if (nopen .gt. 0) then
            g_dcp = ga_create_atom_blocked(geom, basis,
     +              'dft_roks_hessv_2e: dcp')
            g_dpv = ga_create_atom_blocked(geom, basis,
     +              'dft_roks_hessv_2e: dpv')
            nfock          = nfock + 1
            iv_dens(nfock) = g_dcp
            jfac(nfock)    =  1.0d0 ! CP
            kfac(nfock)    = -0.5d0*xc_hfexch() ! CP
            nfock          = nfock + 1
            iv_dens(nfock) = g_dpv
            jfac(nfock) =  0.0d0 ! PV
            kfac(nfock) =  0.5d0*xc_hfexch() ! PV
         endif
c     
c        Compute CV, PV & CP "densities" from argument vector
c     
         call ga_zero(g_tmp1)
         call ga_zero(g_tmp2)
         call ga_zero(g_dcv)
         call ga_vec_to_mat(g_tmp1, 1, nvir, 1, nclosed,
     $                      g_x, 1, ivec)
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
c           PV density
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
c           CP density
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
         call ga_symmetrize(g_dcv)
         if (nopen .gt. 0) then
            call ga_symmetrize(g_dpv)
            call ga_symmetrize(g_dcp)
         endif
c
c        Put the gradients of the energy in place for DFT
c
         if (xc_gotxc()) then
           alo(1) = 1
           alo(2) = 1
           ahi(1) = nbf
           ahi(2) = nbf
           blo(1) = ivec
           blo(2) = 1
           blo(3) = 1
           bhi(1) = ivec
           bhi(2) = nbf
           bhi(3) = nbf
           call nga_copy_patch('n',g_dcv,alo,ahi,g_dd,blo,bhi)
           if (nopen.gt.0) then
             call nga_add_patch(-1.0d0,g_dpv,alo,ahi,2.0d0,
     +                          g_dd,blo,bhi,g_dd,blo,bhi)
             blo(1) = ivec+nvec
             bhi(1) = ivec+nvec
             call nga_copy_patch('n',g_dcv,alo,ahi,g_dd,blo,bhi)
             call nga_add_patch(-1.0d0,g_dcp,alo,ahi,2.0d0,
     +                          g_dd,blo,bhi,g_dd,blo,bhi)
           endif
         endif
c     
c        We actually use dcp+dpv and dcp-dpv to build fock matrix
c     
         if (nopen .gt. 0) then
            call ga_dadd(one, g_dcp, one, g_dpv, g_tmp1)
            call ga_dadd(one, g_dcp,mone, g_dpv, g_tmp2)
            call ga_copy(g_tmp1, g_dcp)
            call ga_copy(g_tmp2, g_dpv)
         endif
      enddo
c
      if (xc_gotxc()) then
        call ga_zero(g_ff)
      endif
c
      if (.not. ga_destroy(g_tmp1))
     &    call errquit('dft_roks_hessv_2e: GA?',0, GA_ERR)
      if (.not. ga_destroy(g_tmp2))
     &    call errquit('dft_roks_hessv_2e: GA?',0, GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_tmp3))
     &       call errquit('dft_roks_hessv_2e: GA?',0, GA_ERR)
      endif
c
c     Fock build
c     
      do ifock = 1, nfock
         iv_fock(ifock) = ga_create_atom_blocked
     $        (geom, basis, 'dft_roks_hessv_2e:h2e fock')
         call ga_zero(iv_fock(ifock))
      enddo
c
      if(util_print('dft timings', print_high))
     &      time1_2e=util_cpusec()   ! start 2e build time
      call fock_2e( geom, basis, nfock, jfac, kfac,
     $     tol2e, oskel, iv_dens, iv_fock,.false.)
      if(util_print('dft timings', print_high))
     &      time2_2e=util_cpusec()   ! end 2e build time
c
c     print fock_2e build time
c
      if(util_print('dft timings', print_high)) then
       if (ga_nodeid().eq.0) then
         write(*,"(4x,'Fock_2e Build Time:',F13.1,'s')")
     &              time2_2e-time1_2e
       endif
      end if
c
      Exc(1) = 0.0d0
      Exc(2) = 0.0d0
      nExc = 2
      if(util_print('dft timings', print_high))
     &      time1_xc=util_cpusec()   ! start xc build time
      if (xc_gotxc()) then
c
c       Swap the calc_type to get fock_xc to do the CPHF LHS at this
c       point. Keep the original setting to restore original state
c       at the end.
c
        rtdb = bgj_get_rtdb_handle()
        if (.not. rtdb_get(rtdb, 'fock_xc:calc_type', mt_int, 1,
     &     calc_type)) then
           calc_type = 1
        endif
        if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, 2))
     &    call errquit('dft_roks_hessv_2e: could not reset calc_type',
     &                 0,UERR)
c
        if (nopen.eq.0) then
          jfock = nfock
        else
          jfock = 8*nfock/3
        endif
c        
c       include hessian contribution only if second derivatives are available
c       and if the calculation is not a rodft calculation
c
        if (.not. rtdb_get(rtdb, 'dft:rodft', mt_log, 1, rodft))
     &      rodft = .false.   ! query the rtdb for the type of calculation

        if(xc_got2nd().and. (.not. rodft)) then
           call fock_xc(geom, nbf, basis, jfock, g_dd, g_ff,
     +                  Exc, nExc,.true.)
        end if
c
        if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1,
     &                     calc_type))
     &    call errquit('dft_roks_hessv_2e: could not reset calc_type',
     &                 0,UERR)
c
        if (nopen.eq.0) then
          do ifock = 1, nfock
            alo(1) = ifock
            alo(2) = 1
            alo(3) = 1
            ahi(1) = ifock
            ahi(2) = nbf
            ahi(3) = nbf
            blo(1) = 1
            blo(2) = 1
            bhi(1) = nbf
            bhi(2) = nbf
            call nga_add_patch(one,g_ff,alo,ahi,
     +                         one,iv_fock(ifock),blo,bhi,
     +                         iv_fock(ifock),blo,bhi)
          enddo
        else
          do ifock = 1, nfock/3
            alo(1) = nfock/3+ifock
            alo(2) = 1
            alo(3) = 1
            ahi(1) = nfock/3+ifock
            ahi(2) = nbf
            ahi(3) = nbf
            blo(1) = 1
            blo(2) = 1
            bhi(1) = nbf
            bhi(2) = nbf
            jfock  = (ifock-1)*3
            call nga_add_patch(half,g_ff,alo,ahi,
     +                         one,iv_fock(jfock+1),blo,bhi,
     +                         iv_fock(jfock+1),blo,bhi)
            call nga_add_patch(one,g_ff,alo,ahi,
     +                         one,iv_fock(jfock+3),blo,bhi,
     +                         iv_fock(jfock+3),blo,bhi)
            alo(1) = 3*nfock/3+ifock
            ahi(1) = 3*nfock/3+ifock
            call nga_add_patch(half,g_ff,alo,ahi,
     +                         one,iv_fock(jfock+1),blo,bhi,
     +                         iv_fock(jfock+1),blo,bhi)
            call nga_add_patch(one,g_ff,alo,ahi,
     +                         one,iv_fock(jfock+2),blo,bhi,
     +                         iv_fock(jfock+2),blo,bhi)
          enddo
        endif
        if (.not.ga_destroy(g_dd))
     +    call errquit('dft_roks_hessv_2e: could not destroy g_dd',
     +                 0,GA_ERR)
        if (.not.ga_destroy(g_ff))
     +    call errquit('dft_roks_hessv_2e: could not destroy g_ff',
     +                 0,GA_ERR)
      endif
      if(util_print('dft timings', print_high))
     &      time2_xc=util_cpusec()   ! end xc build time
c
c     print fock_xc build time
c
      if(util_print('dft timings', print_high)) then
        if (ga_nodeid().eq.0) then
          write(*,"(4x,'Fock_xc Build Time:',F13.1,'s')")
     &                time2_xc-time1_xc
        endif
      end if
c     
      if (oskel) then
         if (oscfps) call pstat_on(ps_sym_sym)
         do ifock = 1, nfock
            call sym_symmetrize(geom, basis, .false., iv_fock(ifock))
         enddo
         if (oscfps) call pstat_off(ps_sym_sym)
      endif
c
      g_tmp1 = ga_create_atom_blocked(geom, basis,
     +         'dft_roks_hessv_2e: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis,
     +         'dft_roks_hessv_2e: tmp2')
      if (nopen .gt. 0) then
         g_tmp3 = ga_create_atom_blocked(geom, basis,
     +            'dft_roks_hessv_2e: tmp3')
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
c        P      =  4(ij|kl) - (ik|jl) - (il|kj)
c         ij,kl
c     
c        K      =  (ik|jl) + (il|kj)
c         ij,kl
c     
c         cv         cv          pv   cp
c        Z   =  2P.[D  ]  +  P.[D  + D  ]
c     
c         pv           cv           cp   pv
c        Z   =  0.5d0*Z   + 0.5*K.[D  - D  ]
c     
c         cp           cv           cp   pv
c        Z   =  0.5d0*Z   - 0.5*K.[D  - D  ]
c     
c        Add the Fock matrices together overwriting the density
c        matrices to form the results above
c     
         call ga_copy(g_tmp1, g_dcv)
         call ga_dscal(g_dcv, two)
         if (nopen .gt. 0) then
            call ga_dadd(one,  g_dcv,  one, g_tmp2, g_dcv)
            call ga_dadd(half, g_dcv, half, g_tmp3, g_dpv)
            call ga_dadd(half, g_dcv,mhalf, g_tmp3, g_dcp)
         endif
c     
c        Closed-Virtual bit
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
c     
         if (nopen .gt. 0) then
c     
c           Open-Virtual bit
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
     $           g_ax, xoff, ivec, four, '+')
c     
c           Closed-Open bit
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
         endif
      enddo
c
      do ifock = 1, nfock
         if (.not. ga_destroy(iv_dens(ifock)))
     $        call errquit('dft_roks_hessv_2e:destroy?',0, GA_ERR)
         if (.not. ga_destroy(iv_fock(ifock)))
     $        call errquit('dft_roks_hessv_2e:destroy?',0, GA_ERR)
      enddo
c
      if (.not. ga_destroy(g_tmp1))
     $     call errquit('dft_roks_hessv_2e:destroy?',0, GA_ERR)
      if (.not. ga_destroy(g_tmp2))
     $     call errquit('dft_roks_hessv_2e:destroy?',0, GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_tmp3))
     $        call errquit('dft_roks_hessv_2e:destroy?',0, GA_ERR)
      endif
c     
      end

