      subroutine dft_uks_hessv_2e(acc, g_x, g_ax, nvec)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
c     
      double precision acc      ! [input] required accuracy of products
      integer g_x               ! [input] handle to input vectors
      integer g_ax              ! [input] handle to output products
      integer nvec              ! [input] No. of vectors
c
c     SAME spin -> Bai,bj Xbj = [4 (ai|bj) - 2 (bi|aj) - 2 (ij|ab)] Xbj
c                             = 4 (Jai - Kai)
c
c     where J and K constructed from symmetrized AO density CXCT and
c     transformed into MO basis of the SAME spin.
c   
c     DIFF spin -> Bai,bj Xbj = 4 (ai|bj) Xbj
c                             = 4 Jai
c     
c     where J constructed from  symmetrized AO density CXCT and J 
c     transformed into MO basis of the OPPOSITE spin.
c
      logical oskel_local
      integer ivec, nfock
      integer g_tmp1, g_tmp2, g_tmp3, g_tmp4
      integer nocc, nvir, g_vecs, nocc_diff, nvir_diff, g_vecs_diff
      integer iset, isetoff, i
      integer istart, istart_diff, iend, iend_diff, vbase, vbase_diff
      integer max_nvec
      parameter (max_nvec=300)
      integer g_dens(4,max_nvec), g_fock(4,max_nvec), g_d, g_j, g_k
      double precision jfac(4,max_nvec), kfac(4,max_nvec)
      double precision tol2e_local
      double precision time1_2e, time2_2e, time1_xc, time2_xc
c
      integer dims(3), chnk(3)
      integer alo(3), ahi(3)
      integer blo(3), bhi(3)
      integer clo(3), chi(3)
      integer g_dd, g_ff
      integer rtdb, calc_type
      integer nExc
      double precision Exc(2)
      integer nocc1, nvir1, nocc_diff1, nvir_diff1
c
      logical xc_gotxc,xc_got2nd
      double precision xc_hfexch
      external xc_gotxc,xc_hfexch,xc_got2nd
c
      double precision itol_floor, itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
c      logical ga_create_atom_blocked
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision zero, one, mone, four
      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0, four=4.0d0)
c
c     If needed create the DFT density and fock matrices.
c
      if (xc_gotxc()) then
        dims(1) = 4*nvec
        dims(2) = nbf
        dims(3) = nbf
        chnk(1) = 4*nvec
        chnk(2) = -1
        chnk(3) = -1
        if (.not.nga_create(MT_DBL,3,dims,"densities",chnk,g_dd))
     +    call errquit("dft_uks_hessv_2e: failed to allocate densities",
     +                 0,GA_ERR)
        call ga_zero(g_dd)
        dims(1) = 4*nvec
        chnk(1) = 4*nvec
        if (.not.nga_create(MT_DBL,3,dims,"focks",chnk,g_ff))
     +    call errquit("dft_uks_hessv_2e: failed to allocate focks",
     +                 0,GA_ERR)
        call ga_zero(g_ff)
      endif
c
c     First allocated arrays and build the AO densities.
c     This will be cleaner when we restructure the fock build to take
c     a single density array argument rather than an array of arrays
c
      do iset = 1, 2
         if (iset .eq. 1) then
            isetoff = 0
            istart = 1
            nocc = nalpha
         else
            isetoff = 2
            istart = nalpha*(nmo-nalpha) + 1
            nocc = nbeta
         endif
         g_vecs = g_movecs(iset)
         nvir  = nmo - nocc
         iend = istart + nocc*nvir - 1
         vbase = nocc+1         ! First virtual
c
c        If we have no electrons of a particular spin the code should
c        still work. E.g. for H atom the beta-spin nocc equals 0.
c
         nocc1 = max(1,nocc)
         nvir1 = max(1,nvir)
*ga:1:0
         if (.not.ga_create(MT_DBL,nvir1,nocc1,'dft_uks_hessv_2e: tmp1',
     $                      0, 0, g_tmp1))
     $            call errquit('dft_uks_hessv_2e: tmp1', 0, GA_ERR)
*ga:1:0
         if (.not.ga_create(MT_DBL,nvir1,nbf,'dft_uks_hessv_2e: tmp2',
     $                      0, 0, g_tmp2))
     $            call errquit('dft_uks_hessv_2e: tmp2', 0, GA_ERR)
         do ivec = 1, nvec
            g_d = ga_create_atom_blocked(geom, basis,
     $                                   'dft_uks_hessv_2e: dens')
            g_j = ga_create_atom_blocked(geom, basis,
     $                                   'dft_uks_hessv_2e: dens')
            g_k = ga_create_atom_blocked(geom, basis,
     $                                   'dft_uks_hessv_2e: dens')
            call ga_zero(g_d)
            call ga_zero(g_j)
            call ga_zero(g_k)
            g_dens(isetoff+1,ivec) = g_d
            g_dens(isetoff+2,ivec) = g_d
            g_fock(isetoff+1,ivec) = g_j
            g_fock(isetoff+2,ivec) = g_k
            if (nocc*nvir .gt. 0) then
               call ga_copy_patch('n', ! Reshape vector into matrix Xbj
     $              g_x, istart, iend, ivec, ivec,
     $              g_tmp1, 1, nvir, 1, nocc)
               call ga_matmul_patch('n', 't', one, zero,
     $              g_tmp1, 1, nvir, 1, nocc,
     $              g_vecs, 1, nocc, 1, nbf,
     $              g_tmp2, 1, nvir, 1, nbf)
               call ga_matmul_patch('n', 'n', one, zero,
     $              g_vecs, 1, nbf, vbase, nmo,
     $              g_tmp2, 1, nvir, 1, nbf,
     $              g_d, 1, nbf, 1, nbf)
               call ga_symmetrize(g_d)
            endif
c
c           Copy density for DFT part
c
            if (xc_gotxc()) then
              alo(1) = 1
              alo(2) = 1
              ahi(1) = nbf
              ahi(2) = nbf
              blo(1) = (iset-1)*nvec+ivec
              blo(2) = 1
              blo(3) = 1
              bhi(1) = (iset-1)*nvec+ivec
              bhi(2) = nbf
              bhi(3) = nbf
              call nga_copy_patch('n',g_d,alo,ahi,g_dd,blo,bhi)
            endif
         enddo
         if (.not. ga_destroy(g_tmp1))
     $       call errquit('dft_uks_hessv_2e: GA?',0, GA_ERR)
         if (.not. ga_destroy(g_tmp2))
     $       call errquit('dft_uks_hessv_2e: GA?',0, GA_ERR)
      enddo
c
c     Now make all the J & K matrices with one pass thru the integrals
c
c     NOTE that symmetrization is not yet implemented except for
c     totally symmetric products.  Assume for the time being that
c     if only 1 RHS is being requested then use symmtery, but disable
c     it if there is more than 1 RHS.
c
      do ivec = 1, nvec
         jfac(1,ivec) = one
         jfac(2,ivec) = zero
         jfac(3,ivec) = one
         jfac(4,ivec) = zero
         kfac(1,ivec) = zero
         kfac(2,ivec) = one*xc_hfexch()
         kfac(3,ivec) = zero
         kfac(4,ivec) = one*xc_hfexch()
      enddo
      oskel_local = oskel .and. (nvec.eq.1)
c
      tol2e_local = min(max(acc,itol_floor),itol_ceil)
c
      nfock = 4*nvec
c
      if(util_print('dft timings', print_high))
     &      time1_2e=util_cpusec()   ! start 2e build time
      call fock_2e(geom, basis, nfock, jfac, kfac,
     $     tol2e_local, oskel_local, g_dens, g_fock, .false.)
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
c     Calculate the DFT parts
c
      Exc(1) = zero
      Exc(2) = zero
      nExc   = 2
      if(util_print('dft timings', print_high))
     &      time1_xc=util_cpusec()   ! start xc build time
      if (xc_gotxc()) then
        rtdb = bgj_get_rtdb_handle()
        if (.not.rtdb_get(rtdb,'fock_xc:calc_type',MT_INT,1,calc_type))
     +      calc_type = 1
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',MT_INT,1,2))
     +      call errquit('dft_uks_hessv_2e: could not reset calc_type',
     +                   0,UERR)
c
c       include hessian contribution only if second derivatives are available
        if(xc_got2nd()) then
           call fock_xc(geom,nbf,basis,nfock,g_dd,g_ff,Exc,nExc,.true.)
        end if
c
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',MT_INT,1,calc_type))
     +      call errquit('dft_uks_hessv_2e: could not reset calc_type',
     +                   0,UERR)
        isetoff = 0
        do iset = 1, 2
          do ivec = 1, nvec
            alo(1) = 2*(iset-1)*nvec+nvec+ivec
            alo(2) = 1
            alo(3) = 1
            ahi(1) = 2*(iset-1)*nvec+nvec+ivec
            ahi(2) = nbf
            ahi(3) = nbf
            blo(1) = 1
            blo(2) = 1
            bhi(1) = nbf
            bhi(2) = nbf
            g_k = g_fock(isetoff+2,ivec)
            call nga_add_patch(one,g_ff,alo,ahi,
     +                         one,g_k,blo,bhi,g_k,blo,bhi)
          enddo
          isetoff = isetoff + 2
        enddo
        if (.not.ga_destroy(g_dd))
     +      call errquit('dft_uks_hessv_2e: could not destroy g_dd',
     $                   0,GA_ERR)
        if (.not.ga_destroy(g_ff))
     +      call errquit('dft_uks_hessv_2e: could not destroy g_ff',
     $                   0,GA_ERR)
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
c     Again, note that this is NOT what you want except for
c     a totally symmetric input vector
c
      if (oskel_local) then
         if (oscfps) call pstat_on(ps_sym_sym)
         do ivec = 1, nvec
            do iset = 1, 4
               call sym_symmetrize(geom, basis, .false.,
     $              g_fock(iset,ivec))
            enddo
         enddo
         if (oscfps) call pstat_off(ps_sym_sym)
      endif
c
      do ivec = 1, nvec         ! Don't need the AO densities any more
         if (.not. ga_destroy(g_dens(1,ivec)))
     $        call errquit('dft_uks_hessv_2e: GA corrupt?',0, GA_ERR)
         if (.not. ga_destroy(g_dens(3,ivec)))
     $        call errquit('dft_uks_hessv_2e: GA corrupt?',0, GA_ERR)
      enddo
c
c     Now transform to the MO basis and assemble the results
c
      do iset = 1, 2
         if (iset .eq. 1) then
            isetoff = 0
            istart = 1
            nocc = nalpha
            g_vecs = g_movecs(1)
            istart_diff = nalpha*(nmo-nalpha) + 1
            nocc_diff = nbeta
            g_vecs_diff = g_movecs(2)
         else
            isetoff = 2
            istart = nalpha*(nmo-nalpha) + 1
            nocc = nbeta
            g_vecs = g_movecs(2)
            istart_diff = 1
            nocc_diff = nalpha
            g_vecs_diff = g_movecs(1)
        endif
        nvir  = nmo - nocc
        nvir_diff = nmo - nocc_diff
        iend = istart + nocc*nvir - 1
        iend_diff = istart_diff + nocc_diff*nvir_diff - 1
        vbase = nocc+1          ! First virtual
        vbase_diff = nocc_diff+1          ! First virtual
c
c       In some special case nocc, nvir or their diff variants might
c       be zero. To ensure that the ga_creates won't fail ensure we
c       always use dimensions greater than zero.
c
        nocc1 = max(1,nocc)
        nvir1 = max(1,nvir)
        nocc_diff1 = max(1,nocc_diff)
        nvir_diff1 = max(1,nvir_diff)
*ga:1:0
        if (.not. ga_create(MT_DBL, nbf, nocc1,
     $                      'dft_uks_hessv_2e: tmp1', 0, 0, g_tmp1))
     $            call errquit('dft_uks_hessv_2e: tmp1', 0, GA_ERR)
*ga:1:0
        if (.not. ga_create(MT_DBL, nvir1, nocc1,
     $                      'dft_uks_hessv_2e: tmp2', 0, 0, g_tmp2))
     $            call errquit('dft_uks_hessv_2e: tmp2', 0, GA_ERR)
*ga:1:0
        if (.not. ga_create(MT_DBL, nbf, nocc_diff1,
     $                      'dft_uks_hessv_2e: tmp3', 0, 0, g_tmp3))
     $            call errquit('dft_uks_hessv_2e: tmp3', 0, GA_ERR)
*ga:1:0
        if (.not. ga_create(MT_DBL,nvir_diff1,nocc_diff1,
     $                      'dft_uks_hessv_2e: tmp4', 0, 0, g_tmp4))
     $            call errquit('dft_uks_hessv_2e: tmp4', 0, GA_ERR)
        do ivec = 1, nvec
           g_j = g_fock(isetoff+1,ivec)
           g_k = g_fock(isetoff+2,ivec)
c     
c          Same spin 2-e contributions
c     
           if (nocc*nvir .gt. 0) then
              call ga_dadd(one, g_j, mone, g_k, g_k) ! K <-- J-K
              call ga_matmul_patch('n', 'n', one, zero,
     $             g_k, 1, nbf, 1, nbf,
     $             g_vecs, 1, nbf, 1, nocc,
     $             g_tmp1, 1, nbf, 1, nocc)
              call ga_matmul_patch('t', 'n', one, zero,
     $             g_vecs, vbase, nmo, 1, nbf,
     $             g_tmp1, 1, nbf, 1, nocc,
     $             g_tmp2, 1, nvir, 1, nocc)
              call ga_dadd_patch(
     $                four, g_tmp2, 1, nvir, 1, nocc,
     $                one, g_ax, istart, iend, ivec, ivec,
     $                g_ax, istart, iend, ivec, ivec)
           endif
c     
c          Different spin
c
           if (nocc_diff*nvir_diff .gt. 0) then
              call ga_matmul_patch('n', 'n', one, zero,
     $             g_j, 1, nbf, 1, nbf,
     $             g_vecs_diff, 1, nbf, 1, nocc_diff,
     $             g_tmp3, 1, nbf, 1, nocc_diff)
              call ga_matmul_patch('t', 'n', one, zero,
     $             g_vecs_diff, vbase_diff, nmo, 1, nbf,
     $             g_tmp3, 1, nbf, 1, nocc_diff,
     $             g_tmp4, 1, nvir_diff, 1, nocc_diff)
              call ga_dadd_patch(
     $             four, g_tmp4, 1, nvir_diff, 1, nocc_diff,
     $             one, g_ax, istart_diff, iend_diff, ivec, ivec,
     $                  g_ax, istart_diff, iend_diff, ivec, ivec)
           endif
 200       continue
           if (.not. ga_destroy(g_j))
     $         call errquit('dft_uks_hessv_2e: GA?',0, GA_ERR)
           if (.not. ga_destroy(g_k))
     $         call errquit('dft_uks_hessv_2e: GA?',0, GA_ERR)
        enddo
        if (.not. ga_destroy(g_tmp1))
     $      call errquit('dft_uks_hessv_2e: GA?',0, GA_ERR)
        if (.not. ga_destroy(g_tmp2))
     $      call errquit('dft_uks_hessv_2e: GA?',0, GA_ERR)
        if (.not. ga_destroy(g_tmp3))
     $      call errquit('dft_uks_hessv_2e: GA?',0, GA_ERR)
        if (.not. ga_destroy(g_tmp4))
     $      call errquit('dft_uks_hessv_2e: GA?',0, GA_ERR)
      enddo
c
      end
c $Id$
