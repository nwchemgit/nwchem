C
C $Id: onedd_cont.F,v 1.28 2003-10-17 22:55:13 carlfahl Exp $
C
      subroutine onedd_cont(rtdb)
c
c  Add the one electron derivative (first and seconds) contributions to 
c  the hessian matrix and the RHS
c  Start with RHF and work our way up.
c
      implicit none
#include "errquit.fh"
c
#include "bas.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
c
      integer rtdb  ! [input] RTDB handle
c
      integer max1e, max2e, mscratch_1e, mscratch_2e
      integer l_buf, k_buf, l_scr, k_scr
      integer matb, lsqatom, l_dens, k_dens, l_wdens, k_wdens
      integer alo(3), ahi(3), blo(3), bhi(3), idummy
      double precision one
      data one /1.0D+00/
      logical doS, doT, doV, doECP
      logical oprint, oneprint
c
      character*255 tmpname
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
c Set up print info
c
      oprint = util_print('hess_follow',print_high)
      oprint = oprint .and. (ga_nodeid().eq.0)
      oneprint = util_print('onedd_cont',print_debug)
c
c Get memory for the integral buffers
c
      if (oprint) write(LuOut,*) 'HESS: getting memory'
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      if (.not. ma_push_get(mt_dbl,max1e,'intdd buffer',l_buf,
     *  k_buf))
     *  call errquit('onedd_cont:could not get intdd buffer',max1e,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,mscratch_1e,'intdd scratch',
     *  l_scr, k_scr))call errquit('onedd_cont: scratch alloc failed',
     *  mscratch_1e, MA_ERR)
c
      if (oneprint) then
        if(ga_nodeid().eq.0)write(LuOut,*)
        if(ga_nodeid().eq.0)write(LuOut,*) 'The weighted density is'
        call ga_print(g_wdens)
        if(ga_nodeid().eq.0)write(LuOut,*)
        if(ga_nodeid().eq.0)write(LuOut,*) 'The density is'
        call ga_print(g_dens(1))
        if (ndens.gt.1) then
          if(ga_nodeid().eq.0)write(LuOut,*)
          if(ga_nodeid().eq.0)write(LuOut,*) 'the beta density is '
          call ga_print(g_dens(2))
        endif
      endif
c
c Set up some local memory to hold pieces of the density and weighted density
c
      if (.not. bas_nbf_ce_max(basis,matb))
     *  call errquit('onedd_cont: could not get matb',555, BASIS_ERR)
      lsqatom = matb * matb
      if (.not. ma_push_get(mt_dbl,ndens*lsqatom,'piece_density',
     *  l_dens,k_dens))
     *  call errquit('onedd_cont:could not allocate l_dens',lsqatom,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'piece_wdensity',l_wdens,
     *     k_wdens)) call errquit
     *    ('onedd_cont:could not allocate l_wdens',lsqatom, MA_ERR)
c
c  first print out the rhs values (should be all zeroes)!
c
      if (oneprint) then
        if(ga_nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*) 'g_rhs before going into onedd_calc'
        write(LuOut,*) 'These values should all be zero'
        endif
        call ga_print(g_rhs)
      endif
      if (oprint.and.(ga_nodeid().eq.0)) 
     *     write(LuOut,*) 'HESS: before the onedd_calc call'
c
      doS = .true. 
      doV = .true.
      doT = .true.
      doECP = .false.
      doECP = bas_get_ecp_handle(basis,idummy)
c
c  Call the actual routine to add the second derivative S, V and T contributions
c
      call onedd_calc(dbl_mb(k_exy),n3xyz,dbl_mb(k_buf),
     *  max1e,dbl_mb(k_scr),mscratch_1e,skelsym,
     *  doS,doT,doV,doECP,g_wdens,g_dens,geom,basis,ncent,ndens,
     *  dbl_mb(k_dens),dbl_mb(k_wdens),matb)
c
      if(skelsym)call sym_hess_symmetrize(geom,dbl_mb(k_exy))
      if (oneprint.and.(ga_nodeid().eq.0)) then
        write(LuOut,*)
        call util_flush(LuOut)
        write(LuOut,*) 'The Hessian after onedd contributions:'
        call hess_hssout(dbl_mb(k_exy),n3xyz,n3xyz,n3xyz)
        write(LuOut,*)
      endif
c
c Get rid of several temporary arrays
c
      if (.not. ma_pop_stack(l_wdens)) 
     *  call errquit('onedd_cont: failed to pop l_wdens',555, MA_ERR)
      if (.not. ma_pop_stack(l_dens)) 
     *  call errquit('onedd_cont: failed to pop l_dens',555, MA_ERR)
c
c Now handle the Hx and Sx terms for the rhs
c First handle Sx term
c
      if (oprint) write(LuOut,*) 'HESS: getting the Sx term'
      doS = .true.
      doV = .false.
      doT = .false.
      call oned_calc(n3xyz,g_rhs,dbl_mb(k_buf),
     *  max1e,dbl_mb(k_scr),mscratch_1e,.false.,
!     *  max1e,dbl_mb(k_scr),mscratch_1e,skelsym,
     *  doS,doT,doV,geom,basis,ncent)
      if (ndens.gt.1) then     ! UHF case
        alo(1) = 1
        ahi(1) = n3xyz
        alo(2) = 1
        ahi(2) = nbf
        alo(3) = 1
        ahi(3) = nbf
        blo(1) = n3xyz + 1
        bhi(1) = 2 * n3xyz
        blo(2) = 1
        bhi(2) = nbf
        blo(3) = 1
        bhi(3) = nbf
        call nga_copy_patch('N', g_rhs, alo, ahi, g_rhs, blo, bhi)
      endif
c
      if (oneprint) then
        if(ga_nodeid().eq.0)write(LuOut,*)
        if(ga_nodeid().eq.0)write(LuOut,*) 'Sx term'
        call ga_print(g_rhs)
      endif
c
c write out Sx terms to disk
c
      call sx_write(g_rhs,n3xyz*ndens,nbf,sx_dra_handle)
c
c Clear g_rhs again
c
      call ga_zero(g_rhs)
c
      if (oprint) write(LuOut,*) 'HESS: getting the Hx term'
      doS = .false.
      doV = .true.
      doT = .true.
      call oned_calc(n3xyz,g_rhs,dbl_mb(k_buf),
     *  max1e,dbl_mb(k_scr),mscratch_1e,.false.,
     *  doS,doT,doV,geom,basis,ncent)
      if (ndens.gt.1) call nga_copy_patch('N', g_rhs, 
     *  alo, ahi, g_rhs, blo, bhi)
c
c  Print out the rhs values if requested
c
      if (oneprint) then
        if(ga_nodeid().eq.0)write(LuOut,*)
        if(ga_nodeid().eq.0)write(LuOut,*) 'g_rhs'
        call ga_print(g_rhs)
      endif
c
c  Get rid of the temporary buffers for the integrals
c
      if (.not. ma_pop_stack(l_scr)) 
     *  call errquit('onedd_cont: failed to pop l_scr',555, MA_ERR)
      if (.not. ma_pop_stack(l_buf)) 
     *  call errquit('onedd_cont: failed to pop l_buf',555, MA_ERR)
c
c  Get the dipole derivatives contribution
c  Reuse g_wdens for the total density
c
      if (ndens.gt.1) then ! want total density
        call ga_dadd (one, g_dens(1), one, g_dens(2), g_wdens)
      else
        call ga_copy (g_dens(1), g_wdens)
      endif
      if (oprint) write(LuOut,*) 
     *  'HESS: getting the dipole derivatives'
      call hess_dipder(geom,basis, dbl_mb(k_dder),
     *  g_wdens,n3xyz,ncent)
      if (oprint) then
        write(luout,*) 'after dipole derivative'
        write(tmpname,*) ' derivative dipole '
        call stpr_print_ddipole(dbl_mb(k_dder),tmpname,ncent,-1.0d00)
      endif
c
      return
      end
