C
C $Id: onedd_cont.F,v 1.2 1999-12-21 23:52:40 windus Exp $
C
      subroutine onedd_cont(rtdb)
c
c  Add the one electron derivative (first and seconds) contributions to 
c  the hessian matrix and the RHS
c  Start with RHF and work our way up.
c
      implicit none
c
#include "bas.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
c
      integer rtdb  ! [input] RTDB handle
c
      integer g_wdens, ndens
      integer max1e, max2e, mscratch_1e, mscratch_2e
      integer l_buf, k_buf, l_scr, k_scr
      integer matb, lsqatom, l_dens, k_dens, l_wdens, k_wdens
      logical doS, doT, doV
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
c Get memory for the integral buffers
c
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      write(6,*) 'max1e and mscratch_1e ',max1e,mscratch_1e
      call flush(6)
      if (.not. ma_push_get(mt_dbl,max1e,'intdd buffer',l_buf,k_buf))
     *  call errquit('onedd_cont:could not get intdd buffer',max1e)
      if (.not. ma_push_get(mt_dbl,mscratch_1e,'intdd scratch',
     *  l_scr, k_scr))call errquit('onedd_cont: scratch alloc failed',
     *  mscratch_1e)
c
c Need to form the weighted density matrix to multiply by the S 
c second derivative integrals.  Also get temporary working space on
c each node.
c
      g_wdens = ga_create_atom_blocked(geom, basis, 'weighted density')
      ndens = 1
      call get_wdens(g_dens,ndens,theory,scftype,rtdb,g_wdens,nbf)
      write(6,*) 'The weighted density is'
      call ga_print(g_wdens)
c
      if (.not. bas_nbf_ce_max(basis,matb))
     *  call errquit('onedd_cont: could not get matb',555)
      lsqatom = matb * matb
      if (.not. ma_push_get(mt_dbl,lsqatom,'piece_density',
     *  l_dens,k_dens))
     *  call errquit('onedd_cont:could not allocate l_dens',lsqatom)
      if (.not. ma_push_get(mt_dbl,lsqatom,'piece_wdensity',l_wdens,
     *     k_wdens)) call errquit
     *    ('onedd_cont:could not allocate l_wdens',lsqatom) 
c
c  Call the actual routine to add the S, V and T contributions
c
      doS = .false.
      doV = .true.
      doT = .true.
      call onedd_calc(dbl_mb(k_exy),ncent*3,g_rhs,ncent*3,dbl_mb(k_buf),
     *  max1e,dbl_mb(k_scr),mscratch_1e,skelsym,
     *  doS,doT,doV,g_wdens,g_dens,geom,basis,ncent,dbl_mb(k_dens),
     *  dbl_mb(k_wdens),matb)
c
c should put the next line in when done with some testing
c     call sym_hess_symmetrize(geom,dbl_mb(k_exy))
      write(6,*) 'The Hessian after derivative S contributions:'
      call hnd_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
c
c Get rid of the weighted density and the other temporary arrays
c
      if (.not. ma_pop_stack(l_wdens)) 
     *  call errquit('onedd_cont: failed to pop l_wdens',555)
      if (.not. ma_pop_stack(l_dens)) 
     *  call errquit('onedd_cont: failed to pop l_dens',555)
      if (.not. ma_pop_stack(l_scr)) 
     *  call errquit('onedd_cont: failed to pop l_scr',555)
      if (.not. ma_pop_stack(l_buf)) 
     *  call errquit('onedd_cont: failed to pop l_buf',555)
      if (.not. ga_destroy(g_wdens)) 
     *  call errquit('onedd_cont: failed to free g_wdens',555)
c
      end
