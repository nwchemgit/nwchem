      subroutine gw_vxc(pars)
      implicit none
#include "errquit.fh" 
#include "mafdecls.fh"
#include "global.fh"      
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"      
#include "cdft.fh"
#include "util.fh"
#include "gw.fh"

      type(gw_params_t) :: pars
      character(len=8), parameter :: pname = 'gw_vxc: '
c
      integer g_dens(2), g_evecs(2), g_vxc(4)
      integer g_wght, g_xyz, g_nq

      logical wght_ga

      integer isp, ivxc_opt, nexc
      integer nscr, l_scr, k_scr, npol, natoms
      integer l_cntoce,l_cntobfr,l_cetobfr,l_rdens_atom
      integer k_cntoce,k_cntobfr,k_cetobfr,k_rdens_atom

      double precision rho_n, ecoul, exc(1)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      call grid_setquad(pars%rtdb)
c
      if (.not.geom_ncent(geom,natoms))
     $  call errquit(pname//'geom_ncent failed',0,GEOM_ERR)
c      
      g_vxc(1) = ga_create_atom_blocked(geom,ao_bas_han,'vxc')
      call ga_zero(g_vxc(1))
      if (.not.ga_duplicate(g_vxc(1),g_dens(1),'alpha density'))
     $  call errquit(pname//'failed to create ga',0,GA_ERR)
      call ga_zero(g_dens(1))
      if (ipol.ne.1) then
        if (.not.ga_duplicate(g_vxc(1),g_vxc(2),'alpha density'))
     $    call errquit(pname//'failed to create ga',0,GA_ERR)
        call ga_zero(g_vxc(2))
        if (.not.ga_duplicate(g_vxc(1),g_dens(2),'alpha density'))
     $    call errquit(pname//'failed to create ga',0,GA_ERR)
        call ga_zero(g_dens(2))
      endif

      if (.not.ma_push_get(mt_int,nshells_ao,'cntoce map',l_cntoce,
     $                     k_cntoce))
     $  call errquit(pname//'failed to allocate cntoce',0,MA_ERR)
      if (.not.ma_push_get(mt_int,nshells_ao*2,'cntobfr map',l_cntobfr,
     $                     k_cntobfr))
     $  call errquit(pname//'failed to allocate cntobfr',0,MA_ERR)
      if (.not.ma_push_get(mt_int,natoms*2,'cetobfr map',l_cetobfr,
     $                     k_cetobfr))
     $  call errquit(pname//'failed to allocate cetobfr',0,MA_ERR)
      call build_maps(ao_bas_han, int_mb(k_cntoce), int_mb(k_cntobfr),
     $                int_mb(k_cetobfr), natoms, nshells_ao)

 
      if(.not.ma_push_get(mt_dbl,2*natoms**2,'rdens atom',l_rdens_atom,
     $                    k_rdens_atom))
     $  call errquit(pname//'failed to allocate rdens_atom',0,MA_ERR)


      do isp=1,ipol
        call ga_dgemm('n','t',nbf_ao,nbf_ao,noc(1),2d0/dble(ipol),
     $                 pars%g_movecs(isp),pars%g_movecs(isp),0d0,
     $                 g_dens(isp))
      enddo
c
      nscr = nbf_ao_mxnbf_ce**2
      if (.not.ma_push_get(mt_dbl,nscr,'scr',l_scr,k_scr))
     $  call errquit(pname//'failed to allocate scr',0,MA_ERR) 
      call util_ga_mat_reduce(nbf_ao, natoms, dbl_mb(k_cetobfr), g_dens,
     $                        ipol, dbl_mb(k_rdens_atom), 'absmax',
     $                        dbl_mb(k_scr),nbf_ao_mxnbf_ce,.true.) 
      if (.not.ma_pop_stack(l_scr))
     $  call errquit(pname//'failed to pop stack',0,MA_ERR)
c
      rho_n = 0.0d0
      cdfit = .true.
      npol = (ipol*(ipol+1))/2
c
      call xc_getv(pars%rtdb, exc, ecoul, 1, 0, 0,
     $             g_dens, g_vxc, .false., g_wght, g_xyz, g_nq,
     $             .false., rho_n, dbl_mb(k_rdens_atom),
     $             dbl_mb(k_cetobfr), natoms) 
c
      end subroutine
