      subroutine gw_vxc(vxc,nmo,pars)
      implicit none
#include "errquit.fh" 
#include "mafdecls.fh"
#include "global.fh"      
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"      
#include "cdft.fh"
#include "util.fh"
#include "cgridfile.fh"
#include "dftps.fh"
#include "gw.fh"

      type(gw_params_t) :: pars
      integer,intent(in) :: nmo
      double precision,intent(out) :: vxc(nmo,*)
c
      character(len=8), parameter :: pname = 'gw_vxc: '
c
      integer g_dens(2), g_evecs(2), g_vxc(4)
      integer g_wght, g_xyz, g_nq

      logical wght_ga, ok

      integer isp, ivxc_opt, nexc, imo
      integer nscr, l_scr, k_scr, npol, natoms
      integer l_cntoce,l_cntobfr,l_cetobfr,l_rdens_atom,l_tmp
      integer k_cntoce,k_cntobfr,k_cetobfr,k_rdens_atom,k_tmp
      integer lcd_coef,icd_coef

      double precision rho_n, ecoul, exc(2)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      logical grid_reopen,xc_rep_close
      external grid_reopen,xc_rep_close
c
c      if (.not.grid_written) then
c        grid_written = grid_reopen(geom)
c      endif
c      if (grid_written) then
c        call grid_setquad(pars%rtdb)
c       call grid_packing_info
c     endif
      grid_written = .false.
      call grid_setquad(pars%rtdb)
      call grid_packing_info

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
        if (.not.ga_duplicate(g_vxc(1),g_vxc(2),'beta vxc'))
     $    call errquit(pname//'failed to create ga',0,GA_ERR)
        call ga_zero(g_vxc(2))
        if (.not.ga_duplicate(g_vxc(1),g_dens(2),'beta density'))
     $    call errquit(pname//'failed to create ga',0,GA_ERR)
        call ga_zero(g_dens(2))
      endif

      if (.not.ma_push_get(mt_int,natoms*2,'cetobfr map',l_cetobfr,
     $                     k_cetobfr))
     $  call errquit(pname//'failed to allocate cetobfr',0,MA_ERR)
      if (.not.ma_push_get(mt_int,nshells_ao,'cntoce map',l_cntoce,
     $                     k_cntoce))
     $  call errquit(pname//'failed to allocate cntoce',0,MA_ERR)
      if (.not.ma_push_get(mt_int,nshells_ao*2,'cntobfr map',l_cntobfr,
     $                     k_cntobfr))
     $  call errquit(pname//'failed to allocate cntobfr',0,MA_ERR)
      call build_maps(ao_bas_han, int_mb(k_cntoce), int_mb(k_cntobfr),
     $                int_mb(k_cetobfr), natoms, nshells_ao)
      if (.not.ma_chop_stack(l_cntoce))
     $  call errquit(pname//'failed to chop stack',0,MA_ERR) 
 
      if(.not.ma_push_get(mt_dbl,ipol*natoms**2,'rdens atom',
     $                    l_rdens_atom,k_rdens_atom))
     $  call errquit(pname//'failed to allocate rdens_atom',0,MA_ERR)
      call dfill(ipol*natoms**2,0.0d0,dbl_mb(k_rdens_atom),1)
c
c     density
c
      do isp=1,ipol
        if (odftps) call pstat_on(ps_dgemm)
        call ga_dgemm('n','t',nbf_ao,nbf_ao,noc(isp),2d0/dble(ipol),
     $                 pars%g_movecs(isp),pars%g_movecs(isp),0d0,
     $                 g_dens(isp))
        if (odftps) call pstat_off(ps_dgemm)
      enddo
c
      nscr = 2*nbf_ao_mxnbf_ce**2
      if (.not.ma_push_get(mt_dbl,nscr,'scr',l_scr,k_scr))
     $  call errquit(pname//'failed to allocate scr',0,MA_ERR) 
      call util_ga_mat_reduce(nbf_ao, natoms, int_mb(k_cetobfr), g_dens,
     $                        ipol, dbl_mb(k_rdens_atom), 'absmax',
     $                        dbl_mb(k_scr),nbf_ao_mxnbf_ce,.true.) 

      if (.not.ma_pop_stack(l_scr))
     $  call errquit(pname//'failed to pop stack',0,MA_ERR)
c
      rho_n = 0.0d0
      cdfit = .true.
      npol = (ipol*(ipol+1))/2
      exc = 0d0
c
c      call xc_rep_init(pars%rtdb,geom,ao_bas_han,g_dens,g_vxc,nbf_ao,
c     $                 ipol,.true.,.false.) 
c
      if (.not.MA_Push_Get(MT_Dbl,4*nbf_cd,'cd_coef',lcd_coef,
     &     icd_coef))
     &    call errquit('gw_vxc: cannot allocate cd_coef',0, MA_ERR)
      call xc_getv(pars%rtdb, exc, ecoul, 2, 0, 0,
     $             g_dens, g_vxc, .false., g_wght, g_xyz, g_nq,
     $             .false., rho_n, dbl_mb(k_rdens_atom),
     $             dbl_mb(k_cetobfr), natoms, dbl_mb(icd_coef)) 
      if (.not.ma_pop_stack(lcd_coef))
     &     call errquit('gw_vxc: cannot pop stack',0, MA_ERR)
      pars%exscf = exc(1)
      pars%ecscf = exc(2)

c      if (.not.xc_rep_close(pars%rtdb,nbf_ao,ipol,g_vxc,g_dens,.false.))
c     $  call errquit(pname//'xc_rep_close failed',0,DISK_ERR) 
c
      if (.not.ma_push_get(mt_dbl,2*nbf_ao,'vxc temp',l_tmp,k_tmp))
     $  call errquit(pname//'failed to allocate temp',0,MA_ERR) 
      do isp=1,ipol
        call ga_dgemm('n','n',nbf_ao,nmo,nbf_ao,1.0d0,g_vxc(isp),
     $                 pars%g_movecs(isp),0.0d0,g_dens(isp))
        do imo=1,nmo
          vxc(imo,isp) = ga_ddot_patch(g_dens(isp),'n',1,nbf_ao,imo,imo,
     $                          pars%g_movecs(isp),'n',1,nbf_ao,imo,imo)
        enddo
      enddo
c
      if (.not.ma_chop_stack(l_cetobfr))
     $  call errquit(pname//'failed to chop stack',0,MA_ERR)
      ok = .true.
      do isp=1,ipol
        ok =  ok .and. ga_destroy(g_dens(isp))
        ok = ok .and. ga_destroy(g_vxc(isp))
      enddo
      if (.not.ok) call errquit(pname//'failed to destroy ga',0,GA_ERR)
      call grid_cleanup(.false.)
c 
      end subroutine
