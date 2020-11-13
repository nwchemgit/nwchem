      subroutine gw_vxc(pars)
      implicit none
#include "gw.fh"

      type(gw_params_t) :: pars
c
      call grid_setquad(pars%rtdb)
      if (hashfx)
        call fock_2e(pars%geom, pars%ao_bas_han, 1, jfac, kfac,
     $               pars%tol2e, .false., g_dens(1), g_vxc(1), .false.)
      endif
      call xc_getvxc(pars%rtdb, Exc, 1, iVxc_opt, 0, g_dens,
     $               g_vxc, .false., rho_n, 0, 0)
c
      end subroutine
