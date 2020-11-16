      subroutine gw_vxc(pars)
      implicit none
#include "cdft.fh"
#include "gw.fh"

      type(gw_params_t) :: pars
c
c      call grid_setquad(pars%rtdb)
c      if (hashfx)
c        call fock_2e(geom, ao_bas_han, 1, jfac, kfac,
c     $               pars%tol2e, .false., g_dens(1), g_vxc(1), .false.)
c      endif
c      call xc_getvxc(pars%rtdb, Exc, 1, iVxc_opt, 0, g_dens,
c     $               g_vxc, .false., rho_n, 0, 0)
c
      end subroutine
