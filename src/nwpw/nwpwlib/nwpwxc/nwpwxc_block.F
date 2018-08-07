C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_block.F
C> Functional initialization
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief One time initialization of the density functional library
C> data structures
C>
C> Initialize the density functional library data structures ensuring
C> that initially there are no functional terms and all weights are
C> zero.
C>
      block data nwpwxc_block
      implicit none
#include "nwpwxcP.fh"
#include "nwpwxc_param.fh"
      data nwpwxc_depends_on_rho   / .false. /
      data nwpwxc_depends_on_gamma / .false. /
      data nwpwxc_depends_on_tau   / .false. /
      data nwpwxc_cam_srhfex       / .false. /
      data nwpwxc_on               / .false. /
      data nwpwxc_num_f      / 0 /
      data nwpwxc_usr_num    / 0 /
      data nwpwxc_wght_hf    / 0.0d0 /
      data nwpwxc_wght_mp2   / 0.0d0 /
      data nwpwxc_wght_cs00  / 0.0d0 /
      data nwpwxc_cam_alpha  / 0.0d0 /
      data nwpwxc_cam_beta   / 0.0d0 /
      data nwpwxc_cam_gamma  / 0.0d0 /
      data nwpwxc_rho_tol    / 1.0d-18 /
      data nwpwxc_disp_id    / NWXC_DISP_NONE /
      data nwpwxc_disp_s6    / 0.0d0 /
      data nwpwxc_disp_s8    / 0.0d0 /
      data nwpwxc_disp_sr6   / 0.0d0 /
      data nwpwxc_disp_sr8   / 0.0d0 /
      data nwpwxc_disp_alpha / 0.0d0 /
      end
C> @}
C> \ingroup nwpwxc_api
C> @{
C>
C> \brief Initialize the density functional library data structures
C>
C> Explicitly initialize the density functional library data structures.
C> This initialization is designed to match that of the `nwpwxc_block`
C> block data code. The only difference is that this routine can be
C> invoked explicitly at a moment of the application programmers choice.
C>
      subroutine nwpwxc_init
      implicit none
#include "nwpwxcP.fh"
#include "nwpwxc_param.fh"
      nwpwxc_depends_on_rho   = .false.
      nwpwxc_depends_on_gamma = .false.
      nwpwxc_depends_on_tau   = .false.
      nwpwxc_cam_srhfex       = .false.
      nwpwxc_on               = .false. 
      nwpwxc_num_f      = 0 
      nwpwxc_usr_num    = 0 
      nwpwxc_wght_hf    = 0.0d0 
      nwpwxc_wght_mp2   = 0.0d0
      nwpwxc_wght_cs00  = 0.0d0
      nwpwxc_cam_alpha  = 0.0d0
      nwpwxc_cam_beta   = 0.0d0
      nwpwxc_cam_gamma  = 0.0d0
      nwpwxc_rho_tol    = 1.0d-18
      nwpwxc_disp_id    = NWXC_DISP_NONE
      nwpwxc_disp_s6    = 0.0d0
      nwpwxc_disp_s8    = 0.0d0
      nwpwxc_disp_sr6   = 0.0d0
      nwpwxc_disp_sr8   = 0.0d0
      nwpwxc_disp_alpha = 0.0d0
      end
C>
C> \brief Reset the density tolerance
C>
C> Reset the tolerance on the electron density. In the functional
C> evaluation densities less than this tolerance will be ignored. 
C> Note that currently this routine resets only the tolerance stored
C> in memory. The updated value will not be reflected on the
C> runtime database.
C>
C> This routine was created to support code where the tolerance is
C> changed during the DFT SCF procedure (presumably to save work).
C>
      subroutine nwpwxc_reset_rho_tolerance(tol_rho)
      implicit none
#include "nwpwxcP.fh"
      double precision tol_rho
      nwpwxc_rho_tol = tol_rho
      end
C> @}
c $Id$
