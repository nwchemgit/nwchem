C> \ingroup nwxc
C> @{
C>
C> \brief One time initialization of the density functional library
C> data structures
C>
C> Initialize the density functional library data structures ensuring
C> that initially there are no functional terms and all weights are
C> zero.
C>
      block data nwxc_block
      implicit none
#include "nwxcP.fh"
      data nwxc_depends_on_rho   / .false. /
      data nwxc_depends_on_gamma / .false. /
      data nwxc_depends_on_tau   / .false. /
      data nwxc_num_f     / 0 /
      data nwxc_usr_num   / 0 /
      data nwxc_wght_hf   / 0.0d0 /
      data nwxc_wght_mp2  / 0.0d0 /
      data nwxc_cam_alpha / 0.0d0 /
      data nwxc_cam_beta  / 0.0d0 /
      data nwxc_cam_gamma / 0.0d0 /
      end
C>
C> \brief Initialize the density functional library data structures
C>
C> Explicitly initialize the density functional library data structures.
C> This initialization is designed to match that of the `nwxc_block`
C> block data code. The only difference is that this routine can be
C> invoked explicitly at a moment of the application programmers choice.
C>
      subroutine nwxc_init
      implicit none
#include "nwxcP.fh"
      nwxc_depends_on_rho   = .false.
      nwxc_depends_on_gamma = .false.
      nwxc_depends_on_tau   = .false.
      nwxc_num_f     = 0 
      nwxc_usr_num   = 0 
      nwxc_wght_hf   = 0.0d0 
      nwxc_wght_mp2  = 0.0d0
      nwxc_cam_alpha = 0.0d0
      nwxc_cam_beta  = 0.0d0
      nwxc_cam_gamma = 0.0d0
      end
C> @}
