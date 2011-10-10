      subroutine dft_roks_hessv( acc, g_x, g_ax )
      implicit none
#include "errquit.fh"
#include "crohf.fh"
#include "cscf.fh"
c     
c     $Id$
c
      integer g_x, g_ax
      double precision acc
c     
      integer gtype,grow,gcol,growp,gcolp
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $  call errquit('dft_roks_hessv: ROKS internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x,gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $  call errquit('dft_roks_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_ax,gtype,growp,gcolp)
      if (growp.ne.crohf_vlen)
     $  call errquit('dft_roks_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      if (gcol.ne.gcolp)
     $  call errquit('dft_roks_hessv: invalid no. of vectors',0,
     &       UNKNOWN_ERR)
c     
c     Call internal routine
c     
      call dft_roks_hessv_xx( basis, geom, nbf, nmo,
     $     nclosed, nopen,
     $     pflg, g_movecs, oskel, 
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
     $     acc, lshift, g_x, g_ax )
c
c     Zap numbers much smaller than acc to ensure hard zeroes 
c     remain unpolluted ... cannot use a threshold larger than the
c     integral accuracy since can break symmetry in non-abelian groups
c     Also must ensure that the threshold tends to zero to permit
c     tight convergence.
c
      call ga_screen(g_ax, min(acc*acc,acc*0.01d0,1d-12))
c
      end
