      subroutine dft_uks_hessv(acc, g_x, g_ax)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "global.fh"

c     
      double precision acc      ! [input] required accuracy of products
      integer g_x               ! [input] handle to input vectors
      integer g_ax              ! [input] handle to output products
c
      integer gtype, vlen, nvec, nvecp, g_tmp
c
c     Multiply a set of vectors by the level-shifted UHF hessian.
c     
c     Check dimensions
c     
      if(.not.cuhf_init_flag)
     $     call errquit('dft_uks_hessv: UKS internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_ax,gtype,vlen,nvec)
      if (vlen.ne.cuhf_vlen)
     $     call errquit('dft_uks_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x,gtype,vlen,nvecp)
      if (vlen.ne.cuhf_vlen)
     $     call errquit('dft_uks_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      if (nvecp .ne. nvec) 
     $     call errquit('dft_uks_hessv: invalid nvecp',nvecp,
     &       UNKNOWN_ERR)
c
      call ga_zero(g_ax)
      if (pflg.gt.2 .or. pflg.le.0) then
        call errquit('dft_uks_hessv: pflg invalid ', pflg,
     $               UNKNOWN_ERR)
      endif
c
      if (pflg .gt. 0) then
         call ga_dadd( lshift, g_x, 0.d0, g_ax, g_ax )
         call uhf_hessv_1e(acc, g_x, g_ax, nvec)
      endif
c     
      if (pflg .gt. 1)then 
         call dft_uks_hessv_2e(acc, g_x, g_ax, nvec)
      endif
c
      end
