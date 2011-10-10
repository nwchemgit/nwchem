      subroutine dft_roks_hessv_xx( basis, geom, nbf, 
     $     nmo, nclosed, nopen, 
     $     pflg,
     $     g_movecs, oskel, g_fcv, g_fpv, g_fcp,
     $     acc, lshift, g_x, g_ax )
C     $Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
c     
      integer basis, geom
      integer nbf, nmo, nclosed, nopen
      integer pflg
      integer g_movecs
      logical oskel
      integer g_fcv, g_fpv, g_fcp
      double precision acc
      double precision lshift
      integer g_x, g_ax, vlen, nvec, g_tmp, gtype
c
      call ga_zero(g_ax)
      if (pflg.gt.2 .or. pflg.le.0) then
        call errquit('dft_roks_hessv_xx: pflg invalid ', pflg,
     $               UNKNOWN_ERR)
      endif
      if (pflg .gt. 0) then
         call ga_dadd( lshift, g_x, 0.d0, g_ax, g_ax )
         call rohf_hessv_1e( basis, geom, nmo, nclosed, nopen,
     $        g_fcv, g_fpv, g_fcp,
     $        g_x, g_ax )
      endif
      if (pflg .gt. 1) then
         call dft_roks_hessv_2e( basis, geom, nbf, nmo, nclosed, nopen,
     $        g_movecs, oskel, g_x, acc, ! was min(1d-6,acc)
     $        g_ax )
      endif
c     
      end
