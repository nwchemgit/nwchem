      subroutine rohf_hessv_xx( basis, geom, nbf, nmo, nclosed, nopen, 
     $     pflg,
     $     g_movecs, oskel, g_fcv, g_fpv, g_fcp,
     $     acc, lshift, g_x, g_ax )
C     $Id: rohf_hxx.F,v 1.2 1996-01-19 07:29:14 d3g681 Exp $
      implicit none
#include "global.fh"
#include "mafdecls.fh"
c     
      integer basis, geom
      integer nbf, nmo, nclosed, nopen
      integer pflg
      integer g_movecs
      logical oskel
      integer g_fcv, g_fpv, g_fcp
      double precision acc
      double precision lshift
      integer g_x, g_ax
c     
      call ga_zero(g_ax)
c
      if (pflg.gt.2 .or. pflg.le.0) then
         call errquit('rohf_hessv_xx: pflg invalid ', pflg)
      endif

      if (pflg .gt. 0) then
         call ga_dadd( lshift, g_x, 0.d0, g_ax, g_ax )
         call rohf_hessv_1e( basis, geom, nmo, nclosed, nopen,
     $        g_fcv, g_fpv, g_fcp,
     $        g_x, g_ax )
      endif

      if (pflg .gt. 1) then
         call rohf_hessv_2e( basis, geom, nbf, nmo, nclosed, nopen,
     $        g_movecs, oskel, g_x, min(1.0d-6,acc),
     $        g_ax )
      endif

c     
      end



