      subroutine rohf_hessv_xx( basis, geom, nbf, nmo, nclosed, nopen, 
     $     pflg,
     $     g_movecs, oskel, g_fcv, g_fpv, g_fcp,
     $     acc, lshift, g_x, g_ax )
C     $Id: rohf_hxx.F,v 1.4 2001-06-21 18:58:07 zzhang Exp $
      implicit none
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
      logical cphf_tddft 
      cphf_tddft = .false.
c
      if (.not. rtdb_get(bgj_get_rtdb_handle(), 
     &     'cphf_solve:cphf_tddft', mt_log, 1, cphf_tddft)) then
         cphf_tddft = .false.
      endif      
      call ga_zero(g_ax)
      if (pflg.gt.2 .or. pflg.le.0) then
         call errquit('rohf_hessv_xx: pflg invalid ', pflg)
      endif
      if(cphf_tddft)then
         call ga_inquire(g_ax,gtype,vlen,nvec)
         if (vlen.ne.((nmo-nclosed)*nclosed))
     $        call errquit('rohf_hxx: invalid vector length',0)
         if (.not. ga_create(MT_DBL, vlen, nvec, 'tddft_precond: tmp',
     $        1, 0, g_tmp)) call errquit('diag', 0)
         call ga_copy(g_x, g_tmp)
         call tddft_scale(g_tmp, g_tmp, .false.)
         goto 200               ! skip 1e part
      endif
      if (pflg .gt. 0) then
         call ga_dadd( lshift, g_x, 0.d0, g_ax, g_ax )
         call rohf_hessv_1e( basis, geom, nmo, nclosed, nopen,
     $        g_fcv, g_fpv, g_fcp,
     $        g_x, g_ax )
      endif
 200  continue
      if (pflg .gt. 1) then
         if(cphf_tddft)then
            call rohf_hessv_2e( basis, geom, nbf, nmo, nclosed, nopen,
     $           g_movecs, oskel, g_tmp, acc, ! was min(1d-6,acc)
     $           g_ax ) 
         else
            call rohf_hessv_2e( basis, geom, nbf, nmo, nclosed, nopen,
     $           g_movecs, oskel, g_x, acc, ! was min(1d-6,acc)
     $           g_ax )
         endif
      endif
      if(cphf_tddft)then
         call tddft_scale(g_tmp, g_ax, .true.)
         if(.not.ga_destroy(g_tmp))
     &        call errquit("uhf_hessv: destroy g_tmp?", 0)
      endif
c     
      end



