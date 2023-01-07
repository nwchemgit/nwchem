      subroutine int_dip_ga(ibas, jbas, g_x, g_y, g_z)
!
! $Id$
!
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
      integer ibas, jbas
      integer g_x, g_y, g_z
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mscratch, max1e
      integer ijshell, ilo, ihi, jlo, jhi
      integer l_buf, l_scr
      integer k_buf, k_scr
      logical do_mirr
      logical util_mirrmat
      external util_mirrmat
      integer g_x_mirr,g_y_mirr,g_z_mirr
      integer g_x_org,g_y_org,g_z_org
!
      double precision center(3)
!
      call dfill(3, 0.0d0, center, 1)
!
      call ga_sync()
      call ga_zero(g_x)
      call ga_zero(g_y)
      call ga_zero(g_z)
      g_x_mirr=0
      g_y_mirr=0
      g_z_mirr=0
      do_mirr=util_mirrmat(1,g_x,g_x_mirr,
     D     .false., .true.)
      do_mirr=do_mirr.and.util_mirrmat(1,g_y,g_y_mirr,
     D     .false., .true.)
      do_mirr=do_mirr.and.util_mirrmat(1,g_z,g_z_mirr,
     D     .false., .true.)
      if(do_mirr) then
         g_x_org=g_x
         g_x=g_x_mirr
         g_y_org=g_y
         g_y=g_y_mirr
         g_z_org=g_z
         g_z=g_z_mirr
            
      else
         if(g_x_mirr.ne.0) then
            if(.not.ga_destroy(g_x_mirr))call errquit(
     $           'could not destroy mirrx handle',1, GA_ERR)
         endif
         if(g_y_mirr.ne.0) then
            if(.not.ga_destroy(g_y_mirr))call errquit(
     $           'could not destroy mirry handle',1, GA_ERR)
         endif
      endif

      if(oscfps) call pstat_on(ps_int_dip)
!
!     grab basis set info type stuff
!
      if(.not.bas_numcont(ibas,nshell_i))call errquit
     $  ('int_dip_ga: bas_numcont failed for ibas',ibas, BASIS_ERR)
      if(.not.bas_numcont(jbas,nshell_j))call errquit
     $  ('int_dip_ga: bas_numcont failed for jbas',jbas, BASIS_ERR)
!
! allocate temp. arrays
!
      call int_mem_dipole(max1e,mscratch,ibas,jbas,1)
!
      if(.not.MA_push_get(MT_DBL,max1e,'int_dip_ga:buf',l_buf,k_buf))
     $  call errquit('int_dip_ga: ma failure',max1e, MA_ERR)
      if(.not.MA_push_get(MT_DBL,mscratch,'int_dip_ga:scr',
     $   l_scr,k_scr))
     $  call errquit('int_dip_ga: ma failure',mscratch, MA_ERR)

!
!     loop over shells
!
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      do jshell = 1,nshell_j
         do ishell = 1,nshell_i
            if(mod(ijshell,nproc).eq.iproc)then
               if(.not.bas_cn2bfr(ibas,ishell,ilo,ihi))
     $              call errquit('int_dip_ga:bas_cn2bfr ?', ibas,
     &       BASIS_ERR)
               if(.not.bas_cn2bfr(jbas,jshell,jlo,jhi))
     $              call errquit('int_dip_ga:bas_cn2bfr ?', jbas,
     &       BASIS_ERR)

!
!     get the integrals we want
!
!               call int_mpole(ibas, ishell, jbas, jshell,
               call int_mpole(jbas, jshell, ibas, ishell,
     $              1,
     $              center,
     $              mscratch, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
!
!     bung integrals into global array
!
               call int_dip_put_in_g(g_x, g_y, g_z, ilo,ihi,jlo,jhi,
     $              dbl_mb(k_buf+(ihi-ilo+1)*(jhi-jlo+1)))
            endif
            ijshell=ijshell+1
         enddo
      enddo
!
      if(.not.MA_pop_stack(l_scr))
     $   call errquit('int_dip_ga:pop failed',0, MA_ERR)
      if(.not.MA_pop_stack(l_buf))
     $   call errquit('int_dip_ga:pop failed',0, MA_ERR)
!
      if (do_mirr) then
         g_x=g_x_org
         call util_mirrmerge(g_x_mirr,g_x)
         if(.not.ga_destroy(g_x_mirr))call errquit(
     $        'could not destroy g_x_mirr handle',1, GA_ERR)
         g_y=g_y_org
         call util_mirrmerge(g_y_mirr,g_y)
         if(.not.ga_destroy(g_y_mirr))call errquit(
     $        'could not destroy g_y_mirr handle',1, GA_ERR)
         g_z=g_z_org
         call util_mirrmerge(g_z_mirr,g_z)
         if(.not.ga_destroy(g_z_mirr))call errquit(
     $        'could not destroy g_x_mirr handle',1, GA_ERR)
      endif
      call ga_sync()
!
      if (util_print('multipole', print_debug)) then
          call ga_print(g_x)
          call ga_print(g_y)
          call ga_print(g_z)
      endif
!
      if(oscfps) call pstat_off(ps_int_dip)
!
      end

      subroutine int_dip_put_in_g(g_x, g_y, g_z,
     1                            ilo, ihi, jlo, jhi, buf)
      implicit none
#include "global.fh"
      integer g_x, g_y, g_z, ilo, ihi, jlo, jhi
!      double precision buf(jlo:jhi,3,ilo:ihi)
      double precision buf(ilo:ihi,3,jlo:jhi)
      integer i, j
c      do i = ilo,ihi
         do j = jlo, jhi
            call ga_put(g_x,ilo,ihi,j,j,buf(ilo,1,j),1)
            call ga_put(g_y,ilo,ihi,j,j,buf(ilo,2,j),1)
            call ga_put(g_z,ilo,ihi,j,j,buf(ilo,3,j),1)
         end do
c      end do
      end

      subroutine int_qdr_ga(ibas, jbas, g_xx, g_xy, g_xz,
     1                                  g_yy, g_yz, g_zz)
!
! $Id$
!
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
      integer ibas, jbas
      integer g_xx, g_xy, g_xz, g_yy, g_yz, g_zz
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mscratch, max1e
      integer ijshell, ilo, ihi, jlo, jhi
      integer l_buf, l_scr
      integer k_buf, k_scr
!
      double precision center(3)
!
      call dfill(3, 0.0d0, center, 1)
!
      call ga_sync()

      call ga_zero(g_xx)
      call ga_zero(g_xy)
      call ga_zero(g_xz)
      call ga_zero(g_yy)
      call ga_zero(g_yz)
      call ga_zero(g_zz)

      if(oscfps) call pstat_on(ps_int_dip)
!
!     grab basis set info type stuff
!
      if(.not.bas_numcont(ibas,nshell_i))call errquit
     $  ('int_qdr_ga: bas_numcont failed for ibas',ibas, BASIS_ERR)
      if(.not.bas_numcont(jbas,nshell_j))call errquit
     $  ('int_qdr_ga: bas_numcont failed for jbas',jbas, BASIS_ERR)
!
! allocate temp. arrays
!
      call int_mem_dipole(max1e,mscratch,ibas,jbas,2)
!
      if(.not.MA_push_get(MT_DBL,max1e,'int_qdr_ga:buf',l_buf,k_buf))
     $  call errquit('int_qdr_ga: ma failure',max1e, MA_ERR)
      if(.not.MA_push_get(MT_DBL,mscratch,'int_qdr_ga:scr',
     $   l_scr,k_scr))
     $  call errquit('int_qdr_ga: ma failure',mscratch, MA_ERR)
!
!     loop over shells
!
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      do jshell = 1,nshell_j
         do ishell = 1,nshell_i
            if(mod(ijshell,nproc).eq.iproc)then
               if(.not.bas_cn2bfr(ibas,ishell,ilo,ihi))
     $              call errquit('int_qdr_ga:bas_cn2bfr ?', ibas,
     &       BASIS_ERR)
               if(.not.bas_cn2bfr(jbas,jshell,jlo,jhi))
     $              call errquit('int_qdr_ga:bas_cn2bfr ?', jbas,
     &       BASIS_ERR)
!
!     get the integrals we want
!
               call int_mpole(ibas, ishell, jbas, jshell,
     $              2,
     $              center,
     $              mscratch, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
!
!     bung integrals into global array
!
               call int_qdr_put_in_g(g_xx, g_xy, g_xz,
     $              g_yy, g_yz, g_zz, ilo,ihi,jlo,jhi,
     $              dbl_mb(k_buf+(ihi-ilo+1)*(jhi-jlo+1)*4))
            endif
            ijshell=ijshell+1
         enddo
      enddo
!
      if(.not.MA_pop_stack(l_scr))
     $   call errquit('int_qdr_ga:pop failed',0, MA_ERR)
      if(.not.MA_pop_stack(l_buf))
     $   call errquit('int_qdr_ga:pop failed',0, MA_ERR)
!
      call ga_sync()
!
      if (util_print('multipole', print_debug)) then
          call ga_print(g_xx)
          call ga_print(g_xy)
          call ga_print(g_xz)
          call ga_print(g_yy)
          call ga_print(g_yz)
          call ga_print(g_zz)
      endif
!
      if(oscfps) call pstat_off(ps_int_dip)
!
      end

      subroutine int_qdr_put_in_g(g_xx, g_xy, g_xz,
     1                            g_yy, g_yz, g_zz,
     2                            ilo, ihi, jlo, jhi, buf)
      implicit none
#include "global.fh"
      integer g_xx, g_xy, g_xz, g_yy, g_yz, g_zz
      integer ilo, ihi, jlo, jhi
      double precision buf(jlo:jhi,6,ilo:ihi)
      integer i, j
      do i = ilo,ihi
         do j = jlo, jhi
            call ga_put(g_xx,i,i,j,j,buf(j,1,i),1)
            call ga_put(g_xy,i,i,j,j,buf(j,2,i),1)
            call ga_put(g_xz,i,i,j,j,buf(j,3,i),1)
            call ga_put(g_yy,i,i,j,j,buf(j,4,i),1)
            call ga_put(g_yz,i,i,j,j,buf(j,5,i),1)
            call ga_put(g_zz,i,i,j,j,buf(j,6,i),1)
         end do
      end do
      end



      subroutine int_oct_ga(ibas, jbas,
     1                       g_xxx, g_xxy, g_xxz, g_xyy, g_xyz,
     2                       g_xzz, g_yyy, g_yyz, g_yzz, g_zzz)
!
! $Id$
!
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
      integer ibas, jbas
      integer g_xxx, g_xxy, g_xxz, g_xyy, g_xyz
      integer g_xzz, g_yyy, g_yyz, g_yzz, g_zzz
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mscratch, max1e
      integer ijshell, ilo, ihi, jlo, jhi
      integer l_buf, l_scr
      integer k_buf, k_scr
!
      double precision center(3)
!
      call dfill(3, 0.0d0, center, 1)
!
      call ga_sync()

      call ga_zero(g_xxx)
      call ga_zero(g_xxy)
      call ga_zero(g_xxz)
      call ga_zero(g_xyy)
      call ga_zero(g_xyz)
      call ga_zero(g_xzz)
      call ga_zero(g_yyy)
      call ga_zero(g_yyz)
      call ga_zero(g_yzz)
      call ga_zero(g_zzz)

      if(oscfps) call pstat_on(ps_int_dip)
!
!     grab basis set info type stuff
!
      if(.not.bas_numcont(ibas,nshell_i))call errquit
     $  ('int_qdr_ga: bas_numcont failed for ibas',ibas, BASIS_ERR)
      if(.not.bas_numcont(jbas,nshell_j))call errquit
     $  ('int_qdr_ga: bas_numcont failed for jbas',jbas, BASIS_ERR)
!
! allocate temp. arrays
!
      call int_mem_dipole(max1e,mscratch,ibas,jbas,3)
!
      if(.not.MA_push_get(MT_DBL,max1e,'int_oct_ga:buf',l_buf,k_buf))
     $  call errquit('int_oct_ga: ma failure',max1e, MA_ERR)
      if(.not.MA_push_get(MT_DBL,mscratch,'int_oct_ga:scr',
     $   l_scr,k_scr))
     $  call errquit('int_oct_ga: ma failure',mscratch, MA_ERR)
!
!     loop over shells
!
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      do jshell = 1,nshell_j
         do ishell = 1,nshell_i
            if(mod(ijshell,nproc).eq.iproc)then
               if(.not.bas_cn2bfr(ibas,ishell,ilo,ihi))
     $              call errquit('int_oct_ga:bas_cn2bfr ?', ibas,
     &       BASIS_ERR)
               if(.not.bas_cn2bfr(jbas,jshell,jlo,jhi))
     $              call errquit('int_oct_ga:bas_cn2bfr ?', jbas,
     &       BASIS_ERR)
!
!     get the integrals we want
!
               call int_mpole(ibas, ishell, jbas, jshell,
     $              3,
     $              center,
     $              mscratch, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
!
!     bung integrals into global array
!
               call int_oct_put_in_g(g_xxx, g_xxy, g_xxz,
     $              g_xyy, g_xyz, g_xzz, g_yyy, g_yyz,
     $              g_yzz, g_zzz, ilo, ihi, jlo, jhi,
     $              dbl_mb(k_buf+(ihi-ilo+1)*(jhi-jlo+1)*10))
            endif
            ijshell=ijshell+1
         enddo
      enddo
!
      if(.not.MA_pop_stack(l_scr))
     $   call errquit('int_oct_ga:pop failed',0, MA_ERR)
      if(.not.MA_pop_stack(l_buf))
     $   call errquit('int_oct_ga:pop failed',0, MA_ERR)
!
      call ga_sync()
!
      if (util_print('multipole', print_debug)) then
           call ga_print(g_xxx)
           call ga_print(g_xxy)
           call ga_print(g_xxz)
           call ga_print(g_xyy)
           call ga_print(g_xyz)
           call ga_print(g_xzz)
           call ga_print(g_yyy)
           call ga_print(g_yyz)
           call ga_print(g_yzz)
           call ga_print(g_zzz)
      endif
!
      if(oscfps) call pstat_off(ps_int_dip)
!
      end

      subroutine int_oct_put_in_g(
     1              g_xxx, g_xxy, g_xxz, g_xyy, g_xyz,
     2              g_xzz, g_yyy, g_yyz, g_yzz, g_zzz,
     3              ilo, ihi, jlo, jhi, buf)
      implicit none
#include "global.fh"
      integer g_xxx, g_xxy, g_xxz, g_xyy, g_xyz
      integer g_xzz, g_yyy, g_yyz, g_yzz, g_zzz
      integer ilo, ihi, jlo, jhi
      double precision buf(jlo:jhi,10,ilo:ihi)
      integer i, j
      do i = ilo,ihi
         do j = jlo, jhi
            call ga_put(g_xxx,i,i,j,j,buf(j, 1,i),1)
            call ga_put(g_xxy,i,i,j,j,buf(j, 2,i),1)
            call ga_put(g_xxz,i,i,j,j,buf(j, 3,i),1)
            call ga_put(g_xyy,i,i,j,j,buf(j, 4,i),1)
            call ga_put(g_xyz,i,i,j,j,buf(j, 5,i),1)
            call ga_put(g_xzz,i,i,j,j,buf(j, 6,i),1)
            call ga_put(g_yyy,i,i,j,j,buf(j, 7,i),1)
            call ga_put(g_yyz,i,i,j,j,buf(j, 8,i),1)
            call ga_put(g_yzz,i,i,j,j,buf(j, 9,i),1)
            call ga_put(g_zzz,i,i,j,j,buf(j,10,i),1)
         end do
      end do
      end



