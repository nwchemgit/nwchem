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
      integer ilo, ihi, jlo, jhi
      integer l_buf, l_scr
      integer k_buf, k_scr
      logical do_mirr
      logical util_mirrmat
      external util_mirrmat
      integer g(3)
      integer i0,i1,j0,j1
      integer iptr0,iptr1,jptr0,jptr1
!
      double precision center(3)
!
      call dfill(3, 0.0d0, center, 1)
!
      g(1)=g_x
      g(2)=g_y
      g(3)=g_z
      call ga_sync()
      call ga_zero(g_x)
      call ga_zero(g_y)
      call ga_zero(g_z)

      if(oscfps) call pstat_on(ps_int_dip)
!
!     grab basis set info type stuff
!
      call ga_distribution(g_x,ga_nodeid(),i0,i1,j0,j1)
cdbg       write(6,'(i5,a,5i5)')
cdbg     I ga_nodeid(),'DDISTRIB',g_x,i0,i1,j0,j1
      if (i0.gt.0 .and. i0.le.i1) then
         if(.not.bas_numcont(ibas,nshell_i))call errquit
     $        ('int_dip_ga: bas_numcont failed for ibas',ibas,BASIS_ERR)
         if(.not.bas_numcont(jbas,nshell_j))call errquit
     $        ('int_dip_ga: bas_numcont failed for jbas',jbas,BASIS_ERR)
!
! allocate temp. arrays
!
         call int_mem_dipole(max1e,mscratch,ibas,jbas,1)
!
         if(.not.MA_push_get(MT_DBL,max1e,'int_dip_ga:buf',l_buf,k_buf))
     $        call errquit('int_dip_ga: ma failure',max1e, MA_ERR)
         if(.not.MA_push_get(MT_DBL,mscratch,'int_dip_ga:scr',
     $        l_scr,k_scr))
     $        call errquit('int_dip_ga: ma failure',mscratch, MA_ERR)

!
!     loop over shells
!
         iproc = ga_nodeid()
         do jshell = 1,nshell_j
            do ishell = 1,nshell_i
               if(.not.bas_cn2bfr(ibas,ishell,ilo,ihi))
     $              call errquit('int_dip_ga:bas_cn2bfr ?', ibas,
     &              BASIS_ERR)
               if(.not.bas_cn2bfr(jbas,jshell,jlo,jhi))
     $              call errquit('int_dip_ga:bas_cn2bfr ?', jbas,
     &              BASIS_ERR)
               iptr0=max(i0,ilo)
               iptr1=min(i1,ihi)
               jptr0=max(j0,jlo)
               jptr1=min(j1,jhi)
               if(iptr1.ge.iptr0.and.jptr1.ge.jptr0) then
!     get the integrals we want
!
                  call int_mpole(jbas, jshell, ibas, ishell,
     $                 1,
     $                 center,
     $                 mscratch, dbl_mb(k_scr), max1e,
     D                 dbl_mb(k_buf))
!     
!     bung integrals into global array
!
                  call int_mpole_put_in_g(g, 3,
     L                 ilo,ihi,jlo,jhi,
     P                 iptr0,iptr1,jptr0,jptr1,
     $                 dbl_mb(k_buf+(ihi-ilo+1)*(jhi-jlo+1)))
               endif
            enddo
         enddo
!
         if(.not.MA_pop_stack(l_scr))
     $        call errquit('int_dip_ga:pop failed',0, MA_ERR)
         if(.not.MA_pop_stack(l_buf))
     $   call errquit('int_dip_ga:pop failed',0, MA_ERR)

      endif 
      call ga_sync()

c      if (.true.) then
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
      subroutine int_mpole_put_in_g(g,
     M     mult_size,
     L     ilo,ihi,jlo,jhi,
     1     iptr0, iptr1, jptr0, jptr1,
     B     buf)
      implicit none
#include "global.fh"
      integer mult_size,mult
      integer g(mult_size)
      integer ilo,ihi,jlo,jhi
      integer iptr0, iptr1, jptr0, jptr1
      double precision buf(ilo:ihi,mult_size,jlo:jhi)
      integer i, j
      do j = jptr0, jptr1
         do mult=1,mult_size
            call ga_put(g(mult),iptr0,iptr1,j,j,buf(iptr0,mult,j),1)
         end do
      end do
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
      integer g(6)
      integer i0,i1,j0,j1
      integer iptr0,iptr1,jptr0,jptr1
      
!
      double precision center(3)
!
      call dfill(3, 0.0d0, center, 1)
!
      call ga_sync()

      g(1)=g_xx
      g(2)=g_xy
      g(3)=g_xz
      g(4)=g_yy
      g(5)=g_yz
      g(6)=g_zz
      call ga_zero(g_xx)
      call ga_zero(g_xy)
      call ga_zero(g_xz)
      call ga_zero(g_yy)
      call ga_zero(g_yz)
      call ga_zero(g_zz)
      call ga_distribution(g_xx,ga_nodeid(),i0,i1,j0,j1)
c      if(ga_nodeid().eq.0) write(6,'(a,4i5)')
cdbg       write(6,'(i5,a,5i5)')
cdbg     I ga_nodeid(),'QDISTRIB',g_xx,i0,i1,j0,j1
      if(oscfps) call pstat_on(ps_int_dip)
!
!     grab basis set info type stuff
!
      if (i0.gt.0 .and. i0.le.i1) then
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
      ijshell = 0
      do jshell = 1,nshell_j
         do ishell = 1,nshell_i
!            if(mod(ijshell,nproc).eq.iproc)then
               if(.not.bas_cn2bfr(ibas,ishell,ilo,ihi))
     $              call errquit('int_qdr_ga:bas_cn2bfr ?', ibas,
     &       BASIS_ERR)
               if(.not.bas_cn2bfr(jbas,jshell,jlo,jhi))
     $              call errquit('int_qdr_ga:bas_cn2bfr ?', jbas,
     &       BASIS_ERR)
               iptr0=max(i0,ilo)
               iptr1=min(i1,ihi)
               jptr0=max(j0,jlo)
               jptr1=min(j1,jhi)
               if(iptr1.ge.iptr0.and.jptr1.ge.jptr0) then
!
!     get the integrals we want
!
               call int_mpole(jbas, ishell, ibas, jshell,
     $              2,
     $              center,
     $              mscratch, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
!
!     bung integrals into global array
!
               call int_mpole_put_in_g(g,6,
     $              ilo,ihi,jlo,jhi,
     P              iptr0,iptr1,jptr0,jptr1,
     $              dbl_mb(k_buf+(ihi-ilo+1)*(jhi-jlo+1)*4))
            endif
            ijshell=ijshell+1
         enddo
      enddo
!
      if(.not.MA_pop_stack(l_scr))
     $   call errquit('int_qdr_ga:pop failed',0, MA_ERR)
      if(.not.MA_pop_stack(l_buf))
     $     call errquit('int_qdr_ga:pop failed',0, MA_ERR)
      endif
!
      call ga_sync()
!
      if (util_print('multipole', print_debug)) then
c      if (.true.)  then
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
#if 0
         do j = jlo, jhi
            call ga_put(g_xx,i,i,j,j,buf(j,1,i),1)
            call ga_put(g_xy,i,i,j,j,buf(j,2,i),1)
            call ga_put(g_xz,i,i,j,j,buf(j,3,i),1)
            call ga_put(g_yy,i,i,j,j,buf(j,4,i),1)
            call ga_put(g_yz,i,i,j,j,buf(j,5,i),1)
            call ga_put(g_zz,i,i,j,j,buf(j,6,i),1)
         end do
#else
            call ga_put(g_xx,jlo,jhi,i,i,buf(jlo,1,i),1)
            call ga_put(g_xy,jlo,jhi,i,i,buf(jlo,2,i),1)
            call ga_put(g_xz,jlo,jhi,i,i,buf(jlo,3,i),1)
            call ga_put(g_yy,jlo,jhi,i,i,buf(jlo,4,i),1)
            call ga_put(g_yz,jlo,jhi,i,i,buf(jlo,5,i),1)
            call ga_put(g_zz,jlo,jhi,i,i,buf(jlo,6,i),1)
#endif
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



