      logical function task_check_sum(rtdb)
C$Id: do_chk_sum.F,v 1.4 1996-03-11 18:17:36 d3e129 Exp $
      implicit none
      integer rtdb
c
      call do_int_chk_sum(rtdb,'geometry','ao basis',1.0d-10)
      task_check_sum = .true.
      end
      subroutine do_int_chk_sum(rtdb,geom_name, basis_name, tol2e)
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c
      integer rtdb
      integer basis, geom
      double precision cpu, wall
      double precision tol2e
      character*(*) geom_name, basis_name
c
      logical int_normalize
*--      logical int_norm_2c
      external int_normalize
*--      external int_norm_2c
c
      write(luout,*)' do_int_chk_sum:rtdb      :',rtdb
      write(luout,*)' do_int_chk_sum:geom_name : <',geom_name,'>'
      write(luout,*)' do_int_chk_sum:basis_name: <',basis_name,'>'
c
      if (.not.rtdb_parallel(.true.))
     &     call errquit('do_int_chk_sum: rtdb_parallel failed?', 0)
      if (.not. geom_create(geom, geom_name))
     $     call errquit('do_int_chk_sum: geom_create failed?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
     $     call errquit('do_int_chk_sum: geom_load failed', 0)
      if (.not. bas_create(basis, basis_name))
     $     call errquit('do_int_chk_sum: basis create failed', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
     $      call errquit('do_int_chk_sum: basis load failed', 0)

      if (ga_nodeid().eq.0) then
        if (.not. geom_print(geom))
     $      call errquit('do_int_chk_sum: geom_print failed', 0)
        if (.not. bas_print(basis))
     $      call errquit('do_int_chk_sum: basis print failed', 0)
        if (.not. gbs_map_print(basis))
     $      call errquit('do_int_chk_sum: gbs_map_print failed', 0)
      endif
*---------------------------------- 
#define NORMFIRST
#if defined(NORMFIRST)
c
c normalize basis set
c
      if (.not.int_normalize(basis))
     &       call errquit('do_int_chk_sum: basis norm. failed', 0)
*--      if (.not.int_norm_2c(basis))
*--     &       call errquit('do_int_chk_sum: basis norm. failed', 0)
*--      if (.not. bas_print(basis))
*--     $    call errquit('do_int_chk_sum: basis print failed', 0)
c
      call int_init(rtdb,1,basis)
      call int_acc_set(tol2e)
#else
      call int_init(rtdb,1,basis)
      call int_acc_set(tol2e)
c
c normalize basis set
c
      if (.not.int_normalize(basis))
     &       call errquit('do_int_chk_sum: basis norm. failed', 0)
*--      if (.not.int_norm_2c(basis))
*--     &       call errquit('do_int_chk_sum: basis norm. failed', 0)
*--      if (.not. bas_print(basis))
*--     $    call errquit('do_int_chk_sum: basis print failed', 0)
c
#endif
      cpu  = util_cpusec()
      wall = util_wallsec()
      call intb_chk_sum(rtdb,basis,.false.)
      call int_chk_sum(rtdb,basis,.false.)
*--      call int_chk_sum(rtdb,basis,.true.)
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      write(luout,'(1x,a,f10.2)')
     &    'checksum cpu  time:',cpu
      write(luout,'(1x,a,f10.2)')
     &    'checksum wall time:',wall
c
      call int_terminate()
      if (.not.(
     &       (bas_destroy(basis)).and.(geom_destroy(geom))
     &       ))
     &       call errquit('rak:error destroying gbs pair',0)
      if (ga_nodeid().eq.0) call MA_summarize_allocated_blocks()
c
      end
