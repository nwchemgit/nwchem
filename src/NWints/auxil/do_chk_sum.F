C> \ingroup task
C> @{
      logical function task_check_sum(rtdb)
C$Id$
      implicit none
      integer rtdb
c
      call do_int_chk_sum(rtdb,'geometry','ao basis',
     &    'ri basis',1.0d-10)
      task_check_sum = .true.
      end
C> @}
      subroutine do_int_chk_sum(rtdb,geom_name, basis_name,
     &    fit_basis_name, tol2e)
      implicit none
#include "geom.fh"
#include "errquit.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c
      integer rtdb
      integer basis, geom, ribasis
      integer bases(2)
      double precision cpu, wall
      double precision tol2e
      character*(*) geom_name, basis_name, fit_basis_name
c
      logical doitcs_all
      logical doitcs_1e
      logical doitcs_ov
      logical doitcs_ke
      logical doitcs_pe
      logical doitcs_h1
      logical doitcs_3ov
      logical doitcs_2e
      logical doitcs_b2e4c
      logical doitcs_2e4c
      logical doitcs_2e3c
      logical doitcs_2e2c

      logical doitcs_any

      logical ribasis_exist
      logical dummyL
c
      logical int_normalize
*--      logical int_norm_2c
      external int_normalize
*--      external int_norm_2c
c
      write(luout,*)' do_int_chk_sum:rtdb          :',rtdb
      write(luout,*)' do_int_chk_sum:geom_name     : <',geom_name,'>'
      write(luout,*)' do_int_chk_sum:basis_name    : <',basis_name,'>'
      write(luout,*)' do_int_chk_sum:fit_basis_name: <',
     &    fit_basis_name,'>'
c
      if (.not.rtdb_parallel(.true.))
     &     call errquit('do_int_chk_sum: rtdb_parallel failed?', 0,
     &             RTDB_ERR)
      if (.not. geom_create(geom, geom_name))
     &     call errquit('do_int_chk_sum: geom_create failed?', 0,
     &                GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
     &     call errquit('do_int_chk_sum: geom_load failed', 0,
     &                GEOM_ERR)
      if (.not. bas_create(basis, basis_name))
     &     call errquit('do_int_chk_sum: basis create failed', 0,
     &            BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
     &      call errquit('do_int_chk_sum: basis load failed', 0,
     &            BASIS_ERR)
      if (.not. bas_create(ribasis, fit_basis_name))
     &     call errquit('do_int_chk_sum: basis create failed', 0,
     &              BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb, geom, ribasis, fit_basis_name)) then
        if (.not.bas_destroy(ribasis)) call errquit
     &      ('do_int_chk_sum:failed to destroy ribasis handle',911,
     &           INT_ERR)
        ribasis = -1 
        ribasis_exist = .false.
      else
        ribasis_exist = .true.
      endif
c
      if (ga_nodeid().eq.0) then
        if (.not. geom_print(geom))
     &      call errquit('do_int_chk_sum: geom_print failed', 0,
     &             GEOM_ERR)
        if (.not. bas_print(basis))
     &      call errquit('do_int_chk_sum: basis print failed', 0,
     &             BASIS_ERR)
        if (.not. gbs_map_print(basis))
     &      call errquit('do_int_chk_sum: gbs_map_print failed', 0,
     &             BASIS_ERR)
        if (ribasis_exist) then
          if (.not. bas_print(ribasis))
     &        call errquit('do_int_chk_sum: fit basis print failed', 0,
     &             BASIS_ERR)
          if (.not. gbs_map_print(ribasis))
     &        call errquit('do_int_chk_sum: gbs_map_print failed', 0,
     &              BASIS_ERR)
        endif
      endif
*---------------------------------- 
#define NORMFIRST
#if defined(NORMFIRST)
c
c normalize basis set
c
      if (.not.int_normalize(rtdb,basis))
     &       call errquit('do_int_chk_sum: basis norm. failed', 0,
     &             BASIS_ERR)
*--      if (.not.int_norm_2c(rtdb,basis))
*--     &       call errquit('do_int_chk_sum: basis norm. failed', 0)
*--      if (.not. bas_print(basis))
*--     &    call errquit('do_int_chk_sum: basis print failed', 0)
c
      bases(1) = basis
      bases(2) = ribasis
c
      if (ribasis_exist) then
        call int_init(rtdb,2,bases)
      else
        call int_init(rtdb,1,bases)
      endif
      call int_acc_set(tol2e)
#else
      bases(1) = basis
      bases(2) = ribasis
c
      if (ribasis_exist) then
        call int_init(rtdb,2,bases)
      else
        call int_init(rtdb,1,bases)
      endif
      call int_acc_set(tol2e)
c
c normalize basis set
c
      if (.not.int_normalize(rtdb,basis))
     &       call errquit('do_int_chk_sum: basis norm. failed', 0)
*--      if (.not.int_norm_2c(rtdb,basis))
*--     &       call errquit('do_int_chk_sum: basis norm. failed', 0)
*--      if (.not. bas_print(basis))
*--     &    call errquit('do_int_chk_sum: basis print failed', 0)
c
#endif
*
* check flags on rtdb to see what is to be done.
* *default is to do nothing
* computational flags:
* intcsum:all
* intcsum:1e 
* intcsum:ov
* intcsum:ke
* intcsum:pe
* intcsum:h1
* intcsum:3ov
* intcsum:2e 
* intcsum:b2e4c
* intcsum:2e4c
* intcsum:2e3c
* intcsum:2e2c
*
* print flags
* intcsum:ovprint
* intcsum:keprint
* intcsum:peprint
* intcsum:h1print
* intcsum:3ovprint
* intcsum:b2e4cprint
* intcsum:2e4cprint
* intcsum:2e3cprint
* intcsum:2e2cprint
*
c determine computational flag setup
      doitcs_all   = .false.
      doitcs_1e    = .false.
      doitcs_ov    = .false.
      doitcs_ke    = .false.
      doitcs_pe    = .false.
      doitcs_h1    = .false.
      doitcs_3ov   = .false.
      doitcs_2e    = .false.
      doitcs_b2e4c = .false.
      doitcs_2e4c  = .false.
      doitcs_2e3c  = .false.
      doitcs_2e2c  = .false.
c
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:all',MT_LOG,1,dummyL)) then
        doitcs_all = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:1e',MT_LOG,1,dummyL)) then
        doitcs_1e = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:2e',MT_LOG,1,dummyL)) then
        doitcs_2e = dummyL
      endif
c
      if (doitcs_all) then
        doitcs_ov    = .true.
        doitcs_ke    = .true.
        doitcs_pe    = .true.
        doitcs_h1    = .true.
        doitcs_3ov   = .true.
        doitcs_b2e4c = .true.
        doitcs_2e4c  = .true.
        doitcs_2e3c  = .true.
        doitcs_2e2c  = .true.
      endif
      if (doitcs_1e) then
        doitcs_ov    = .true.
        doitcs_ke    = .true.
        doitcs_pe    = .true.
        doitcs_3ov   = .true.
      endif
      if (doitcs_2e) then
        doitcs_b2e4c = .true.
        doitcs_2e4c  = .true.
        doitcs_2e3c  = .true.
        doitcs_2e2c  = .true.
      endif
c      
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:ov',MT_LOG,1,dummyL)) then
        doitcs_ov = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:ke',MT_LOG,1,dummyL)) then
        doitcs_ke = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:pe',MT_LOG,1,dummyL)) then
        doitcs_pe = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:h1',MT_LOG,1,dummyL)) then
        doitcs_h1 = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:3ov',MT_LOG,1,dummyL)) then
        doitcs_3ov = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:b2e4c',MT_LOG,1,dummyL)) then
        doitcs_b2e4c = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:2e4c',MT_LOG,1,dummyL)) then
        doitcs_2e4c = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:2e3c',MT_LOG,1,dummyL)) then
        doitcs_2e3c = dummyL
      endif
      dummyL = .false.
      if (rtdb_get(rtdb,'intcsum:2e2c',MT_LOG,1,dummyL)) then
        doitcs_2e2c = dummyL
      endif
c
      doitcs_any = doitcs_ov.or.doitcs_ke.or.doitcs_pe.or.doitcs_h1
      doitcs_any = doitcs_any.or.doitcs_3ov.or.doitcs_b2e4c
      doitcs_any = doitcs_any.or.doitcs_2e4c.or.doitcs_2e3c
      doitcs_any = doitcs_any.or.doitcs_2e2c
      if (.not.doitcs_any) then
        write(luout,*)' no specified tasks for checksum'
        write(luout,*)' add one of the following set directives ',
     &      'to your input deck'
        write(luout,*)' '
        write(luout,*)' computational flags:'
        write(luout,*)' set intcsum:all logical true'
        write(luout,*)' set intcsum:1e logical true'
        write(luout,*)' set intcsum:2e logical true'
        write(luout,*)' set intcsum:ov logical true'
        write(luout,*)' set intcsum:ke logical true'
        write(luout,*)' set intcsum:pe logical true'
        write(luout,*)' set intcsum:h1 logical true'
        write(luout,*)' set intcsum:3ov logical true'
        write(luout,*)' set intcsum:b2e4c logical true'
        write(luout,*)' set intcsum:2e4c logical true'
        write(luout,*)' set intcsum:2e3c logical true'
        write(luout,*)' set intcsum:2e2c logical true'
        write(luout,*)' '
        write(luout,*)' '
        write(luout,*)' print flags:'
        write(luout,*)' set intcsum:ovprint logical true'
        write(luout,*)' set intcsum:keprint logical true'
        write(luout,*)' set intcsum:peprint logical true'
        write(luout,*)' set intcsum:h1print logical true'
        write(luout,*)' set intcsum:3ovprint logical true'
        write(luout,*)' set intcsum:b2e4cprint logical true'
        write(luout,*)' set intcsum:2e4cprint logical true'
        write(luout,*)' set intcsum:2e3cprint logical true'
        write(luout,*)' set intcsum:2e2cprint logical true'
        write(luout,*)' '
        write(luout,*)' '
      endif
c
      cpu  = util_cpusec()
      wall = util_wallsec()
      
      if (doitcs_ov) then
        call int_chk_sum_ov(rtdb,basis,.false.)
      endif
      if (doitcs_ke) then
        call int_chk_sum_ke(rtdb,basis,.false.)
      endif
      if (doitcs_pe) then
        call int_chk_sum_pe(rtdb,basis,.false.)
      endif
      if (doitcs_h1) then
        call int_chk_sum_h1(rtdb,basis,.false.)
      endif
      if (doitcs_3ov) then
        call int_chk_sum_3ov(rtdb,basis,ribasis,.false.)
      endif
      if (doitcs_2e2c) then
        call int_chk_sum_2e2c(rtdb,basis,ribasis,.false.)
      endif
      if (doitcs_2e3c) then
        call int_chk_sum_2e3c(rtdb,basis,ribasis,.false.)
      endif
      if (doitcs_2e4c) then
        call int_chk_sum_2e4c(rtdb,basis,.false.)
      endif
      if (doitcs_b2e4c) then
        call intb_chk_sum(rtdb,basis,.false.)
      endif
*      call int_chk_sum(rtdb,basis,.false.)
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
c
      if (ribasis_exist) then
        if (.not.(bas_destroy(ribasis)))
     &      call errquit('rak:error destroying ribasis ',0,
     &              BASIS_ERR)
      endif
      if (.not.(bas_destroy(basis)))
     &    call errquit('rak:error destroying basis',0,
     &            BASIS_ERR)
      if (.not.(geom_destroy(geom)))
     &    call errquit('rak:error destroying geometry',0, BASIS_ERR)
      if (ga_nodeid().eq.0) call MA_summarize_allocated_blocks()
c
      end
