*
* $Id: task_vscf.F,v 1.1 2000-02-25 00:40:53 windus Exp $
*
      logical function task_vscf(rtdb)
c     
c     Also will be hooked up to analytic methods as they are available.
c
c     Since this routine is directly invoked by application modules
c     no input is processed in this routine.
c
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
#include "global.fh"
#include "eaf.fh"
#include "util.fh"
c
      logical task_freq
      external task_freq
c     
      integer rtdb
c     
      integer nat, geom
      integer ngrid, ncoup, nresmode(4), iexc, nread
      integer l_freq, k_freq, l_nmode, k_nmode, nc
      integer offset, fd
      character*255 filename
      double precision vcfct
      logical status, dmdr, readv
c     
      double precision cpu, wall
c 
      logical      ignore
c
c  set up a couple of variables that I will eventually let be user controlled
c
      dmdr = .false.
      readv = .false. ! temporarily don't let anyone read in previous results
      nread = 0
c    
C      call ecce_print_module_entry('task frequencies')
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_vscf: failed to invalidate status',0)
      if (ga_nodeid().eq.0 .and.
     $    util_print('task_vscf', print_low)) then
        write(LUOut,*)
        write(LUOut,*)
        call util_print_centered(LUOut,
     $      'VSCF Analysis',
     $      40,.true.)
        write(LUOut,*)
      endif

      status = task_freq(rtdb)
      if (.not.status) call errquit
     &    ('task_vscf: task_freq failed',911)
c     ignore = rtdb_parallel(.false.)
c
c     create/load reference geometry to get the number of atoms
c
      if (.not.geom_create(geom,'geometry')) call errquit
     $      ('task_vscf:geom_create failed?',1)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     $      call errquit
     $      ('task_vscf:geom_rtdb_load failed?',2)
      if (.not. geom_ncent(geom,nat)) call errquit
     $      ('task_vscf:geom_ncent failed?',3)
      nc = nat*3
      if (.not. geom_destroy(geom)) call errquit
     $      ('task_vscf:geom_destroy failed?',911)
c
c     load parameters from rtdb or assign defaults
c
      if (.not. rtdb_get(rtdb,'vscf:ngrid',mt_int,1,ngrid))
     $      ngrid = 16
      if (.not. rtdb_get(rtdb,'vscf:iexc',mt_int,1,iexc))
     $      iexc = 1
      if (.not. rtdb_get(rtdb,'vscf:ncoup',mt_int,1,ncoup))
     $      ncoup = 2
      if (.not. rtdb_get(rtdb,'vscf:vcfct',mt_dbl,1,vcfct))
     $      vcfct = 1.0d+0
      if(.not.rtdb_put(rtdb,'vscf:nresmode',mt_int,4,nresmode))then
         nresmode(1) = 1
         nresmode(2) = 1
         nresmode(3) = 1
         nresmode(4) = 1
      endif
c
c  read in frequencies and normal modes
c
      if (ga_nodeid().eq.0) then
        call util_file_name('nmode',.false.,.false.,filename)
        if (eaf_open(filename,eaf_rw,fd).ne.0)
     &    call errquit('task_vscf:eaf_open failed',555)
      if (.not. ma_push_get(mt_dbl, nc,'frequencies', l_freq, k_freq))
     &  call errquit('task_vscf: could not allocate l_freq',nc)                 
      if (.not. 
     &  ma_push_get(mt_dbl, nc*nc,'normal modes', l_nmode, k_nmode))
     &   call errquit('task_vscf: could not allocate l_nmode',nc*nc)
        offset = 0
        if (eaf_read(fd,offset,dbl_mb(k_freq),nc).ne.0)
     &    call errquit('task_vscf:error reading frequencies',555)
        offset = offset + nc
        if (eaf_read(fd,offset,dbl_mb(k_nmode),nc*nc).ne.0)
     &    call errquit('task_vscf:error reading modes',555)
        if (eaf_close(fd).ne.0)
     &    call errquit('task_vscf:file close problem',555)
      endif                                                                        
c
c need to allocate freq and nmode for other nodes and broadcast from
c node 0
c

c
c     Actually do the deed
c
      call vscfm(ngrid,readv,nresmode,dmdr,ncoup,nread,vcfct,iexc,
     &     rtdb,nat,dbl_mb(k_nmode),dbl_mb(k_freq))
c
c     get rid of extra memory
c
      if (.not. ma_pop_stack(l_nmode))
     $  call errquit('task_vscf:ma_pop of l_nmode failed',555)
      if (.not. ma_pop_stack(l_freq))
     $  call errquit('task_vscf:ma_pop of l_freq failed',555)
c     ignore = rtdb_parallel(.true.)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - cpu
c
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_vscf: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_vscf: failed storing walltime',0)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit('task_vscf: failed to set status',0)
c
c
c      call ecce_print1('cpu time', mt_dbl, cpu, 1)
c      call ecce_print1('wall time', mt_dbl, wall, 1)
      task_vscf = status
C      if (task_vscf) then
C        call ecce_print_module_exit('task frequencies', 'ok')
C      else
C        call ecce_print_module_exit('task frequencies', 'failed')
C      endif
c
c
      end
