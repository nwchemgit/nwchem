*
* $Id$
*
C> \ingroup task
C> @{
C>
C> \brief Generic interface to perform Vibrational SCF (VSCF) calculations
C>
C> VSCF calculations explore the anharmonicity of the potential energy
C> surface using finite displacements along two normal modes simultanously.
C>
C> \return Return .true. if successfull, and .false. otherwise.
C>
      logical function task_vscf(rtdb)
c     
c     Also will be hooked up to analytic methods as they are available.
c
c     Since this routine is directly invoked by application modules
c     no input is processed in this routine.
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
#include "global.fh"
#include "eaf.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
c
      logical task_freq
      external task_freq
c     
      integer rtdb !< [Input] The RTDB handle
c     
      integer nat, geom, i
      integer ngrid, ncoup, iexc
      integer l_freq, k_freq, l_nmode, k_nmode, nc
      integer fd
      double precision offset
      character*255 filename
      character*32 theory
      double precision vcfct
      logical status, dmdr, freq_done
c     
      double precision cpu, wall
c 
c  set up a couple of variables that I will eventually let be user controlled
c
      dmdr = .false.
c    
      call ecce_print_module_entry('task vscf')
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_vscf: failed to invalidate status',0,
     &       RTDB_ERR)
      if (ga_nodeid().eq.0 .and.
     $    util_print('task_vscf', print_low)) then
        write(LuOut,*)
        write(LuOut,*)
        call util_print_centered(LuOut,'VSCF Analysis',40,.true.)
        write(LuOut,*)
      endif
c
c     First do a frequency analysis. Check if restart and frequencies
c     were already done
c
      if (.not. rtdb_get(rtdb,'vscf:freq_done',mt_log,1,freq_done))
     $      freq_done=.false.
      status = .true.
      if (.not.freq_done) status = task_freq(rtdb)
      if (.not.status) then
         call errquit('task_vscf: task_freq failed',911, CALC_ERR)
      else
         if (.not. rtdb_put(rtdb,'vscf:freq_done',mt_log,1,freq_done))
     $   call errquit('task_vscf: rtdb_put freq_done',911, RTDB_ERR)
      endif
c
c     Create/load reference geometry to get the number of atoms
c
      if (.not.geom_create(geom,'geometry')) call errquit
     $      ('task_vscf:geom_create failed?',1, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     $      call errquit
     $      ('task_vscf:geom_rtdb_load failed?',2, RTDB_ERR)
      if (.not. geom_ncent(geom,nat)) call errquit
     $      ('task_vscf:geom_ncent failed?',3, GEOM_ERR)
      nc = nat*3
      if (.not. geom_destroy(geom)) call errquit
     $      ('task_vscf:geom_destroy failed?',911, GEOM_ERR)
c
c     Load parameters from rtdb or assign defaults
c
      if (.not. rtdb_get(rtdb,'vscf:ngrid',mt_int,1,ngrid))
     $      ngrid = 16
      if (.not. rtdb_get(rtdb,'vscf:iexc',mt_int,1,iexc))
     $      iexc = 1
      if (.not. rtdb_get(rtdb,'vscf:ncoup',mt_int,1,ncoup))
     $      ncoup = 2
      if (.not. rtdb_get(rtdb,'vscf:vcfct',mt_dbl,1,vcfct))
     $      vcfct = 1.0d+0
c
c  allocate space for freq and normal modes
c
      if (.not. ma_push_get(mt_dbl, nc,'frequencies', l_freq, k_freq))
     &    call errquit('task_vscf: could not allocate l_freq',nc,MA_ERR)
      if (.not. 
     &  ma_push_get(mt_dbl, nc*nc,'normal modes', l_nmode, k_nmode))
     &   call errquit('task_vscf: could not allocate l_nmode',nc*nc,
     &       MA_ERR)
c
c  read in frequencies and normal modes
c
      if (ga_nodeid().eq.0) then
        call util_file_name('nmode',.false.,.false.,filename)
        if (eaf_open(filename,eaf_rw,fd).ne.0)
     &    call errquit('task_vscf:eaf_open failed',555, DISK_ERR)
        offset = 16 
        if (eaf_read(fd,offset,dbl_mb(k_freq),nc*8).ne.0)
     &     call errquit('task_vscf:error reading frequencies',555,
     &       DISK_ERR)
        offset = offset + nc*8
        if (eaf_read(fd,offset,dbl_mb(k_nmode),nc*nc*8).ne.0)
     &    call errquit('task_vscf:error reading modes',555, DISK_ERR)
        if (eaf_close(fd).ne.0)
     &    call errquit('task_vscf:file close problem',555, DISK_ERR)
      endif                                                       
c
c need to allocate freq and nmode for other nodes and broadcast from
c node 0
c
      call ga_brdcst(msg_frq,dbl_mb(k_freq),
     &               nc*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(msg_eigv,dbl_mb(k_nmode),
     &               nc*nc*ma_sizeof(mt_dbl,1,mt_byte),0)
c
c     Actually do the deed
c
      call vscfm(ngrid,dmdr,ncoup,vcfct,iexc,
     &     rtdb,nat,dbl_mb(k_nmode),dbl_mb(k_freq))
c
c     get rid of extra memory
c
      if (.not. ma_pop_stack(l_nmode))
     $  call errquit('task_vscf:ma_pop of l_nmode failed',555, MA_ERR)
      if (.not. ma_pop_stack(l_freq))
     $  call errquit('task_vscf:ma_pop of l_freq failed',555, MA_ERR)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - cpu
c
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_vscf: failed storing cputime',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_vscf: failed storing walltime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit('task_vscf: failed to set status',0, RTDB_ERR)
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      task_vscf = status
      if (task_vscf) then
        call ecce_print_module_exit('task vscf', 'ok')
      else
        call ecce_print_module_exit('task vscf', 'failed')
      endif
c
c
      end
C> @}
