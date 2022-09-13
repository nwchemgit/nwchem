c
c     Task handler for FCIDUMP
c
      logical function task_fcidump(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
#include "geom.fh"
#include "global.fh"
c
      integer rtdb
c
      logical status
      double precision cputime, walltime
      character*32 theory
      logical molecorb
      integer occ_frozen,vir_frozen,geom
      character*255 filename
      integer iunit,ios
c
c     ! Time at the start
      cputime  = util_cpusec()
      walltime = util_wallsec()
c     hack forcing dft
c      if (.not. rtdb_cput(rtdb, 'task:theory', 1, 'dft'))
c     &  call errquit('task_fcidump: theory not specified',0, RTDB_ERR)
      
c
c     only dft works for now
c      
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     &  call errquit('task_fcidump: theory not specified',0, RTDB_ERR)
      if (theory.eq.'dft') then
         status = .true.
      else
         status = .false.
         write(luout,*) ' fcidump for theory ',
     &        theory(1:inp_strlen(theory)), ' not implemented yet'
         call errquit(' fcidump stop',0,0)
      end if

c     ! Perform the task
c      call fcidump_input(rtdb)
      if (.not. rtdb_get(rtdb, 'fcidump:mos', mt_log, 1, molecorb))
     &        molecorb=.true.
      occ_frozen=0
      vir_frozen=0
      if(molecorb) then
c     uses moints_full
c     get frozen definitions
c
c     occupied
c         
      if (.not. geom_create(geom,'geometry'))
     $     call errquit('moints: geom failed',0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     $     call errquit('moints: geometry not found',0,
     &       RTDB_ERR)
      if(.not. geom_freeze(rtdb, geom, 'fcidump', occ_frozen)) then
        call errquit('fcidump: unable to freeze core orbitals', 0,
     &       INPUT_ERR)
      endif
      if (.not. geom_destroy(geom))
     $     call errquit('fcidump: geom corrupt',0, GEOM_ERR)

      if (.not. rtdb_put(rtdb,'fourindex:occ_frozen',MT_INT,1,
     $     occ_frozen))
     $   call errquit('fcidump: rtdb failure for occ_frozen',555,
     &     RTDB_ERR)
c     virtual
      if (.not.rtdb_get(rtdb, 'fcidump:frozen virtual', MT_INT, 1,
     $     vir_frozen)) then
         if (.not.rtdb_get(rtdb,'vir_frozen',MT_INT,1,vir_frozen))
     $        vir_frozen = 0
      endif
      if (.not.rtdb_put(rtdb,'fourindex:vir_frozen',MT_INT,1,
     $     vir_frozen))
     $   call errquit('fcidump: rtdb failure for occ_virtual',555,
     &     RTDB_ERR)
      endif
      if (.not. rtdb_put(rtdb,'fourindex:osym',MT_LOG,1,
     $     .false.))
     $   call errquit('fcidump: rtdb failure for :osym',557,
     &     RTDB_ERR)
      call util_file_name('fcidump', .false., .false.,filename)
      iunit=1948
      if(ga_nodeid().eq.0) then
         open(unit=iunit,
     1        file=filename(1:inp_strlen(filename)),
     2        status='unknown',
     3        form='formatted',
     4        iostat=ios)
         if(ios.ne.0) then
            write(luout,*) ' opening ',filename(1:inp_strlen(filename)),
     c           ' failed '
            task_fcidump=.false.
            return
         else
            write(luout,*) ' writing FCIDUMP file ',
     c           filename(1:inp_strlen(filename))
            write(luout,*) '-------------------------------------------'
         endif
      endif
cTODO need to broadcast io failure      

         call fcidump_header(rtdb,occ_frozen,vir_frozen,iunit)
      if(ga_nodeid().eq.0) then
         close(unit=iunit,
     1        status='keep',
     2        iostat=ios)
      endif
      
      if(molecorb) then
         if (.not. rtdb_cput(rtdb,'fourindex:file format',1, 'fcidump'))
     $   call errquit('fcidump: rtdbcput failure ',5, RTDB_ERR)
         if (.not.rtdb_cput(rtdb,'fourindex:moints file',1,
     F        filename(1:inp_strlen(filename))))
     $   call errquit('fcidump: rtdbcput failure ',2, RTDB_ERR)
         call moints_full(rtdb)
      else
         call errquit(' atomic not implemented yet ',0,CAPMIS_ERR)
      endif
      
      task_fcidump=status
      
      return
 911  continue
      end
