c
c $Id: qmmm.F,v 1.38 2004-12-04 02:59:27 marat Exp $
c
c
      function qmmm_init(irtdb)
      implicit none
c
#include "inp.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "global.fh"
      integer irtdb
      logical qmmm_init
c
      character*32 pname
c
      call qmmm_print_header(irtdb)
c
      pname = "qmmm_init"

      call util_print_push() 
      call util_print_rtdb_load(irtdb, 'qmmm')
c
      if(qmmm_print_debug()) 
     >  write(*,*) ga_nodeid(),"in "//pname
    
      qmmm_init = .true.

      call qmmm_data_rdinp(irtdb)
c     initialize mm part
c     (note it will also want operation)
c     --------------------------------------------
      call mm_init(irtdb)
      call qmmm_init_qm(irtdb)

      if (.not. rtdb_put(irtdb,'qmmm:init',mt_log,1,.true.))
     $     call errquit('qmmm: failed qmmm:init', 0, RTDB_ERR)


      if(qmmm_print_debug()) 
     >  write(*,*) ga_nodeid(),"out "//pname

      call util_print_pop() 

      return
      end
c
      subroutine qmmm_init_qm(irtdb)
      implicit none
c
#include "inp.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
      integer irtdb
c
c     allocate and initialize qmmm data structure
c     -------------------------------------------
      call qmmm_data_init(irtdb)

c     initialize bq settings 
      call qmmm_bq_data_init(irtdb)

c     create initial qm geometry
c     --------------------------
      call qmmm_geom_init(irtdb)

      return
      end
c
      function qmmm_end(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      integer rtdb
      logical qmmm_end
      qmmm_end = .true.
c
      call mm_end()
      call qmmm_data_release()
      call qmmm_bq_data_dealloc()

      qmmm_end = .true.
      return
      end
c
      function qmmm_energy(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
#include "global.fh"
c
      integer rtdb
      logical qmmm_energy
c     local variables
      logical status
      character*255 theory
      double precision energy
      double precision qm_energy
      double precision mm_energy
      double precision eatoms
c     external functions
      logical ignore_mm
      logical ignore_qm
      logical  task_energy_doit
      external task_energy_doit

      call util_print_push() 
      call util_print_rtdb_load(rtdb, 'qmmm')
      
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

       if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)
      call qmmm_print_current_geom(rtdb)

c
c     get qm coords out of geometry object into qmmm module
c     ----------------------------------------------------
      call qmmm_get_qm_coord_qm(rtdb)
c
c     push coords into mm module
c     --------------------------
      call qmmm_push_qm_coord_mm(rtdb)
      call mm_task_sp()
c
c     create full geometry
c     -------------------
      call qmmm_create_full_geom(rtdb)

c     QM energy calculation
c     --------------------
      status = task_energy_doit(rtdb,theory,qm_energy)
c
      call qmmm_restore_qm_geom(rtdb)

c     combine QM and MM energies
c     -------------------------
      mm_energy=0.0d0
      if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)
      mm_energy = mm_energy/cau2kj
      
      if (.not. rtdb_get(rtdb,'qmmm:ignore_mm',mt_log,1,ignore_mm)) then
         ignore_mm = .false.
      end if

      if (.not. rtdb_get(rtdb,'qmmm:ignore_qm',mt_log,1,ignore_qm)) then
         ignore_qm = .false.
      end if

      if(ignore_mm) mm_energy = 0.0d0
      if(ignore_qm) qm_energy = 0.0d0

      energy = qm_energy + mm_energy-eatoms

      if (.not. rtdb_put(rtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'qmmm:mm_energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)


      call qmmm_print_energy(rtdb)

      qmmm_energy=status

      call util_print_pop()
      call mm_write_restart()

      end



