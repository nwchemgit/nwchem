c
c $Id$
c
C> \ingroup task
C> @{
      function task_qmmm_generic(rtdb,operation)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
c
      integer rtdb
      character*(*) operation
      logical  task_qmmm_generic
c
      logical  task_qmmm_gradient
      external task_qmmm_gradient
c
      logical  task_qmmm_energy
      external task_qmmm_energy
c
      logical  task_qmmm_optimize
      external task_qmmm_optimize
c
      logical  task_qmmm_dynamics
      external task_qmmm_dynamics
c     local variables
      logical ostatus
      
      ostatus = .false.
      if(operation.eq."energy") then
        ostatus = task_qmmm_energy(rtdb) 
      else if (operation.eq."gradient") then
        ostatus = task_qmmm_gradient(rtdb)
      else if (operation.eq."dynamics") then
        ostatus = task_qmmm_dynamics(rtdb)
      end if 
      task_qmmm_generic = ostatus
      end
C> @}
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
      logical checkgrad
      character*32 pname
      character*32 theory
      character*30 operation
      character*30 operation0
      logical aux_geom
      logical oprint
      logical oinit
      logical qmmm_active_region
      logical numerical
c
c     do nothing if already initialized
c     ---------------------------------
      if (.not. rtdb_get(irtdb,'qmmm:init',mt_log,1,qmmm_init))
     $     qmmm_init = .false.
      if(qmmm_init) return
c
      call qmmm_print_header()
c
      pname = "qmmm_init"

      call util_print_push() 
      call util_print_rtdb_load(irtdb, 'qmmm')
c
      oprint = util_print("debug_info", print_debug)
      oprint = oprint .and. (ga_nodeid().eq.0)
c
      if(oprint) 
     >  write(*,*)"in "//pname
c
   
      qmmm_init = .true.

      if (.not. rtdb_cget(irtdb, 'task:operation', 1, operation))
     $     operation = 'energy'

      if(operation.eq."fep") then
        operation0 = "dynamics" 
      else
        operation0=operation
      end if

      if (.not. rtdb_cput(irtdb, 'task:operation', 1, operation0))
     $    call errquit(pname,0, RTDB_ERR)

      call qmmm_data_rdinp(irtdb)
c     initialize mm part
c     (note it will also want operation)
c     --------------------------------------------
      call mm_init(irtdb)

      call qmmm_cons_init(irtdb)
      call qmmm_links_init(irtdb)

      call qmmm_data_init(irtdb)

      call qmmm_bq_data_init(irtdb)


      if(.not.rtdb_cget(irtdb,'task:theory',1,theory))
     + call errquit('task: no task input for theory?',0, RTDB_ERR)

      qmmm_active_region=.false.
      numerical=.false.
      if(operation.eq."freq") then
        aux_geom = .true.
        numerical = .true.
        qmmm_active_region=.true.
      else if(operation.eq."hessian") then
        aux_geom = .true.
        numerical = .true.
        qmmm_active_region=.true.
      else if(operation.eq."gradient") then
        aux_geom = .false.
      else if(operation.eq."energy") then
        aux_geom = .false.
      else if(operation.eq."property") then
        aux_geom = .false.
      else if(operation.eq."neb") then
        aux_geom = .false.
      else
        aux_geom = .false.
      end if

      if (inp_compare(.false.,'esp',theory)) aux_geom=.false.

      if(oprint) 
     >  write(*,*) "aux_geom",aux_geom
c
c     this will allow frozen/esp density treatment of fixed QM region
c     ---------------------------------------------------------------
      if(qmmm_active_region) then
      if (.not. rtdb_put(irtdb,'qmmm:active_region',mt_log,1,.true.))
     $     call errquit('qmmm: failed qmmm:init', 0, RTDB_ERR)
      end if
      if (numerical) then
       if (.not. rtdb_put(irtdb, 'task:numerical', mt_log, 1, 
     $        numerical)) call errquit('qmmm: rtdb
     $             is corrupt',0, RTDB_ERR)
      endif


      call qmmm_geom_init(irtdb,aux_geom)
      if(operation.eq."neb") then
         call qmmm_geom_init_neb(irtdb)
      end if

      call md_sp()

      if (.not. rtdb_put(irtdb,'qmmm:init',mt_log,1,.true.))
     $     call errquit('qmmm: failed qmmm:init', 0, RTDB_ERR)

      if (.not. rtdb_cput(irtdb, 'task:operation', 1, operation))
     $    call errquit(pname,0, RTDB_ERR)

      if(oprint) 
     >  write(*,*) "out "//pname


      if (.not. rtdb_get(irtdb,'qmmm:checkgrad',mt_log,1,checkgrad))
     $     checkgrad = .false.

      if(checkgrad) then
       call qmmm_check_forces(irtdb)
      end if

      call util_print_pop() 

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
#include "global.fh"

      integer rtdb
      logical qmmm_end
c
      character*32 pname
      logical oinit
      logical oprint
c
      logical qmmm_print_debug
      external qmmm_print_debug
c
      pname = "qmmm_end"
      if(qmmm_print_debug())
     >  write(*,*)"in "//pname
      qmmm_end = .true.
c
      if (.not. rtdb_get(rtdb,'qmmm:init',mt_log,1,oinit))
     $     oinit = .false.
c
      if(oinit) then
        call mm_end()
        call qmmm_geom_end()
        call qmmm_bq_data_dealloc()
        call qmmm_data_release()
        call qmmm_link_data_release()
        call qmmm_cons_end()
      end if

      
      if (.not. rtdb_put(rtdb,'qmmm:init',mt_log,1,.false.))
     $     call errquit('qmmm: failed qmmm:init', 0, RTDB_ERR)

      qmmm_end = .true.
      if(qmmm_print_debug())
     >  write(*,*)"out "//pname
      call ga_sync()
      return
      end
c



