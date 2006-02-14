c
c $Id: qmmm.F,v 1.67 2006-02-14 01:08:56 marat Exp $
c
c
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

    
      qmmm_init = .true.

      if (.not. rtdb_cget(irtdb, 'task:operation', 1, operation))
     $     operation = 'energy'

      if(operation.eq."thermodynamics") then
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

      if(operation.eq."freq") then
        aux_geom = .true.
      else if(operation.eq."hessian") then
        aux_geom = .true.
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

      integer rtdb
      logical qmmm_end
      qmmm_end = .true.
c
      call mm_end()
      call qmmm_geom_end()
      call qmmm_bq_data_dealloc()
      call qmmm_data_release()
      call qmmm_link_data_release()
      call qmmm_cons_end()

      if (.not. rtdb_put(rtdb,'qmmm:init',mt_log,1,.false.))
     $     call errquit('qmmm: failed qmmm:init', 0, RTDB_ERR)

      qmmm_end = .true.
      return
      end
c



