c
c $Id: qmmm.F,v 1.62 2005-09-12 01:28:27 marat Exp $
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
      
      if(operation.eq."energy") then
        ostatus = task_qmmm_energy(rtdb) 
      else if (operation.eq."gradient") then
        ostatus = task_qmmm_gradient(rtdb)
      else if (operation.eq."dynamics") then
        ostatus = task_qmmm_dynamics(rtdb)
      end if 
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
      logical ignore
      logical checkgrad
      character*32 pname
      character*32 theory
      character*30 operation
      character*255 geom
      logical aux_geom
      logical loadgeom
      logical ouselbfgs
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

      ignore = MA_set_hard_fail(.true.)
    
      qmmm_init = .true.

      call qmmm_data_rdinp(irtdb)
c     initialize mm part
c     (note it will also want operation)
c     --------------------------------------------
      call mm_init(irtdb)

      call qmmm_cons_init(irtdb)
      call qmmm_links_init(irtdb)

      call qmmm_data_init(irtdb)

      call qmmm_bq_data_init(irtdb)

      if (.not. rtdb_cget(irtdb, 'task:operation', 1, operation))
     $     operation = 'energy'

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
      else
        aux_geom = .false.
      end if

      if (inp_compare(.false.,'esp',theory)) aux_geom=.false.

c      if (rtdb_cget(irtdb,"qmmm:geom",1,geom)) then
c        loadgeom = .true.
c        aux_geom = .true.
c      else
c        loadgeom = .false.
c      end if

      if(oprint) 
     >  write(*,*) "aux_geom",aux_geom

      call qmmm_geom_init(irtdb,aux_geom)

      if(loadgeom) call md_sp()

      if (.not. rtdb_put(irtdb,'qmmm:init',mt_log,1,.true.))
     $     call errquit('qmmm: failed qmmm:init', 0, RTDB_ERR)


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



