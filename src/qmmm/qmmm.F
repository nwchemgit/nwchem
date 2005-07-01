c
c $Id: qmmm.F,v 1.51 2005-07-01 01:20:55 marat Exp $
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
      logical ignore
      logical checkgrad
      character*32 pname
      character*30 operation
      character*30 optimization
      character*255 geom
      logical aux_geom
      logical loadgeom
c
      call qmmm_print_header()
c
      pname = "qmmm_init"

      call util_print_push() 
      call util_print_rtdb_load(irtdb, 'qmmm')
c
      if(qmmm_print_debug()) 
     >  write(*,*) ga_nodeid(),"in "//pname

      ignore = MA_set_hard_fail(.true.)
    
      qmmm_init = .true.

      call qmmm_data_rdinp(irtdb)
c     initialize mm part
c     (note it will also want operation)
c     --------------------------------------------
      call mm_init(irtdb)

      call qmmm_cons_init()
      call qmmm_links_init(irtdb)

      call qmmm_data_init(irtdb)

      call qmmm_bq_data_init(irtdb)

      if (.not. rtdb_cget(irtdb, 'task:operation', 1, operation))
     $     operation = 'energy'

      
      if(operation.eq."optimize") then
        aux_geom = .true.
      else if(operation.eq."freq") then
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
        aux_geom = .true.
      end if

      if (.not. rtdb_cget(irtdb, 'qmmm:optimization', 1, optimization))
     $     call errquit('qmmm_optimize: failed ',0,0)

      if (optimization.ne.'qm') aux_geom = .false.
      

      if (rtdb_cget(irtdb,"qmmm:geom",1,geom)) then
        loadgeom = .true.
        aux_geom = .true.
      else
        loadgeom = .false.
      end if

      call qmmm_geom_init(irtdb,aux_geom)

      if(loadgeom) call md_sp()

      if (.not. rtdb_put(irtdb,'qmmm:init',mt_log,1,.true.))
     $     call errquit('qmmm: failed qmmm:init', 0, RTDB_ERR)


      if(qmmm_print_debug()) 
     >  write(*,*) ga_nodeid(),"out "//pname


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
      call qmmm_geom_data_release()
      call qmmm_bq_data_dealloc()
      call qmmm_data_release()
      call qmmm_link_data_release()

      if (.not. rtdb_put(rtdb,'qmmm:init',mt_log,1,.false.))
     $     call errquit('qmmm: failed qmmm:init', 0, RTDB_ERR)

      qmmm_end = .true.
      return
      end
c



