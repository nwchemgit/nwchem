c
c $Id$
c
C> \ingroup task
C> @{
      function task_qmmm_property(rtdb)
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
      logical task_qmmm_property
c
      logical  task_property
      external task_property
c
      logical status
      character*255 vectors_in
      character*30 pname 
c
      pname = "task_qmmm_property"
c
c     create full geometry
c     -------------------
      call qmmm_geom_create_full(rtdb)
c
c     update coordinates of classical charges
c     ---------------------------------------
      call qmmm_bq_coord_update()
c
c     no checking if external vectors
c     -------------------------------
      if (.not. rtdb_cget(rtdb,'prop:vectors',1,vectors_in))
     &    vectors_in = " "

      if(vectors_in.ne." ") then
        if (.not. rtdb_put(rtdb,'prop:nocheck',
     >                   mt_log,1,.true.))
     &    call errquit(pname//'prop:nocheck failed',911,RTDB_ERR)
      end if

      status = task_property(rtdb)
      task_qmmm_property = status
      return
      end
C> @}



