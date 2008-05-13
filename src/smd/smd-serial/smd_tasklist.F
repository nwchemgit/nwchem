      subroutine smd_tasklist_init()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_tasklist_data.fh"
      character*30 pname 
      integer i
      integer t(100)
c
      pname = "smd_tasklist_init"      
c
      smd_ntask = 10
      if(.not.nga_create(MT_INT, 1, smd_ntask, 
     >    'smd_ga_task_params',0, smd_ga_task_params)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)
c
      if(.not.nga_create(MT_INT, 1, smd_ntask, 
     >    'smd_ga_task_status',0, smd_ga_task_status)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)

      do i=1,10
       t(i) = i
      end do

      call nga_put(smd_ga_task_params,1, 10,t,1)
      call ga_print(smd_ga_task_params)

      do i=1,10
       t(i) = t_avail
      end do

      call nga_put(smd_ga_task_status,1, 10,t,1)
      call ga_print(smd_ga_task_status)


      return
      end

      function smd_tasklist_get(i)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_tasklist_data.fh"
      character*30 pname 
      integer i
      integer t_p(100)
      integer t_s(100)
      logical smd_tasklist_get
c
      pname = "smd_tasklist_init"      
c
      smd_ntask = 10

      call nga_get(smd_ga_task_status,1, 10,t_s,1)
      call ga_sync()

      smd_tasklist_get =.false.
      do i=1,smd_ntask
        if(t_s(i).eq.t_avail) goto 100
      end do
      return
100   continue
      t_s(i) = t_busy
      call ga_sync()
      call nga_put(smd_ga_task_status,i,i,t_s,1)
      smd_tasklist_get =.true.
      return
      end

