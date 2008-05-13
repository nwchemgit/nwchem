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
      integer t_params(100)
c
      pname = "smd_tasklist_init"      
c
      smd_ntask = 11
      if(.not.nga_create(MT_INT, 1, smd_ntask, 
     >    'smd_ga_task_params',0, smd_ga_task_params)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)

      t_params(1) = 1
      do i=2,11
       t_params(i) = i-1
      end do

      call nga_put(smd_ga_task_params,1, 11,t_params,1)
      call ga_print(smd_ga_task_params)


      return
      end



