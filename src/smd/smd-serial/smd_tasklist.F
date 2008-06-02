      subroutine smd_tasklist_init()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_tasklist_data.fh"
      character*30 pname, theory
      integer i
      integer istart,iend
      integer t(100)
      logical result
      character*72 tag
c
      pname = "smd_tasklist_init"      
c
      tag = "smd:theory"
      call smd_rtdb_get_string(tag,1,theory,result)
      if(.not.result)
     > call errquit(
     >       pname//'no '//tag,0, RTDB_ERR)
c
      tag = "smd:fragment_istart"
      call smd_rtdb_get_int(tag,1,istart,result)
      if(.not.result)
     > call errquit(
     >       pname//'no '//tag,0, RTDB_ERR)
c
      tag = "smd:fragment_iend"
      call smd_rtdb_get_int(tag,1,iend,result)
      if(.not.result)
     > call errquit(
     >       pname//'no '//tag,0, RTDB_ERR)
c
c
      smd_ntask = iend-istart+1
c
      call smd_rtdb_get_log("smd:steal",1,osteal,result)
      if(.not.result) osteal = .false.

      if(.not.nga_create(MT_INT, 1, smd_ntask, 
     >    'smd_ga_task_params',0, smd_ga_task_params)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)
c
      if(.not.nga_create(MT_INT, 1, smd_ntask, 
     >    'smd_ga_task_status',0, smd_ga_task_status)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)
c
c      if(.not.nga_create(MT_INT, 1, smd_ntask, 
c     >    'smd_ga_group_pool',0, smd_ga_group_pool)) 
c     >  call errquit(
c     >       pname//'nga_create error ',0, RTDB_ERR)


      do i=1,smd_ntask 
        t(i) = istart+i-1
      end do

      call nga_put(smd_ga_task_params,1,smd_ntask,t,1)
      call ga_print(smd_ga_task_params)

      do i=1,smd_ntask
       t(i) = t_avail
      end do

      call nga_put(smd_ga_task_status,1,smd_ntask,t,1)
      do i=1,smd_ntask
       t(i) = t_busy
      end do

c      call nga_put(smd_ga_group_pool,1,smd_ntask,t,1)
c      call ga_print(smd_ga_task_status)
c      call ga_print(smd_ga_task_status)
 
      return
      end

      function smd_tasklist_get1(i)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_tasklist_data.fh"
#include "smd_group_data.fh"
      character*30 pname 
      integer i
      integer t_p(100)
      integer t_s(100)
      logical smd_tasklist_get1
c
      pname = "smd_tasklist_get"      
c
c      do i=1,smd_group_id
      call ga_lock(smd_group_id)
      call nga_get(smd_ga_task_status,1, smd_ntask,t_s,1)
      call ga_sync()
c      end do

      smd_tasklist_get1 =.false.
      do i=1,smd_ntask
        if(t_s(i).eq.t_avail) goto 100
      end do
      return
100   continue
      t_s(i) = t_busy
      call ga_sync()
      call nga_put(smd_ga_task_status,i,i,t_s,1)
      call ga_unlock(smd_group_id)
      smd_tasklist_get1 =.true.
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
#include "smd_group_data.fh"
      character*30 pname 
      integer t
      integer i,ig,j,ng,it
      integer t_p(100)
      integer t_s(100)
      logical smd_tasklist_get
c
      pname = "smd_tasklist_get"      
c
      ig = smd_group_id
      ng = smd_ngroup
c      write(90+ig,*) "getting task"
c      write(*,*) "ig,ng,ntask",ig,ng,smd_ntask
c      call nga_put(smd_ga_group_pool,ig,ig,t_avail,1)
      call nga_get(smd_ga_task_status,1, smd_ntask,t_s,1)
      call ga_sync()

c      write(90+ig,*) "task array"
c      do i=1,smd_ntask 
c      write(90+ig,*) i,t_s(i).eq.t_avail
c      end do
 
      smd_tasklist_get =.false.
c      write(90+ig,*) "looking for unfinished task"
c     honestly look for unfinished task
      do i=ig,smd_ntask,ng
        if(t_s(i).eq.t_avail) goto 100
      end do
c
      if(.not.osteal) return
c     start stealing work from left cpus
c      call nga_get(smd_ga_task_status,1, smd_ntask,t_s,1)
      call ga_sync()
      do j=1,ig-1
       it = 0
       do i=j,smd_ntask,ng
        if(t_s(i).eq.t_avail) then
         it = it + 1
         if(it.eq.2) goto 100
        end if
       end do
      end do
c     start stealing work from right cpus
c      call nga_get(smd_ga_task_status,1, smd_ntask,t_s,1)
      call ga_sync()
      do j=ig+1,ng
       it = 0
       do i=j,smd_ntask,ng
        if(t_s(i).eq.t_avail) then
         it = it + 1
         if(it.eq.2) goto 100
        end if
       end do
      end do
      return
100   continue
c      call nga_put(smd_ga_group_pool,ig,ig,t_busy,1)
      call nga_put(smd_ga_task_status,i,i,t_busy,1)
      call ga_sync()
      smd_tasklist_get =.true.
      return
      end

      subroutine smd_tasklist_set()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_tasklist_data.fh"
#include "smd_group_data.fh"
      character*30 pname 
      integer t
      integer i,ig,j,ng,it
      integer t_p(100)
      integer t_s(100)
      integer p_world,p_default
c
      pname = "smd_tasklist_set"      
c
      p_world   = ga_pgroup_get_world()
      p_default = ga_pgroup_get_default()
      call ga_pgroup_sync(p_world)
c
      call ga_pgroup_set_default(p_world)
c
      if(ga_nodeid().eq.0) then
      do i=1,smd_ntask
        t_s(i)=t_avail
      end do
      call nga_put(smd_ga_task_status,1, smd_ntask,t_s,1)
      end if

      call ga_pgroup_set_default(p_default)

      return
      end
