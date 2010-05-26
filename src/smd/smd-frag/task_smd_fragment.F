      function task_smd_fragment(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical task_smd_fragment
      logical result
c
      character*30 pname
      character*30 operation
      double precision energy
c
      logical task_smd_fragment_optimize
      external task_smd_fragment_optimize
      logical task_smd_fragment_espcharge
      external task_smd_fragment_espcharge
      pname = "task_smd_fragment"
c
      if (.not. rtdb_cget(rtdb,'task:operation', 1,operation))
     $      call errquit(pname,0, RTDB_ERR)
c
      call smd_init("fragment",rtdb)
      call ga_sync()
c
c     close main rtdb as we may be going into groups
c     ----------------------------------------------
      call smd_rtdb_close("keep")
c
      if(inp_compare(.false., 'espcharge', operation)) then
         result=task_smd_fragment_espcharge()
         goto 10
      else if(inp_compare(.false., 'optimize', operation)) then
         result=task_smd_fragment_optimize()
         goto 10
      end if
      energy = 0.0d0
      if(ga_nodeid().eq.0) 
     +    call util_io_enable()
      result = .true.
      call smd_task_iterate()
10    continue
      call smd_rtdb_open(rtdb)
      energy =0.0d0
      if (.not. rtdb_put(rtdb,'fragment:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
      task_smd_fragment = .true.

      call util_message("done with task_smd_fragment_energy")
      if(ga_nodeid().eq.0) then
       call smd_charge_write_default()
      end if
      call smd_global_end()
      call smd_end()
      return
      END

      function task_smd_fragment_espcharge()
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer rtdb
      logical task_smd_fragment_espcharge
      logical result
c
      character*30 pname
      character*30 operation,tag1
c
      double precision energy
      double precision eps,eps1
      integer nc
      integer i
c
      pname = "task_smd_fragment_esp_charge"
c      
      tag1="smd:fragment:chargetol"
      call smd_rtdb_get_dbl(tag1,1,eps,result)
      if(.not.result) eps = 1.0D-4
c      
      tag1="smd:fragment:ncycles"
      call smd_rtdb_get_int(tag1,1,nc,result)
      if(.not.result) nc=10
c
      
      call task_fragment_espcharge_populate("espcharge","t_espcharge")
      call smd_system_set_component("task","t_espcharge")
      do i=1,nc
      call smd_task_iterate()
      call smd_fragment_charge_diff(eps1)
      if(ga_nodeid().eq.0) then
       write(luout,*) "current charge difference is",eps1
       call smd_charge_write_default()
      end if
      if(eps1.le.eps) then
        if(ga_nodeid().eq.0) 
     +        write(luout,*) "achieved charge convergence"
         goto 10
      end if
      end do
10    continue
      return
      END

      function task_smd_fragment_optimize()
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer rtdb
      logical task_smd_fragment_optimize
      logical result
c
      character*30 pname
      character*30 operation,tag1
c
      double precision energy
      double precision eps,eps1
      integer nc
      integer i
      integer gh
c
      pname = "task_smd_fragment_esp_charge"
c      
      tag1="smd:fragment:opttol"
      call smd_rtdb_get_dbl(tag1,1,eps,result)
      if(.not.result) eps = 1.0D-4
c      
      tag1="smd:fragment:ncycles"
      call smd_rtdb_get_int(tag1,1,nc,result)
      if(.not.result) nc=10
c
      call smd_global_create_copy("coords",gh)      
      call task_fragment_optim_populate("optimize","t_opt")
      call smd_system_set_component("task","t_opt")

      do i=1,nc
      call smd_task_iterate()
      call smd_global_diff("coords",gh,eps1)      
      if(ga_nodeid().eq.0) then
       write(luout,*) "current coord difference is",eps1
       call smd_coord_write_default()
      end if
      call smd_coords_save_numbered(i,.false.)
      if(eps1.le.eps) then
        if(ga_nodeid().eq.0) 
     +        write(luout,*) "achieved coord convergence"
         goto 10
      end if
      end do
10    continue
      if(ga_nodeid().eq.0) then
       call smd_coord_write_default()
      end if
       return
      END

      subroutine task_fragment_populate()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "smd_fids.fh"
c     
      character*32 pname
      character*80 tag
      character*80 namespace
      character*80 tag1
      character*30 operation
      integer nt,i1,i2,i,j
      integer i_ifr,i_sop,i_ifu
      logical result
c
      pname = "task_fragment_populate"
c
      call util_message("in "//pname)
c
c     find which monomers to do as given in the input file
c     ----------------------------------------------------
      tag1="smd:fragment:istart"
      call smd_rtdb_get_int(tag1,1,i1,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
      tag1="smd:fragment:iend"
      call smd_rtdb_get_int(tag1,1,i2,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
      nt = i2-i1+1
c
      tag1="task:operation"
      call smd_rtdb_get_string(tag1,1,operation,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
c
c     create task namespace
c     ---------------------
      namespace = "fragment task"
      call smd_namespace_create(namespace)
c
c     create task data structures
c     ---------------------------
      tag = "task:ntasks"
      call smd_data_put_int(namespace,tag,1,nt)
      tag = "task:func_id"
      call smd_data_create_get(namespace,tag,nt,MT_INT,i_ifu)
      tag = "task:frag_id"
      call smd_data_create_get(namespace,tag,nt,MT_INT,i_ifr)
      tag = "task:operation"
      call smd_data_create_get(namespace,tag,nt,MT_STR,i_sop)
c
c    fill up fragment id
c    -------------------
      j=0
      do i=i1,i2
        int_mb(i_ifr+j) = i
        j=j+1
      enddo
c
      do i=1,nt
        int_mb(i_ifu+i-1) = i_smd_monomer
      enddo
c
c     fill in operation
c     -----------------
      call smd_set_char_to_ma2(operation,nt,byte_mb(i_sop))
      call smd_system_set_component("task",namespace)

      call util_message("out "//pname)

      return
      end

      subroutine task_fragment_espcharge_populate(operation,sname)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "smd_fids.fh"
      character*(*) operation
      character*(*) sname
c     
      character*32 pname
      character*80 tag
      character*80 namespace
      character*80 tag1
      integer nt,i1,i2,i,j
      integer i_ifr,i_sop,i_ifu
      logical result
c
      pname = "task_fragment_espcharge_populate"
c
      call util_message("in "//pname)
c
c     find which monomers to do as given in the input file
c     ----------------------------------------------------
      tag1="smd:fragment:istart"
      call smd_rtdb_get_int(tag1,1,i1,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
      tag1="smd:fragment:iend"
      call smd_rtdb_get_int(tag1,1,i2,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
      nt = i2-i1+1
c
c     create task namespace
c     ---------------------
      namespace = "fragment task"
      call smd_namespace_create(sname)
c
c     create task data structures
c     ---------------------------
      tag = "task:ntasks"
      call smd_data_put_int(sname,tag,1,nt)
      tag = "task:func_id"
      call smd_data_create_get(sname,tag,nt,MT_INT,i_ifu)
      tag = "task:frag_id"
      call smd_data_create_get(sname,tag,nt,MT_INT,i_ifr)
      tag = "task:operation"
      call smd_data_create_get(sname,tag,nt,MT_STR,i_sop)
c
c    fill up fragment id
c    -------------------
      j=0
      do i=i1,i2
        int_mb(i_ifr+j) = i
        j=j+1
      enddo
c
      do i=1,nt
        int_mb(i_ifu+i-1) = i_smd_monomer
      enddo
c
c     fill in operation
c     -----------------
      call smd_set_char_to_ma2(operation,nt,byte_mb(i_sop))

      call util_message("out "//pname)

      return
      end

      subroutine task_fragment_optim_populate(operation,sname)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "smd_fids.fh"
      character*(*) operation
      character*(*) sname
c     
      character*32 pname
      character*80 tag
      character*80 namespace
      character*80 tag1
      integer nt,i1,i2,i,j
      integer i_ifr,i_sop,i_ifu
      logical result
c
      pname = "task_fragment_espcharge_populate"
c
      call util_message("in "//pname)
c
c     find which monomers to do as given in the input file
c     ----------------------------------------------------
      tag1="smd:fragment:istart"
      call smd_rtdb_get_int(tag1,1,i1,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
      tag1="smd:fragment:iend"
      call smd_rtdb_get_int(tag1,1,i2,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
      nt = i2-i1+1
c
c     create task namespace
c     ---------------------
      namespace = "fragment task"
      call smd_namespace_create(sname)
c
c     create task data structures
c     ---------------------------
      tag = "task:ntasks"
      call smd_data_put_int(sname,tag,1,nt)
      tag = "task:func_id"
      call smd_data_create_get(sname,tag,nt,MT_INT,i_ifu)
      tag = "task:frag_id"
      call smd_data_create_get(sname,tag,nt,MT_INT,i_ifr)
      tag = "task:operation"
      call smd_data_create_get(sname,tag,nt,MT_STR,i_sop)
c
c    fill up fragment id
c    -------------------
      j=0
      do i=i1,i2
        int_mb(i_ifr+j) = i
        j=j+1
      enddo
c
      do i=1,nt
        int_mb(i_ifu+i-1) = i_smd_monomer
      enddo
c
c     fill in operation
c     -----------------
      call smd_set_char_to_ma2(operation,nt,byte_mb(i_sop))

      call util_message("out "//pname)

      return
      end
