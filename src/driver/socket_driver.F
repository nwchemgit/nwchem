      logical function socket_driver(rtdb)
      implicit none
      integer          rtdb

#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
#include "msgtypesf.h"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,oprint,done
      character*80  stype,socket_ip,socket_movecs
      character*80 statebuffer
      character*80 buffer
      character*20 msg,ip,pp,bead_list,geom_name
      integer ii,it,port,iter,n1,n2,inet
      integer option,tmpint,nbytes,sock
      integer rion(2),fion(2),nion,nion0
      real*8 unita(3,3),invunita(3,3),stress(3,3),energy
      real*8 cpu1,cpu2

*     **** external functions ****
      logical  task_gradient
      external task_gradient
      integer  nion_bead_list,size_bead_list
      external nion_bead_list,size_bead_list
      real*8   energy_bead_list
      external energy_bead_list

      value = .true.

      oprint = (ga_nodeid().eq.0)

      !call util_file_prefix('movecs',socket_movecs)
      socket_movecs = "socket.movecs"
      bead_list = 'asebead'
      geom_name = 'geometry'
      call init_bead_list(rtdb,bead_list,socket_movecs)
      call add_bead_list(bead_list,socket_movecs,geom_name)
      nion = nion_bead_list(bead_list,1)
      n1 = 9*8
      n2 = 3*8*nion

      !**** allocate rion, fion ****
      value = MA_push_get(mt_dbl,6*nion,'rion',rion(2),rion(1))
      value = value.and.
     >        MA_push_get(mt_dbl,6*nion,'fion',fion(2),fion(1))
      if (.not.value)
     >   call errquit('socket_driver:stack failed',1,RTDB_ERR)


      if (oprint) then
         write(luout,*)
         write(luout,*)
         call util_print_centered(luout,
     >        'NWChem Socket Interface',
     >        40,.true.)
         write(luout,*)
         write(luout,*)
      endif

      if (.not.rtdb_cget(rtdb,'driver:socket_type',1,stype))
     >  call errquit('driver_input: rtdb get failed',7,RTDB_ERR)

      if (.not.rtdb_cget(rtdb,'driver:socket_ip',1,socket_ip)) then
          if (inp_compare(.false.,stype,'unix')) then
             socket_ip = "nwchem"
         else
             socket_ip = "127.0.0.1:31415"
         end if
      end if
      
      if (inp_compare(.false.,stype,'unix')) then
          ip = socket_ip
          inet = 0
          ii = index(socket_ip,' ')
          port = 0
      else
          inet = 1
          ii = index(socket_ip,':')
          it = index(socket_ip,' ') - 1
          pp(1:it-ii) = socket_ip(ii+1:it)
          pp(it-ii+1:it-ii+5) = '    '
          read(pp,*,err=300) port
          go to 301
 300      port = 31415
 301      continue
          ip(1:ii-1) = socket_ip(1:ii-1)
          ip(ii:ii+4) = "    "
      end if

      taskid = ga_nodeid()

      statebuffer(1:13) = "READY        "

      if (taskid.eq.MASTER) then
         write(luout,*)
         write(luout,*) "== i-PI Socket Client Driver =="
         write(luout,'(" Connected to    = ",A)') socket_ip
         write(luout,'(" Number of atoms =",I8)') nion
         write(luout,*)
         call util_talker(ip,inet,ii-1,port,sock)
         write(luout,*)
         !call nwpw_talker("127.0.0.1",9,port,sock)
      end if
      done = .false.
      iter = 0
      do while (.not.done)
         msg = "            "

         !**** read option ****
         option = 0
         if (taskid.eq.MASTER) then
             call util_talker_read(sock,msg,12)
             if (inp_compare(.false.,msg(1:4),'EXIT'))  then
                option = 0
             else if (inp_compare(.false.,msg(1:6),'STATUS'))  then
                option = 1
             else if (inp_compare(.false.,msg(1:7),'POSDATA'))  then
                option = 2
             else if (inp_compare(.false.,msg(1:8),'GETFORCE'))  then
                option = 3
             else if (inp_compare(.false.,msg(1:4),'INIT'))  then
                option = 4
             end if
         end if
         call ga_brdcst(9+MSGINT,option,1,MASTER)

         !*** EXIT ***
         if (option.eq.0)  then
            done = .true.

         !*** STATUS - writing to socket ***
         else if (option.eq.1) then
            if (taskid.eq.MASTER) then
               call util_talker_write(sock,statebuffer,12)
            end if

         !*** POSDATA - reading from socket ***
         else if (option.eq.2) then
            value = rtdb_parallel(.false.)
            if (taskid.eq.MASTER) then
               call util_talker_read(sock,unita,n1)
               call util_talker_read(sock,invunita,n1)
               call util_talker_read(sock,tmpint,4)
               nion0 = tmpint
               call util_talker_read(sock,dbl_mb(rion(1)),n2)
               call coords_set_bead_list(bead_list,1,dbl_mb(rion(1)))
            end if
            value = rtdb_parallel(.true.)
            !call coords_set_bead_list(bead_list,1,dbl_mb(rion(1)))

            !**** compute things here ****
            call socket_current_second(cpu1)
            call  run_bead_list(bead_list,1,task_gradient)
            energy = energy_bead_list(bead_list,1)
            call gradient_get_bead_list(bead_list,1,dbl_mb(fion(1)))
            call dscal(3*nion,-1.0d0,dbl_mb(fion(1)),1)
            call dcopy(9,0.0d0,0,stress,1)
            call socket_current_second(cpu2)

            iter = iter + 1
            if (taskid.eq.MASTER) then
               write(luout,*)
               write(luout,*) "== i-PI Socket Client Computation =="
               write(luout,'(" Connected to    = ",A)') socket_ip
               write(luout,'(" Number of atoms =",I8, 
     >                       " (natoms last read =",I8,")")') nion,nion0
               write(luout,'(" Iteration       =",I8)') iter
               write(luout,'(" Iteration Time  =",F8.3," seconds")')
     >            cpu2-cpu1
               write(luout,'(" Current Energy  =",F16.6)') energy
               write(luout,*)
            end if

            statebuffer(1:13) = "HAVEDATA     "

         !*** GETFORCE - writing to socket ***
         else if (option.eq.3) then
            statebuffer(1:13) = "FORCEREADY   "
            if (taskid.eq.MASTER) then
               call util_talker_write(sock,statebuffer,12)
               call util_talker_write(sock,energy,8)
               tmpint = nion
               call util_talker_write(sock,tmpint,4)
               call util_talker_write(sock,dbl_mb(fion(1)),n2)
               call util_talker_write(sock,stress,n1)
               tmpint = 1
               call util_talker_write(sock,tmpint,4)
               call util_talker_write(sock,tmpint,1)
            end if
            statebuffer(1:13) = "NEEDINIT     "

         !*** INIT - reading from socket ***
         else if (option.eq.4) then
            if (taskid.eq.MASTER) then
               call util_talker_read(sock,tmpint,4)
               call util_talker_read(sock,tmpint,4)
               nbytes = tmpint
               call util_talker_read(sock,buffer,nbytes)
            end if
            statebuffer(1:13) = "READY        "
         end if
      end do
      if (taskid.eq.MASTER) then
         call util_talker_close(sock)
      end if


      call delete_bead_list(bead_list,1)

      !**** delete from stack ****
      value =           MA_pop_stack(fion(2))
      value = value.and.MA_pop_stack(rion(2))
      if (.not.value)
     >   call errquit('socket_driver:stack failed',2,RTDB_ERR)


      socket_driver = value
      return
      end


      SUBROUTINE socket_current_second(T)
      real*8 T
      real*4 dummy,etime
      real*4 s(2)
      real   t2(2)
         
      dummy = etime(s)
      T = dble(s(1))

      RETURN
      END
