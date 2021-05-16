         !call nwpw_talker("127.0.0.1",9,50001,np_i)
         !call nwpw_talker_write(np_i,"Eric says hello world! ",23)
         !call nwpw_talker_write(np_i,"hello world again ericdone",26)
         !call nwpw_talker_close(np_i)

      subroutine runsocket_ipi(nion,energy,rion,fion)
      implicit none
      integer nion
      real*8 energy
      real*8 rion(3,nion),fion(3,nion)

#include "inp.fh"
#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical done
      integer sock,nbytes,nion0,option,ii,it,port,inet
      integer*4 tmpint
      real*8 cpu1,cpu2

      character*80 statebuffer
      character*80 buffer,ipport,stype
      character*20 msg,ip,pp

*     **** external functions ****
      real*8   cgsd_energy
      external cgsd_energy

      real*8 unita(3,3),iunita(3,3),stress(3,3)
      real*8 rion0(3,20)

      call control_socket_type(stype)
      call control_socket_ip(ipport)
      if (inp_compare(.false.,stype,'unix')) then
         ip = ipport
         inet = 0
         ii = index(ipport,' ')
         port = 0
      else
         ii = index(ipport,':')
         it = index(ipport,' ') - 1
         pp(1:it-ii) = ipport(ii+1:it)
         pp(it-ii+1:it-ii+5) = '    '
         read(pp,*,err=300) port
         go to 301
 300     port = 31415
 301     continue
         ip(1:ii-1) = ipport(1:ii-1)
         ip(ii:ii+4) = "    "
      end if

      call Parallel_taskid(taskid)

      statebuffer(1:13) = "READY        "

      if (taskid.eq.MASTER) then
         call nwpw_talker(ip,inet,ii-1,port,sock)
         !call nwpw_talker("127.0.0.1",9,port,sock)
      end if
      done = .false.
      it = 0
      do while (.not.done)
         msg = "            "

         !**** read option ****
         option = 0
         if (taskid.eq.MASTER) then
             call nwpw_talker_read(sock,msg,12)
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
         call Parallel_Brdcst_ivalue(MASTER,option)

         !*** EXIT ***
         if (option.eq.0)  then
            done = .true.

         !*** STATUS - writing to socket ***
         else if (option.eq.1) then
            if (taskid.eq.MASTER) then
               call nwpw_talker_write(sock,statebuffer,12)
            end if

         !*** POSDATA - reading from socket ***
         else if (option.eq.2) then
            if (taskid.eq.MASTER) then
               call nwpw_talker_read(sock,unita,9*8)
               call nwpw_talker_read(sock,iunita,9*8)
               call nwpw_talker_read(sock,tmpint,4)
               nion0 = tmpint
               call nwpw_talker_read(sock,rion0,3*8*nion0)
            end if
            call Parallel_Brdcst_ivalue(MASTER,nion0)
            call Parallel_Brdcst_values(MASTER,9,unita)
            call Parallel_Brdcst_values(MASTER,9,iunita)
            call Parallel_Brdcst_values(MASTER,3*nion,rion0)

            do ii=1,nion0
               rion(1,ii) = rion0(1,ii)
               rion(2,ii) = rion0(2,ii)
               rion(3,ii) = rion0(3,ii)
            end do

            !**** compute things here ****
            call current_second(cpu1)
            call ewald_phafac()
            energy = cgsd_energy(.false.)
            call cgsd_energy_gradient(fion)
            !call dscal(3*nion,(-1.0d0),fion,1)
            call dcopy(9,0.0d0,0,stress,1)
            call current_second(cpu2)

            if (taskid.eq.MASTER) then
               it = it + 1
               write(*,*) "== i-PI Socket Client Computation =="
               write(*,'(" Connected to   = ",A)') ipport
               write(*,'(" Iteration      =",I8)') it
               write(*,'(" Iteration Time =",F8.3," seconds")')cpu2-cpu1
               write(*,'(" Current Energy = ",F16.6)') energy
c               write(*,*) "energy=",energy,nion
c               write(*,*) "rion,fion=" 
c               do ii=1,nion
c                  write(*,*) ii,rion(1,ii),rion(2,ii),rion(3,ii),
c     >                          fion(1,ii),fion(2,ii),fion(3,ii) 
c               end do
            end if

            statebuffer(1:13) = "HAVEDATA     "

         !*** GETFORCE - writing to socket ***
         else if (option.eq.3) then
            statebuffer(1:13) = "FORCEREADY   "
            if (taskid.eq.MASTER) then
               call nwpw_talker_write(sock,statebuffer,12)
               call nwpw_talker_write(sock,energy,8)
               tmpint = nion
               call nwpw_talker_write(sock,tmpint,4)
               call nwpw_talker_write(sock,fion,3*nion*8)
               call nwpw_talker_write(sock,stress,9*8)
               tmpint = 1
               call nwpw_talker_write(sock,tmpint,4)
               call nwpw_talker_write(sock,tmpint,1)
            end if
            statebuffer(1:13) = "NEEDINIT     "

         !*** INIT - reading from socket ***
         else if (option.eq.4) then
            if (taskid.eq.MASTER) then
               call nwpw_talker_read(sock,tmpint,4)
               call nwpw_talker_read(sock,tmpint,4)
               nbytes = tmpint
               call nwpw_talker_read(sock,buffer,nbytes)
            end if
            statebuffer(1:13) = "READY        "
         end if
      end do
      if (taskid.eq.MASTER) then
         call nwpw_talker_close(sock)
      end if

      return
      end
