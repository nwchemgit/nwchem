C$Id: util_sgroup.F,v 1.1 2004-11-29 21:37:17 edo Exp $
      logical function util_sgstart()
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "util_sgroup.fh"
      integer nprocs, nsize, midproc, proclist(1000)
      integer i
c
c    Create groups
c   
      write(luout,*) 'Creating groups'
      nprocs = ga_nnodes()
      midproc = nprocs/2
c
c     hardwired value of 2 subgroups
c
      ngroups=2
      do i = 1, midproc
        proclist(i) = i - 1
      end do
      nsize = midproc
      if (ga_nodeid().eq.0) then
        do i = 1, nsize
          write(luout,*) 'Group 1 ',i,proclist(i)
        end do
      endif
      grouplist(1) = ga_pgroup_create(proclist,nsize)
      do i  = midproc + 1, nprocs
        proclist(i-midproc) = i - 1
      end do
      nsize = nprocs - midproc
      if (ga_nodeid().eq.0) then
        do i = 1, nsize
          write(luout,*) 'Group 2 ',i,proclist(i)
        end do
      endif
      grouplist(2) = ga_pgroup_create(proclist,nsize)
      write(luout,*) 'Groups are created'
      if (ga_nodeid().lt.midproc) then
        call ga_pgroup_set_default(grouplist(1))
c
c     active
c
        util_sgstart=.true.
        mygroup=1
        write(0,*) ga_nodeid(),' active'
      else
        call ga_pgroup_set_default(grouplist(2))
c
c     inactive
c        
        util_sgstart=.false.
        mygroup=2
        write(0,*) ga_nodeid(),' inactive'
      endif
      sgactive=util_sgstart
      return
      end
      subroutine util_sgend(rtdb,rtdb_name)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
      integer rtdb
      character*(*) rtdb_name
c
      if(sgactive) then
         write(0,*) ga_nodeid(),' deleting cloned rtdb '
         if (.not. rtdb_close(rtdb, 'delete'))call errquit
     (        (' sgend: rtdb_close and delete failed ',0,0)
      endif
      write(luout,*) 'end of subgroups call '
      call ga_pgroup_set_default(ga_pgroup_get_world())
c
c     reopen world rtdb
c
      if (.not. rtdb_open(rtdb_name, 'old', rtdb))
     $     call errquit('sgend: rtdb_open old failed ', 0, 0)
      
      return
      end
      subroutine util_sgrtdb(rtdb,rtdbname)
      implicit none
      integer rtdb
      character*(*) rtdbname
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
      character*256 grtdbname
      integer ii
c
      write(luout,*) 'util_sgrtdb start'
      call util_pname0(rtdbname,grtdbname,ii,ngroups)
c     
c     creates rtdb file grtdbname=rtdbname.ii
c     
      if(.not.rtdb_clone(rtdb,grtdbname))
     .     call errquit('oo ',0,0)
      call ga_sync()
      if (.not. rtdb_open(grtdbname, 'old', rtdb)) call errquit
     (     (' rtdb_open old failed ',0,0)
      return
      end
