c
c     $Id: paulsmatrix.F 27778 2015-10-22 21:25:19Z bylaska $
c

*     *********************************
*     *                               *
*     *       task_paulsmatrix        *
*     *                               *
*     *********************************
      logical function task_paulsmatrix(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"
#include "bafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)

      logical value,oprint
      character*80 filename
      character*255 full_filename,full_filename2
      integer norder,l,msglen
      integer tmp(2)
      integer i,j,ilo,ihi,jlo,jhi
      real*8 cpu1,cpu2,cpu3,cpu4,cpu5,diag_time,read_time

*     
      value = .true.
      taskid = ga_nodeid()
      np     = ga_nnodes()
      oprint = (taskid.eq.MASTER)

      if (btdb_cget(rtdb,'paulsmatrix_filename',1,filename)) then
         call current_second(cpu1)
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                                full_filename)
         if (oprint) then
            write(luout,*) "np = ",np
            write(luout,*) "matrix filename:      ", filename
            write(luout,*) "matrix full_filename: ", full_filename
         end if


         if (taskid.eq.MASTER) then
            l = index(full_filename,' ') - 1
            call openfile(5,full_filename,l,'r',l)
            call iread(5,norder,1)
            call closefile(5)
         end if
         if (np.gt.1) then
            msglen = 1
            call ga_brdcst(9+MSGINT,norder,mitob(msglen),MASTER)
         end if

         if (oprint) then
             write(luout,*) "matrix order = ",norder
             write(luout,*) 
         end if

         !*** allocate tmp array ***
         if (.not.BA_alloc_get(mt_dbl,(norder),'tmp',tmp(2),tmp(1)))
     >      call errquit("task_paulsmatrix:Cannot create tmp matrix",
     >                   0,MA_ERR)


        call current_second(cpu2)
c        call paulsmatrix_george(full_filename,
c     >                          norder,dbl_mb(tmp(1)),
c     >                          read_time,diag_time)

        call paulsmatrix_eric(full_filename,
     >                          norder,dbl_mb(tmp(1)),
     >                          read_time,diag_time)
        call current_second(cpu3)

        if (oprint) then
           write(luout,*) 
           write(luout,*) "10 lowest Eigenvalues:"
           do j=1,10
              write(luout,101) j,dbl_mb(tmp(1)+j-1)
           end do
        end if
 101  format(" Eigenvalue(",I6,")=",F20.6)



         if (.not.BA_free_heap(tmp(2)))
     >      call errquit("task_paulsmatrix:Cannot destroy tmp matrix",
     >                   0,MA_ERR)

        call current_second(cpu4)
        if (taskid.eq.MASTER) then
           write(luout,*) 
           write(luout,*) "Timings:"
           write(luout,*) "---------------------------"
           write(luout,'(A,E15.6)') " prologue:   ",cpu2-cpu1
           write(luout,'(A,E15.6)') " read matrix:",read_time
           write(luout,'(A,E15.6)') " diagonalize:",diag_time
           write(luout,'(A,E15.6)') " epilogue:   ",cpu4-cpu3
        end if

      else
         if (oprint) then
            write(luout,*) 
     >      "Filename not entered: set paulsmatrix_filename datfile"
         end if
      end if

      task_paulsmatrix = value
      return
      end


      subroutine paulsmatrix_george(full_filename,norder,eigs,
     >                              read_time,diag_time)
      implicit none
      character*(*) full_filename
      integer norder
      real*8  eigs(*)
      real*8 read_time,diag_time

#include "btdb.fh"
#include "bafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)

      logical value,oprint
      integer l,msglen
      integer hmatrix,vmatrix
      integer i,j,ilo,ihi,jlo,jhi
      real*8 cpu1,cpu2,cpu3

      value = .true.
      taskid = ga_nodeid()
      np     = ga_nnodes()
      oprint = (taskid.eq.MASTER)

      !*** create the ga matrix hmatrix ***
      if (.not.ga_create(mt_dbl,norder,norder,'hmatrix1',
     >               norder,0,hmatrix))
     >   call errquit("paulsmatrix_george:Cannot create ga matrix",
     >                0,GA_ERR)

      !*** create the ga matrix vmatrix ***
      if (.not.ga_create(mt_dbl,norder,norder,'vmatrix1',
     >               norder,0,vmatrix))
     >   call errquit("paulsmatrix_george:Cannot create ga matrix",
     >                1,GA_ERR)

      call current_second(cpu1)
      if (oprint) 
     > write(luout,*) "Into Reading matrix"

      if (taskid.eq.MASTER) then
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,norder,1)
         do j=1,norder
            call dread(5,eigs,norder)
            ilo = 1
            ihi = norder
            jlo = j
            jhi = j
            call ga_put(hmatrix,ilo,ihi,jlo,jhi,eigs,norder)
         end do
         call closefile(5)
      end if
      call current_second(cpu2)

      if (oprint)
     >write(luout,*) "Into diagonalize matrix,current time=",cpu2-cpu1

      call ga_diag_std(hmatrix,vmatrix,eigs)
      call current_second(cpu3)

        !*** destroy the ga matrix hmatrix ***
        if (.not.ga_destroy(hmatrix))
     >      call errquit("paulsmatrix_george:Cannot destroy ga matrix",
     >                   0,GA_ERR)

        !*** destroy the ga matrix vmatrix ***
        if (.not.ga_destroy(vmatrix))
     >      call errquit("paulsmatrix_george:Cannot destroy ga matrix",
     >                   1,GA_ERR)

      read_time = cpu2-cpu1
      diag_time = cpu3-cpu2
      return
      end


      subroutine paulsmatrix_eric(full_filename,norder,eigs,
     >                            read_time,diag_time)
      implicit none
      character*(*) full_filename
      integer norder
      real*8  eigs(*)
      real*8 read_time,diag_time

#include "btdb.fh"
#include "bafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)

      logical oprint
      real*8 cpu1,cpu2,cpu3,cpu4
      integer i,j,np_i,np_j,mapping,l
      integer hml(2)

*     **** external functions ****
      logical  Dneall_m_allocate,Dneall_m_free
      external Dneall_m_allocate,Dneall_m_free


      call current_second(cpu1)
      call dcopy(norder,0.0d0,0,eigs,1)

*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      np_i = np
      np_j = 1
      do j=1,np
         do i=1,np
            if ((i*j).eq.np) then
               if (abs(i-j).lt.abs(np_i-np_j)) then
                  np_i = i
                  np_j = j
               end if
            end if
         end do
      end do

      call Parallel2d_Init(np_j)
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

      mapping = 1                       !0-cyclic, 1-block mapping
      call Dne_init(1,norder,mapping)  


      if (oprint) then
         write(luout,*)
         write(luout,*) "Erics Eigenvalue Solver:"
         write(luout,*) "np   = ",np
         write(luout,*) "np_i = ",np_i
         write(luout,*) "np_j = ",np_j
         write(luout,*)
      end if
      if (.not.Dneall_m_allocate(1,hml))
     > call errquit("paulsmatrix_eric:Cannot allocate hml",0,MA_ERR)


      if (oprint) write(luout,*) "Into Reading matrix"

      if (taskid.eq.MASTER) then
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,norder,1)
      end if

      do j=1,norder
         if (taskid.eq.MASTER) call dread(5,eigs,norder)
         call Parallel_Brdcst_values(0,norder,eigs)
         do i = 1,norder
            call Dneall_m_set_value(eigs(i),1,1,i,j,dbl_mb(hml(1)))
         end do
      end do

      if (taskid.eq.MASTER) then
         call closefile(5)
         write(*,*) "Into diagonalize"
      end if
      call current_second(cpu2)
      call Dneall_m_diagonalize(1,dbl_mb(hml(1)),eigs,.true.)

      call current_second(cpu3)

      if (.not.Dneall_m_free(hml))
     > call errquit("paulsmatrix_eric:Cannot free hml",0,MA_ERR)

      call Dne_end()
      call Parallel2d_Finalize()
      call Parallel_Finalize()

      read_time = cpu2-cpu1
      diag_time = cpu3-cpu2

      return
      end
