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
      integer hmatrix,vmatrix,tmp(2)
      integer i,j,ilo,ihi,jlo,jhi
      real*8 cpu1,cpu2,cpu3,cpu4,cpu5

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
            !call closefile(5)
         end if
         if (np.gt.1) then
            msglen = 1
            call ga_brdcst(9+MSGINT,norder,mitob(msglen),MASTER)
         end if


         if (oprint) then
             write(luout,*) "matrix order = ",norder
         end if
         !write(luout,*) "MATRIX ORDER = ",norder,taskid

         !*** allocate tmp array ***
         if (.not.BA_alloc_get(mt_dbl,(norder),'tmp',tmp(2),tmp(1)))
     >      call errquit("task_paulsmatrix:Cannot create tmp matrix",
     >                   0,MA_ERR)

         !*** create the ga matrix hmatrix ***
         if (.not.ga_create(mt_dbl,norder,norder,'hmatrix1',
     >                  norder,0,hmatrix))
     >      call errquit("task_paulsmatrix:Cannot create ga matrix",
     >                   0,GA_ERR)

         !*** create the ga matrix vmatrix ***
         if (.not.ga_create(mt_dbl,norder,norder,'vmatrix1',
     >                  norder,0,vmatrix))
     >      call errquit("task_paulsmatrix:Cannot create ga matrix",
     >                   1,GA_ERR)

        call current_second(cpu2)
        if (taskid.eq.MASTER) then
        write(luout,*) "Reading matrix, time=",cpu2-cpu1
        do j=1,norder
            call dread(5,dbl_mb(tmp(1)),norder)

            ilo = 1
            ihi = norder
            jlo = j
            jhi = j
            call ga_put(hmatrix,ilo,ihi,jlo,jhi,dbl_mb(tmp(1)),norder)
        end do
        call closefile(5)
        end if

        call current_second(cpu3)
        if (taskid.eq.MASTER)
     >    write(luout,*) "Diagonalize matrix, time=",cpu3-cpu1
        call ga_diag_std(hmatrix,vmatrix,dbl_mb(tmp(1)))
        call current_second(cpu4)
        if (taskid.eq.MASTER)
     >    write(luout,*) "Finished Diagonalizing, time=",cpu4-cpu1
        if (taskid.eq.MASTER) then
           write(luout,*) 
           write(luout,*) "Eigenvalues:"
           do j=1,norder
              write(luout,101) j,dbl_mb(tmp(1)+j-1)
           end do
        end if
 101  format("EIG(",I6,")=",F20.6)


        !*** destroy the ga matrix hmatrix ***
        if (.not.ga_destroy(hmatrix))
     >      call errquit("task_paulsmatrix:Cannot destroy ga matrix",
     >                   0,GA_ERR)

        !*** destroy the ga matrix vmatrix ***
        if (.not.ga_destroy(vmatrix))
     >      call errquit("task_paulsmatrix:Cannot destroy ga matrix",
     >                   1,GA_ERR)

         if (.not.BA_free_heap(tmp(2)))
     >      call errquit("task_paulsmatrix:Cannot destroy tmp matrix",
     >                   0,MA_ERR)

        call current_second(cpu5)
        if (taskid.eq.MASTER) then
           write(luout,*) 
           write(luout,*) "Timings:"
           write(luout,*) "---------------------------"
           write(luout,'(A,E15.6)') " prologue:   ",cpu2-cpu1
           write(luout,'(A,E15.6)') " read matrix:",cpu3-cpu2
           write(luout,'(A,E15.6)') " diagonalize:",cpu4-cpu3
           write(luout,'(A,E15.6)') " epilogue:   ",cpu5-cpu4
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

