      subroutine grid_tasks(rtdb, coord, qsh, nqlen,
     &                      wqbuf, npts_per_task,
     &                      aqbuf, npts_per_ang)
c
c$Id: grid_tasks.F,v 1.1 1999-04-23 16:30:36 d3h449 Exp $
c
      implicit none
c
c     npts_per_task = target grid points per task
c
      integer rtdb, qsh, nqlen, npts_per_task, npts_per_ang
      double precision coord, wqbuf, aqbuf
      dimension coord(3,ncenters)
      dimension qsh(3,nqlen)
      dimension wqbuf(4,npts_per_task)
      dimension aqbuf(4,npts_per_ang)
c
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      logical oreadfile, owritefile
      integer filesize, memsize
      character*255 gridfilename
c
      integer iqsh, irsh, ictr, iang, npts_iqsh
      integer nq
      logical leb
      logical grid_file_open
      external grid_file_open
      logical grid_file_close
      external grid_file_close
      logical grid_test_mem
      external grid_test_mem
c
c     input   
c        qsh(1,iqsh) = irsh | where within nrad
c        qsh(2,iqsh) = ictr | current center
c        qsh(3,iqsh) = iang | angular quad for given irsh radial
c    output
c        iqtask(1,iqt) = points in task iqt
c        iqtask(2,iqt) = integration center for task iqt
c
c
c     Open grid pts file
c
      call grid_packing_info
c
c     temp set gridfilename to gridpts
c
c     Grid_pts caching/file information (from scf_get_info)
c
      if (.not. rtdb_get(rtdb, 'grid:filesize',
     $     MT_INT, 1, filesize)) filesize = 0
      if (.not. rtdb_get(rtdb, 'grid:memsize',
     $     MT_INT, 1, memsize)) memsize = 0
c
c     The opening routine will put the .pid on the grid_pts filename
c     (hence even tho' parallel file, open as sequential)
c
      if (.not. rtdb_cget(rtdb, 'grid:filename', 1, gridfilename)) 
     $   call util_file_name('gridpts',.true.,.false.,gridfilename)
      if (.not. rtdb_cput(rtdb, 'grid:filename', 1, gridfilename)) 
     $   call errquit('grid_tasks: trouble putting gridfilename', 0)
c
      oreadfile = .false.
      owritefile= .false.
c
      if (filesize.gt.0 .or. memsize.gt.0) then
         if (.not. grid_file_open(gridfilename, memsize, filesize))
     $        call errquit('grid_tasks: opening gridpts?', 0)
      endif
c
      owritefile= grid_test_mem (rtdb, geom)
c
      leb = .false.
      if (iangquad.ne.0) leb = .true.
      do 100 iqsh = 1, nqshells
c
         irsh = qsh(1,iqsh)         
         ictr = qsh(2,iqsh)
         iang = qsh(3,iqsh)
         nq = npts_per_ang
c
c        Generate all angular quadrature points for this radial shell
c
         call grid_qpts(coord(1,ictr), ictr, irsh, iang, aqbuf, nq, leb)
         npts_iqsh = nq
c
c        now pack this aqbuf into the write buffer (write and update task counter as needed)
c
         call grid_file_write(npts_iqsh, ictr, aqbuf)
c
  100 continue
c
c     Write last buffer of grid pts and close file     
c
      nq = 1
      n_tot_tasks = 0 ! sets nctr to 0; writes last buf and returns #tasks
      call grid_file_write(nq, n_tot_tasks, aqbuf)
      write(6,*)' n_tot_tasks = ', n_tot_tasks
c  
c      if (.not. grid_file_close(.true.))
c     $   call errquit('grid_tasks: closing gridpts?', 0)
c
      return
      end 
