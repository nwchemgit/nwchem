      subroutine grid_tasks(rtdb, coord, qsh, nqlen,
     &                      aqbuf, npts_per_ang,nq_tot)
c
c$Id: grid_tasks.F,v 1.15 2000-02-02 02:37:31 edo Exp $
c
      implicit none
c
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
      integer rtdb, qsh, nqlen, npts_per_ang, nq_tot
      integer lrij, irij, lrq, irq, lp, ip, lexpo, iexpo
      integer ndbl
      double precision coord, aqbuf, rad
      dimension coord(3,ncenters)
      dimension qsh(3,nqlen)
      dimension aqbuf(4,npts_per_ang)
c
      logical oreadfile, owritefile
      integer filesize, memsize
      character*255 gridfilename
      character*4 whichf
c
      integer iqsh, irsh, ictr, iang, npts_iqsh
      integer nq,nproc,me,g_buck,g_buckind
      logical trunc_neigh
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
      owritefile= grid_test_mem (rtdb, geom, nq_tot)
c
c
c     get partitioning function
c

      if(lssw) do_interp = xc_interp(whichf)
c
c     Allocate some scratch space for construction of the spatial weights
c
      ndbl = (ncenters*(ncenters+1))/2
      if (.not.MA_Push_Get(MT_Dbl,ndbl,'Rij',lRij,iRij))
     &   call errquit('  in grid_tasks: cannot allocate Rij',0)
C      call a_dist(coord, Dbl_MB(iRij), ncenters)
      if (.not.MA_Push_get(mt_dbl,ncenters*npts_per_ang,'rq',lrq,irq))
     &   call errquit('grid_tasks: cannot allocate rq',0)
      if (.not.MA_Push_get(mt_dbl,ncenters*npts_per_ang,'p',lp,ip))
     &   call errquit('grid_tasks: cannot allocate rq',0)
      if (.not.MA_Push_get(mt_dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_tasks: cannot allocate rq',4)
c      
      nproc=ga_nnodes()
      me=ga_nodeid()
c
      if (.not.rtdb_get(rtdb, 'dft:trunc_neigh', mt_log, 1,
     &   trunc_neigh))trunc_neigh=.true.

      call grid_buckinit(g_buck,g_buckind)

      do 100 iqsh = 1+me, nqshells,nproc
c

         irsh = qsh(1,iqsh)         
         ictr = qsh(2,iqsh)
         iang = qsh(3,iqsh)
         nq = npts_per_ang
c
c        Generate all angular quadrature points for this radial shell
c
         call dfill(4*nq, 0.0d0, aqbuf, 1)
c
         call grid_qpts(coord(1,ictr), ictr, irsh, iang, aqbuf, nq, 
     ,        rad)
         npts_iqsh = nq
c
c        modify the weight to include the spatial weight (and symmetry factor)
c
         call grid_setspac(rtdb, coord, aqbuf, dbl_mb(irij), 
     &                     dbl_mb(irq), dbl_mb(ip), ictr, npts_iqsh,
     ,        dbl_mb(iexpo), rad, trunc_neigh)

         if (npts_iqsh.eq.0) then
            go to 100
         endif
         
c
c        since only generating symmetry unique quadrature 
c        scale weights with symmetry factor if desired
c
         call grid_sym(ictr, aqbuf, npts_iqsh)
c            
c        now pack this aqbuf into the write buffer (write and update task counter as needed)
c
         call grid_file_write(npts_iqsh, ictr, rad, aqbuf)
c

  100 continue

      call grid_buckend(g_buck,g_buckind)

c
c     Write last buffer of grid pts and close file     
c
      nq = 1
      n_tot_tasks = 0 ! sets nctr to 0; writes last buf and returns #tasks
      call grid_file_write(nq, n_tot_tasks, rad, aqbuf)
c      write(6,*)' n_tot_tasks = ', n_tot_tasks
c
c     clean up
c
      if (.not. MA_Pop_Stack(lexpo))
     &   call errquit('grid_tasks: pop stack failed.',1)
      if (.not. MA_Pop_Stack(lp))
     &   call errquit('grid_tasks: pop stack failed.',2)
      if (.not. MA_Pop_Stack(lrq))
     &   call errquit('grid_tasks: pop stack failed.',3)
      if (.not. MA_Pop_Stack(lrij))
     &   call errquit('grid_tasks: pop stack failed.',4)
      return
      end 
