      subroutine grid_tasks(rtdb, coord, qsh, nqlen,
     &                      wqbuf, npts_per_task,
     &                      aqbuf, npts_per_ang)
c
c$Id: grid_tasks.F,v 1.6 1999-06-14 22:39:27 d3h449 Exp $
c
      implicit none
c
c     npts_per_task = target grid points per task
c
      integer rtdb, qsh, nqlen, npts_per_task, npts_per_ang
      integer lrij, irij, lrq, irq, lp, ip, lexpo, iexpo
      integer ndbl
      double precision coord, wqbuf, aqbuf, rad
      dimension coord(3,ncenters)
      dimension qsh(3,nqlen)
      dimension wqbuf(4,npts_per_task)
      dimension aqbuf(4,npts_per_ang)
c
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
      logical oreadfile, owritefile
      integer filesize, memsize
      character*255 gridfilename
      character*4 whichf
c
      integer iqsh, irsh, ictr, iang, npts_iqsh
      integer nq,n1,n2,nproc
Clrchi,irchi,lrdelchi,irdelchi,
c     ,     lrchi_loc,irchi_loc,lrde_loc,irde_loc
c      integer nest_tsk,nrchi_task
      integer nxtask
      external nxtask
      logical leb,lssw
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
c      nrchi_task=1
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
c
c     get partitioning function
c
      lssw=.true.
      whichf='dumm'
      if (.not.rtdb_cget(rtdb, 'dft:ssweights',  1, whichf))then
        lssw=.false.
      endif
      if(lssw) call xc_interp(whichf)
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
c     estimate of ntask/node
c
c      nproc=ga_nnodes()
cc
cc     is this too conservative or WRONG?
cc
c      nest_tsk=ncenters*1.2*nqshells/nproc
cc
cc     allocate rchi
cc
c      if (.not.MA_alloc_Get(MT_Dbl,nest_tsk,'rchi_atom',
c     &     lrchi,irchi))
c     &     call errquit('dft_scf: cannot allocate rchi_atom',0)
cc
c      if (.not.MA_alloc_Get(MT_Dbl,nest_tsk,'rdchi_atom',
c     &     lrdelchi,irdelchi))
c     &     call errquit('dft_scf: cannot allocate rchi_atom',0)
c      if (.not.MA_Push_Get(MT_Dbl,ncenters,'rchi_atom',
c     &     lrchi_loc,irchi_loc))
c     &     call errquit('dft_scf: cannot allocate rchi_atom',0)
cc
c      if (.not.MA_Push_Get(MT_Dbl,ncenters,'rdchi_atom',
c     &     lrde_loc,irde_loc))
c     &     call errquit('dft_scf: cannot allocate rchi_atom',0)
c      call dfill(nest_tsk, 0.0d0, dbl_mb(irchi), 1)
c      call dfill(nest_tsk, 0.0d0, dbl_mb(irdelchi), 1)
c      
      nproc=ga_nnodes()
c
      n1 = 0
      n2 = nxtask(nproc,1)
      do 100 iqsh = 1, nqshells
c
         if (n1.eq.n2 )then

         irsh = qsh(1,iqsh)         
         ictr = qsh(2,iqsh)
         iang = qsh(3,iqsh)
         nq = npts_per_ang
c
c        Generate all angular quadrature points for this radial shell
c
         call dfill(4*nq, 0.0d0, aqbuf, 1)
c
         call grid_qpts(coord(1,ictr), ictr, irsh, iang, aqbuf, nq, leb,
     ,        rad)
         npts_iqsh = nq
c
c        modify the weight to include the spatial weight (and symmetry factor)
c
         call grid_setspac(rtdb, coord, aqbuf, dbl_mb(irij), 
     &                     dbl_mb(irq), dbl_mb(ip), ictr, npts_iqsh,
     ,        lssw, whichf, dbl_mb(iexpo), rad)

         if (npts_iqsh.eq.0) then
            n1 = n1 + 1
            n2 = nxtask(nproc,1)
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
cc
cc     compute rchi and store it in memory
cc
cc coord redefined in setspac?
cc     WARNING: grid_rchi reorders aqbuf from (4,*) to (3,*)
cc
c         call grid_rchi(nrchi_task,npts_per_task,
c     ,        npts_iqsh,ictr, ncenters,
c     ,        dbl_mb(irchi),dbl_mb(irdelchi), 
c     ,        dbl_mb(irchi_loc),dbl_mb(irde_loc), 
c     ,        rad, aqbuf, coord, 
c     ,        dbl_mb(iexpo), dbl_mb(irq),dbl_mb(irij))
c
c
         n1 = n1 + 1
         n2 = nxtask(nproc,1)
c     
      else
c     
         n1 = n1 + 1
c     
      endif
  100 continue

      n1 = nxtask(-nproc,1)
c
c     Write last buffer of grid pts and close file     
c
      nq = 1
      n_tot_tasks = 0 ! sets nctr to 0; writes last buf and returns #tasks
      call grid_file_write(nq, n_tot_tasks, rad, aqbuf)
c      call grid_rchi(nrchi_task,npts_per_task,
c     ,     0,ictr, ncenters,
c     ,     dbl_mb(irchi),dbl_mb(irdelchi), 
c     ,     dbl_mb(irchi_loc),dbl_mb(irde_loc), 
c     ,     rad, aqbuf, coord, 
c     ,     dbl_mb(iexpo), dbl_mb(irq),dbl_mb(irij))
      write(6,*)' n_tot_tasks = ', n_tot_tasks
c  
c      if (.not. grid_file_close(.true.))
c     $   call errquit('grid_tasks: closing gridpts?', 0)
c
c
c     clean up
c
c      if (.not. MA_Pop_Stack(lrde_loc))
c     &   call errquit('grid_tasks: pop stack failed.',0)
c      if (.not. MA_Pop_Stack(lrchi_loc))
c     &   call errquit('grid_tasks: pop stack failed.',0)
cc      if (.not. MA_Pop_Stack(lrdelchi))
cc     &   call errquit('grid_tasks: pop stack failed.',0)
cc      if (.not. MA_Pop_Stack(lrchi))
cc     &   call errquit('grid_tasks: pop stack failed.',0)
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
