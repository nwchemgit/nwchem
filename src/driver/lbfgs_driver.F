      logical function lbfgs_driver(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
      integer irtdb
      double precision eps,xtol,gtol,stpmin,stpmax
      integer iprint(2),iflag,icall,n,m,mp,lp,j,iter
      logical diagco
      integer master
c
c     the driver for LBFGS must always declare LB2 as EXTERNAL
c
      external lb2
      common /lb3/mp,lp,gtol,stpmin,stpmax
      common /lb3p/master

      integer igeom
      integer ncent
      integer i_c,h_c
      integer i_g,h_g
      integer i_scr,h_scr
      integer i_d,h_d
      integer nscr
      integer maxiter,nrest
      double precision energy,stp1
      logical task_gradient
      external task_gradient
      character*32 pname

      pname = "lbfgs_driver"

      master=ga_nodeid()
      m=5
      diagco= .false.
      eps= 1.0d-5
      xtol= 1.0d-16
      icall=0

      if(.not.geom_create(igeom,'geometry'))
     + call errquit(pname//'Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit(pname//'Failed to load geometry',0, GEOM_ERR)

      if(.not. geom_ncent(igeom,ncent) ) 
     >    call errquit("qmmm:geom_ncent",0,0)

      if(.not.ma_alloc_get(mt_dbl,3*ncent,'c',h_c,i_c))
     + call errquit( pname//'Failed to allocate memory for c',
     + 3*ncent, MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,3*ncent,'g',h_g,i_g))
     + call errquit( pname//'Failed to allocate memory for c',
     + 3*ncent, MA_ERR)

      nscr = 3*ncent*(2*m+1)+2*m
      if(.not.ma_alloc_get(mt_dbl,nscr,'scratch',h_scr,i_scr))
     + call errquit( pname//'Failed to allocate memory for scr',
     + 3*ncent, MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,3*ncent,'diag',h_d,i_d))
     + call errquit( pname//'Failed to allocate memory for scr',
     + 3*ncent, MA_ERR)

       if(.not.geom_cart_coords_get(igeom,dbl_mb(i_c)))
     + call errquit('qmmm: Failed to get coord ',0, GEOM_ERR)

      if(.not.rtdb_get(irtdb,'lbfgs:maxiter',mt_int,1,maxiter))
     +  maxiter = 20   

c
      stp1 = 0.0d0
10    continue
      n=100
      iprint(1)= 1
      iprint(2)= 0
      iflag=0
c
      lbfgs_driver = .true.
 20   continue
      if (.not. task_gradient(irtdb))
     $   call errquit('driver: task_gradient failed',0, GEOM_ERR)
      if (.not. rtdb_get(irtdb,'task:energy', mt_dbl, 1, energy))
     $   call errquit('driver: could not get energy',0, RTDB_ERR)
      if (.not. rtdb_get(irtdb, 'task:gradient', mt_dbl, 3*ncent,
     $  dbl_mb(i_g))) call errquit('driver: could not get gradient',0,0)

      write(*,*) "stp1, energy",stp1, energy
      call lbfgs(3*ncent,m,dbl_mb(i_c),energy,
     >           dbl_mb(i_g),diagco,
     >           dbl_mb(i_d),iprint,eps,xtol,
     >           dbl_mb(i_scr),iflag,iter,stp1)

c
c     if line search exited due to may function evaluations
c     try restarting 
      if(iflag.eq.-1) then
       if(master.eq.0) then
          write(6,*) "@ too many function evaluations"
          write(6,*) "@ restarting the optimization"
          call util_flush(6) 
       end if
c       iflag=0
c       call dfill(nscr,0.0d0,dbl_mb(i_scr),1)
c       stp1=1.0
       go to 20
      end if
      if(iflag.le.0) then
          lbfgs_driver=.false.
          go to 50
      end if
      icall=icall + 1
      if(iter.gt.maxiter) go to 50
      if(.not.geom_cart_coords_set(igeom,dbl_mb(i_c)))
     + call errquit('qmmm: Failed to get coord ',0, GEOM_ERR)
      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit(pname//'Failed to load geometry',0, GEOM_ERR)
      go to 20
  50  continue

      if(.not.ma_free_heap(h_d))
     & call errquit(pname//'Failed to deallocate stack c',ncent,
     &       MA_ERR)


      if(.not.ma_free_heap(h_scr))
     & call errquit(pname//'Failed to deallocate stack c',ncent,
     &       MA_ERR)

      if(.not.ma_free_heap(h_g))
     & call errquit(pname//'Failed to deallocate stack c',ncent,
     &       MA_ERR)


      if(.not.ma_free_heap(h_c))
     & call errquit(pname//'Failed to deallocate stack c',ncent,
     &       MA_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      end
c
c     ** last line of simple driver (sdrive) **

c $Id$
