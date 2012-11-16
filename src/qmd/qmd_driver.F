c
c     main qmd driver
c
      logical function qmd_driver(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "stdio.fh"
c
#include "qmd_common.fh"
c
      integer rtdb
c
      integer geom                  ! geometry object
c
      integer l_r, k_r              ! coordinates
      double precision r(3)
c
      integer l_v, k_v              ! velocities
      double precision v(3)
c
      integer l_g, k_g              ! gradients
      integer l_m, k_m              ! masses
c
      integer nat                   ! number of atoms
      integer nsteps_nucl           ! total nuclear steps
      integer isteps_nucl           ! nuclear md step
      double precision dt_nucl      ! nuclear md time step
      double precision esys         ! system energy (from task_gradient)
      double precision ekin         ! kinetic energy from the md part
c
      logical status
c
      character*32 pname
c
      logical task_gradient
      external task_gradient
c
      status = .true.
      pname = "qmd_driver: "
c
c     initialize parameters
      call qmd_init(rtdb,nsteps_nucl,dt_nucl)
c
c     create geometry object 
      if (.not. geom_create(geom,'geometry'))
     &   call errquit(pname//'geom_create', 911, GEOM_ERR)
c
c     load in the geometry
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     &   call errquit(pname//'geom_rtdb_load', geom, RTDB_ERR)
c
c     get number of atoms
      if (.not. geom_ncent(geom,nat))
     $   call errquit(pname//'geom_ncent',geom, GEOM_ERR)
      if (nat.eq.0) return
c
c     local memory 
      if(.not.ma_push_get(mt_dbl,nat,'mass',l_m,k_m))  ! mass
     $   call errquit(pname//'ma_push_get mass', nat,MA_ERR)
      if(.not.ma_push_get(mt_dbl,3*nat,'coords',l_r,k_r)) ! coords
     $   call errquit(pname//'ma_push_get coords', 3*nat,MA_ERR)
      if(.not.ma_push_get(mt_dbl,3*nat,'vels',l_v,k_v))  ! vels
     $   call errquit(pname//'ma_push_get vels', 3*nat,MA_ERR)
      if(.not.ma_push_get(mt_dbl,3*nat,'grads',l_g,k_g))  ! grads
     $   call errquit(pname//'ma_push_get grads',3*nat,MA_ERR)
c
c     get the masses, initial coordinates and velocities
      if (.not.geom_masses_get(geom,nat,dbl_mb(k_m)))  ! masses
     $   call errquit(pname//'geom_masses_get',geom,GEOM_ERR)
      if (.not.geom_coords_vels_get(geom,dbl_mb(k_r),dbl_mb(k_v)))
     $   call errquit(pname//'geom_coords_vels_get',0,GEOM_ERR)
c
c     perform energy gradient calculation on initial coordinates
      esys = 0.d0
      if (.not. task_gradient(rtdb))
     &   call errquit(pname//'task_gradient',0,CALC_ERR)
      if(.not.rtdb_get(rtdb,'task:gradient',mt_dbl,3*nat,dbl_mb(k_g)))
     &   call errquit(pname//'rtdb_get gradient',0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,'task:energy',mt_dbl,1,esys))
     &   call errquit(pname//'rtdb_get energy',0,RTDB_ERR)
c
c     nuclear molecular dynamics loop
      ndof = 3*nat - 6.d0  ! 3N-6  need to determine this...
      do isteps_nucl = 1,nsteps_nucl
c
        call qmd_nucl_step(rtdb,geom,nat,dt_nucl,dbl_mb(k_m),
     &     dbl_mb(k_r),dbl_mb(k_v),dbl_mb(k_g),esys,ekin)   ! nuclear md step
c
        if (ga_nodeid().eq.0) write(luout,*) "ekin: ",ekin
        if (ga_nodeid().eq.0) write(luout,*) "esys: ",esys
        if (ga_nodeid().eq.0) write(luout,*) "etotal: ",ekin+esys
        if (ga_nodeid().eq.0) write(luout,*) "targ_temp: ",targ_temp
        if (ga_nodeid().eq.0) write(luout,*) "final_temp: ",final_temp
c
      end do  ! isteps_nucl
c
c     clear memory
      if(.not.ma_pop_stack(l_g))
     + call errquit(pname//'ma_pop_stack g',0,MA_ERR)
      if(.not.ma_pop_stack(l_v))
     + call errquit(pname//'ma_pop_stack v',0,MA_ERR)
      if(.not.ma_pop_stack(l_r))
     + call errquit(pname//'ma_pop_stack r',0,MA_ERR)
      if(.not.ma_pop_stack(l_m))
     + call errquit(pname//'ma_pop_stack m',0,MA_ERR)
      if(.not.ma_verify_allocator_stuff())
     + call errquit(pname//'ma_verify_allocator_stuff',0,MA_ERR)
      if(.not.geom_destroy(geom))
     + call errquit(pname//'geom_destroy',0,GEOM_ERR)
c
      qmd_driver = status
c
      end
