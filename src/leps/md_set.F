      logical function md_set(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
c      external md_driver
c
      integer geom,natom 
      integer md_data_get_i_c, md_data_get_i_m, md_data_get_i_q
      integer md_data_get_i_t, md_data_get_i_v, md_data_get_i_a
      integer md_data_get_i_b, md_data_get_i_d, md_data_get_i_f
      integer md_data_get_natom

      integer coor,mass,bq,tag,vel,acc,bterm,dterm,frc

      integer i 

c      set input parameters 
c      write(*,*) ' In MD_SET '
      call md_data_inp_set(rtdb)
c      write(*,*) ' In MD_SET -1'
cc     create geometry object 
      if (.not. geom_create(geom, 'geometry'))
     & call errquit('md_set_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
cC     load geometry into the object
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     & call errquit('md_set_getgeom: no geometry',
     & 0,RTDB_ERR)

cc    Set the natom, it is needed for all other allocations
      natom=md_data_get_natom()

C     Get the handles for everything
      coor=md_data_get_i_c()
      mass=md_data_get_i_m()
      bq=md_data_get_i_q()
      tag=md_data_get_i_t()
      vel=md_data_get_i_v()
      acc=md_data_get_i_a()
      bterm=md_data_get_i_b()
      dterm=md_data_get_i_d()
      frc=md_data_get_i_f()
c
c
cC ----- Fill up the local copies with data(mass,coords,charge,tags) from geometry object
cC------- accessed by the geom handle
      call md_data_geom_copy(rtdb,geom) 

c      write(*,*)  coor,mass,bq,tag,vel,acc,bterm,dterm,frc   
       
C ----- get initial velocities and accelerations and other inital conditions
c       if (.not. grad(rtdb)) 
c     & call errquit('md_gradient: grad not calculated', 0, RTDB_ERR)
c  ---- get the gradient to fill up the force array        
c      if (.not. rtdb_get(rtdb,theory,mt_dbl,3*natom,dbl_mb(frc)))
c     & call errquit('md_gradient:  no grad found', 0, RTDB_ERR)
      do i=1,3*natom
         dbl_mb(coor+i-1)=0.0
         dbl_mb(vel+i-1)=0.0
         dbl_mb(acc+i-1)=0.0
         dbl_mb(bterm+i-1)=0.0
         dbl_mb(dterm+i-1)=0.0
      enddo

      if (.not. geom_destroy(geom))
     & call errquit('md_set:  geom_destroy failed',
     & 911,GEOM_ERR)
c      write(*,*) 'OUt of MD_set'
      md_set=.true.
      end



c $Id$
