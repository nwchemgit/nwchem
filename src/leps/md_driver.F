      function md_driver(rtdb)
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
      logical md_driver
c
c      integer geom
c      integer natom
c      integer h_c,i_c
c      integer h_m,i_m
c      integer h_q,i_q
c      integer h_t,i_t
c      double precision eleps, md_driver1
       md_driver=.true.
c
c      write(*,*) "doing leps energy"
c      eleps=0.0 
cc     create geometry object 
c      if (.not. geom_create(geom, 'geometry'))
c     & call errquit('leps_getgeom: failed creating geometry',
c     & 0,GEOM_ERR)
cC     load geometry into the object
c      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
c     & call errquit('leps_getgeom: no geometry',
c     & 0,RTDB_ERR)
cC     get the number of centers 
c      if (.not. geom_ncent(geom,natom))
c     & call errquit('leps_getgeom: ncent error',
c     & 0,GEOM_ERR)      
cc     print tgeometry
c      if (.not. geom_print(geom))
c     & call errquit('leps_getgeom: print error',
c     & 0,GEOM_ERR)      
cc    
cc     allocate memory for coordinates
c      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'coordinate' ,
c     & h_c,i_c))
c     & call errquit('leps_getgeom: coordinate not allocated',
c     & 3*natom,MA_ERR)
c      if (.not. ma_alloc_get(mt_dbl, natom, 'charge', h_q,i_q))
c     & call errquit('leps_getgeom: charge not allocated',
c     & natom,MA_ERR)
c      if (.not.ma_alloc_get(mt_dbl, natom, 'mass', h_m,i_m))
c     & call errquit('leps_getgeom: mass not allocated',
c     & natom,MA_ERR)
c      if (.not.ma_alloc_get(mt_byte, 16*natom, 'tags', h_t,i_t))
c     & call errquit('leps_getgeom: tags not allocated',
c     & natom,MA_ERR)
c
cC-----Finished allocating memory for info we get from database
c
cC ----- Fill up the local copies with data from geometry object
cC------- accessed by the geom handle
c      if (.not. geom_cart_get(geom, natom, byte_mb(i_t),
c     & dbl_mb(i_c), dbl_mb(i_q)))
c     & call errquit('leps_getgeom: coord error',0,GEOM_ERR)
c      if (.not. geom_masses_get(geom, natom, dbl_mb(i_m)))
c     & call errquit('leps_getgeom:  mass error',0,GEOM_ERR)
cC--- now compute the leps energy
c      eleps=md_driver1(natom,dbl_mb(i_c),eleps)
cC--- end of leps energy computation
c
cC-----Print what you get out of the database
c      call md_driver0(natom,dbl_mb(i_c),dbl_mb(i_m),byte_mb(i_t))
cC-----Calculate the energy
c      if (.not. rtdb_put(rtdb,'leps:energy',mt_dbl,1,eleps))
c     & call errquit('md_driver:  no energy ', 0,
c     & RTDB_ERR)
c
cC     GET RID OF MEMORY WITH MA_CHOP
cC     DEALLOCATE THE FIRST/LOWEST PART OF THE STACK
cC     THIS ACTION DEALLOCATES ALL ADDITIONAL
cC     MEMORY ADDED TO THE STACK AFTER LOWEST PART
c
c
cC     DESTROY GEOM. OBJECT
c
cC      WRITE(6,*) 'GEOM before destory', geom
c
c      if (.not. geom_destroy(geom))
c     & call errquit('leps_getgeom:  geom_destroy failed',
c     & 911,GEOM_ERR)


      return
      END

