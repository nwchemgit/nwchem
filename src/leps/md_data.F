      subroutine md_data_allocate()
      implicit none
#include "md_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"

      integer i
c      natom=n
      if (natom.eq.0)
     &    call errquit('md_data: Natom is not set in md_data.
     &    Allocation failed.',natom,MA_ERR) 
 
      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'coordinate' ,
     & h_c,i_c))
     & call errquit('md_data_getgeom: coordinate not allocated',
     & 3*natom,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, natom, 'charge' ,
     & h_q,i_q))
     & call errquit('md_data_getgeom: charge not allocated',
     & natom,MA_ERR)

      if (.not.ma_alloc_get(mt_dbl, natom, 'mass', h_m,i_m))
     & call errquit('md_data_getgeom: mass not allocated',
     & natom,MA_ERR)

      if (.not.ma_alloc_get(mt_byte, 16*natom, 'tags', h_t,i_t))
     & call errquit('md_data_getgeom: tags not allocated',
     & natom,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'velocities' ,
     & h_v,i_v))
     & call errquit('md_data_getgeom: velocities not allocated',
     & natom,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'accelerations' ,
     & h_a,i_a))
     & call errquit('md_data_getgeom: accelerations not allocated',
     & natom,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'b-term' ,
     & h_b,i_b))
     & call errquit('md_data_getgeom: b-term not allocated',
     & natom,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'd-term' ,
     & h_d,i_d))
     & call errquit('md_data_getgeom: d-term not allocated',
     & natom,MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'force' ,
     & h_f,i_f))
     & call errquit('md_data_getgeom: force not allocated',
     & natom,MA_ERR)

C --- zero out all the alocated memory
      do i=1,3*natom
       dbl_mb(i_c+i-1)=0.0
       dbl_mb(i_v+i-1)=0.0
       dbl_mb(i_a+i-1)=0.0
       dbl_mb(i_b+i-1)=0.0
       dbl_mb(i_d+i-1)=0.0
      enddo
c      write(*,*) 'Done with all the allocation'
c      write(*,*) 'i_c,i_m,i_q,i_v,i_t,i_a,i_b,i_d,i_f'
c      write(*,*) i_c,i_m,i_q,i_v,i_t,i_a,i_b,i_d,i_f

      end 

      subroutine md_data_free_all()
      implicit none
#include "md_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"

      if(.not.ma_free_heap(h_c))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_m))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_q))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_t))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_v))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_a))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_b))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_d))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)
      if(.not.ma_free_heap(h_f))
     &     call errquit('memory deallocation:coord',0,RTDB_ERR)

      end

      function md_data_get_i_c()
      implicit none
#include "md_data.fh"
      integer md_data_get_i_c

      md_data_get_i_c = i_c

      end

      function md_data_get_i_v()
      implicit none
#include "md_data.fh"
      integer md_data_get_i_v

      md_data_get_i_v = i_v

      end

      function md_data_get_i_a()
      implicit none
#include "md_data.fh"
      integer md_data_get_i_a

      md_data_get_i_a = i_a

      end

      function md_data_get_i_b()
      implicit none
#include "md_data.fh"
      integer md_data_get_i_b

      md_data_get_i_b = i_b

      end

      function md_data_get_i_d()
      implicit none
#include "md_data.fh"
      integer md_data_get_i_d

      md_data_get_i_d = i_d

      end

      function md_data_get_i_f()
      implicit none
#include "md_data.fh"
      integer md_data_get_i_f

      md_data_get_i_f = i_f

      end

      function md_data_get_i_m()
      implicit none
#include "md_data.fh"
      integer md_data_get_i_m

      md_data_get_i_m = i_m

      end

      function md_data_get_i_q()
      implicit none
#include "md_data.fh"
      integer md_data_get_i_q
      md_data_get_i_q = i_q
      end

      function md_data_get_i_t()
      implicit none
#include "md_data.fh"
      integer md_data_get_i_t
      md_data_get_i_t = i_t
      end

      subroutine md_data_set_natom(n)
#include "md_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer n 

c      if (.not. geom_ncent(geom,n))
c     & call errquit('md_data: ncent error',
c     & 0,GEOM_ERR)

      natom=n

      end
 
      function md_data_get_natom()
      implicit none
#include "md_data.fh"
       integer md_data_get_natom
       md_data_get_natom = natom

      end

      subroutine md_data_geom_copy(rtdb,geom)
#include "md_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer rtdb,geom


      if (.not. geom_masses_get(geom,natom,dbl_mb(i_m) ))
     & call errquit('md_data_getgeom:  mass error',0,GEOM_ERR)

      if (.not. geom_cart_get(geom, n, byte_mb(i_t),
     &                         dbl_mb(i_c),dbl_mb(i_q) ))
     & call errquit('md_data_getgeom: coord error',0,GEOM_ERR)
c       write(*,*) 'Geom_copy:Filled cartesians and masses  '


      end

      subroutine md_data_geom_set(rtdb)
#include "md_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer rtdb,geom

      if (.not. geom_create(geom, 'geom1'))
     & call errquit('md_driver_getgeom: failed creating geometry',
     & 0,GEOM_ERR)

c      write(*,*) 'Open a geometry' 
C     load geometry into the geom object from database
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     & call errquit('md_data: no geometry',
     & 0,RTDB_ERR)

      if (.not. geom_ncent(geom,natom))
     & call errquit('md_data: ncent error',
     & 0,GEOM_ERR)

      if (.not. geom_cart_set(geom, natom, byte_mb(i_t),
     &   dbl_mb(i_c), dbl_mb(i_q)))
     &   call errquit('md_data_putgeom: coord error',0,GEOM_ERR)

      if(.not.geom_rtdb_store(rtdb,geom,"geometry"))
     & call errquit('md_driver:Failed to store geom to rtdb',0,RTDB_ERR)

      if (.not. geom_destroy(geom))
     & call errquit('md_driver_getgeom: failed geom_destroy',
     & 0,GEOM_ERR)

      end
   
      subroutine md_data_inp_set(rtdb)
#include "md_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer rtdb

      if (.not.rtdb_get(rtdb,"mymd:temp",mt_dbl,1,temp))
     >      call errquit('failed to get temp',0,
     >       RTDB_ERR)
      if (.not.rtdb_get(rtdb,"mymd:timestep",mt_dbl,1,dt))
     >      call errquit('failed to get timestep',0,
     >       RTDB_ERR)
      if (.not.rtdb_get(rtdb,"mymd:nsteps",mt_int,1,nsteps))
     >      call errquit('failed to get nsteps',0,
     >       RTDB_ERR)

      end
      
      function md_data_get_timestep()
      implicit none
#include "md_data.fh"
      double precision md_data_get_timestep

      md_data_get_timestep = dt

      end
 
      function md_data_get_nsteps()
      implicit none
#include "md_data.fh"
      integer md_data_get_nsteps

      md_data_get_nsteps = nsteps

      end

      function md_data_get_temp()
      implicit none
#include "md_data.fh"
      double precision md_data_get_temp

      md_data_get_temp = temp

      end

      function md_data_set_temp(t)
      implicit none
#include "md_data.fh"
      logical md_data_set_temp
      double precision t   
      temp=t
      md_data_set_temp=.true.
      end

c $Id$
