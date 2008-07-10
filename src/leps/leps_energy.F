      function leps_energy(rtdb)
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
      logical leps_energy
c
      integer geom
      integer natom
      integer h_c,i_c
      integer h_m,i_m
      integer h_q,i_q
      integer h_t,i_t
c
      write(*,*) "doing leps energy"
c     load geometry
      if (.not. geom_create(geom, 'geometry'))
     & call errquit('leps_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     & call errquit('leps_getgeom: no geometry',
     & 0,RTDB_ERR)
      if (.not. geom_ncent(geom,natom))
     & call errquit('leps_getgeom: ncent error',
     & 0,GEOM_ERR)      

      if (.not. geom_print(geom))
     & call errquit('leps_getgeom: ncent error',
     & 0,GEOM_ERR)      
c    
c     allocate memory for coordinates
      if (.not. ma_alloc_get(mt_dbl, 3*natom, 'coordinate' ,
     & h_c,i_c))
     & call errquit('leps_getgeom: coordinate not allocated',
     & 3*natom,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, natom, 'charge', h_q,i_q))
     & call errquit('leps_getgeom: charge not allocated',
     & natom,MA_ERR)
      if (.not.ma_alloc_get(mt_dbl, natom, 'mass', h_m,i_m))
     & call errquit('leps_getgeom: mass not allocated',
     & natom,MA_ERR)
      if (.not.ma_alloc_get(mt_byte, 16*natom, 'tags', h_t,i_t))
     & call errquit('leps_getgeom: tags not allocated',
     & natom,MA_ERR)

      if (.not. geom_cart_get(geom, natom, byte_mb(i_t),
     & dbl_mb(i_c), dbl_mb(i_q)))
     & call errquit('leps_getgeom: coord error',0,GEOM_ERR)
      if (.not. geom_masses_get(geom, natom, dbl_mb(i_m)))
     & call errquit('leps_getgeom:  mass error',0,GEOM_ERR)

      call leps_energy0(natom,dbl_mb(i_c),dbl_mb(i_m),byte_mb(i_t))
 
      if (.not. rtdb_put(rtdb,'leps:energy',mt_dbl,1,0.0))
     & call errquit('leps_energy:  no energy ', 0,
     & RTDB_ERR)
      return
      END

      function leps_energy0(n,c,m,t)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"

      logical leps_energy0
      integer n
      double precision c(3,n)
      double precision m(n)
      character*16 t(n)
c
      integer i

      do i=1,n
       write(*,*) i,c(1,i),c(2,i),c(3,i),m(i),t(i)
      end do    


      return
      END
