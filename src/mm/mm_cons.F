      subroutine mm_cons_reaction(irtdb,do_grad,etot,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mm_geom_data.fh"
      integer irtdb
      logical do_grad
      double precision g(3,*)
      integer n
      double precision etot

      character*32 pname
      integer igeom 
      integer ncent
      integer i_c,h_c
c     integer nfg
c     integer i_g,i_ifg
      integer ns
      integer i_is,h_is
      integer i_as,h_as
      integer i_rs,h_rs
      double precision k0,r0
      double precision econs,deds,s
      integer master

      master = 0
c
      if (.not. rtdb_get(irtdb,'qmmm:n_spring',mt_int,1,ns))
     $     ns = 0

      if(ns.eq.0) return

     
      pname = "mm_cons_reaction"

c     if(ga_nodeid().eq.master) 
c    + write(*,*) "in",pname

c     nfg = qmmm_geom_get_nfg()

      if(.not.ma_alloc_get(mt_int,2*ns,'qmmm rcons ind',h_is,i_is))
     + call errquit(pname//'Failed to allocate memory for is',
     + 2*ns, MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,ns,'qmmm rcons r',h_rs,i_rs))
     + call errquit(pname//'Failed to allocate memory for rs',
     + ns, MA_ERR)
        call dfill(ns,0.0d0,dbl_mb(i_rs),1)

      if(.not.ma_alloc_get(mt_dbl,ns,'qmmm rcons coeff',h_as,i_as))
     + call errquit(pname//'Failed to allocate memory for as',
     + ns, MA_ERR)


      if (.not. rtdb_get(irtdb,'qmmm:i_spring',mt_int,
     +                    2*ns,int_mb(i_is)))
     + call errquit(pname//'Failed to allocate memory for is',
     + 3*ns, MA_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:a_spring',mt_dbl,
     +                    ns,dbl_mb(i_as)))
     +  call dfill(ns,1.0d0,dbl_mb(i_as),1)

      if (.not. rtdb_get(irtdb,'qmmm:k_spring',mt_dbl,
     +                    1,k0))
     + call errquit(pname//'Failed to get k0',
     + 1, RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:r_spring',mt_dbl,
     +                    1,r0))
     + call errquit(pname//'Failed to get r0',
     + 1, RTDB_ERR)


c    --------------------------------------
c    get coordinates out of geometry
c    --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('mm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('mm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_ncent(igeom, ncent) ) 
     >    call errquit("mm:geom_ncent",0,0)

      if(nact.ne.ncent) 
     + call errquit(pname//': geometry mismatch',0, GEOM_ERR)

      if(.not.ma_alloc_get(mt_dbl,3*ncent,'ctmp',h_c,i_c))
     + call errquit( pname//'Failed to allocate memory for c',
     + 3*ncent, MA_ERR)

      if(.not.geom_cart_coords_get(igeom,dbl_mb(i_c)))
     + call errquit(pname//': Failed to get geometry',0, GEOM_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit(pname//'Failed to destroy geometry',0, GEOM_ERR)

      call qmmm_cons_spring_egrad(ns,nact,
     >                       do_grad,
     >                       k0,
     >                       r0,
     >                       econs,
     >                       int_mb(i_is),
     >                       dbl_mb(i_as),
     >                       dbl_mb(i_rs),
     >                       dbl_mb(i_c),
     >                       deds,
     >                       s,
     >                       g)
      etot  = etot + econs
c
      if (.not. rtdb_put(irtdb,'qmmm:spring_e',mt_dbl,
     +                    1,econs))
     + call errquit(pname//'Failed to set spring_e',
     + 1, RTDB_ERR)

      if (.not. rtdb_put(irtdb,'qmmm:spring_s',mt_dbl,
     +                    1,s))
     + call errquit(pname//'Failed to set spring_s',
     + 1, RTDB_ERR)


      if (.not. rtdb_put(irtdb,'qmmm:spring_deds',mt_dbl,
     +                    1,deds))
     + call errquit(pname//'Failed to set spring_deds',
     + 1, RTDB_ERR)

      if(ga_nodeid().eq.master) then
        write(*,*) "spring e,deds,s",econs,deds,s
      end if
 
      if(.not.ma_free_heap(h_c))
     + call errquit(pname//'Failed to deallocate memory for is',
     + 3*ns, MA_ERR)

      if(.not.ma_free_heap(h_as))
     + call errquit(pname//'Failed to deallocate memory for is',
     + 3*ns, MA_ERR)

      if(.not.ma_free_heap(h_rs))
     + call errquit(pname//'Failed to deallocate memory for rs',
     + 3*ns, MA_ERR)

      if(.not.ma_free_heap(h_is))
     + call errquit(pname//'Failed to deallocate memory for is',
     + 3*ns, MA_ERR)

      return
      end
