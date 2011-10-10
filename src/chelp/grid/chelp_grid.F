      subroutine chelp_grid_gen(rtdb)
      implicit none
#include "errquit.fh"
#include "chelp_grid_data.fh"
      integer rtdb
      call chelp_grid_gen_from geom(rtdb,"geometry")
      return
      end
      subroutine chelp_grid_gen_from geom(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "chelp_grid_data.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
      character*(*) geom
c     
      integer igeom
      integer nt
      integer h_c,i_c
      integer h_t,i_t
      integer h_q,i_q
      integer h_r,i_r
      character*30 pname
      double precision rpad,rad
      double precision rmax(3),rmin(3),delr
      integer nr(3)
      integer nrt
      logical oprune,oext
c
      pname = "chelp_grid_gen_from geom"

c
c     get geometry information
c     ------------------------
      if(.not.geom_create(igeom,"chelp geometry"))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)
c      
      if(.not.geom_rtdb_load(rtdb,igeom,geom))
     + call errquit(pname//'Failed to load geometry',0, GEOM_ERR)
c      
      if(.not.geom_ncent(igeom, nt))
     $           call errquit(pname,0,RTDB_ERR)
c    
      if(.not.ma_push_get(mt_dbl,3*nt,'tmp h_c',h_c,i_c))
     + call errquit(pname//'Failed to allocate memory for h_ftmp',0,0)
c      
      if(.not.ma_push_get(mt_dbl,nt,'tmp h_q',h_q,i_q))
     + call errquit('qmmm: Failed to allocate memory for q',nt,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*nt,'tmp h_t',h_t,i_t))
     + call errquit('qmmm: Failed to allocate memory for t',nt,
     &       MA_ERR)

      if(.not.geom_cart_get(igeom,nt,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q)))
     + call errquit(pname//'Failed to get geometry',0, GEOM_ERR)
c
c     compute bounds
c     --------------
      if (.not.rtdb_get(rtdb,"chelp:rpad",
     >                   mt_dbl,1,rpad)) 
     >   rpad = 1.0d0
       call chelp_grid_bounds(nt,dbl_mb(i_c),rpad,rmin,rmax)
c
c     compute grid size
c     -----------------
      if (.not.rtdb_get(rtdb,"chelp:delr",
     >                   mt_dbl,1,delr)) 
     >   delr =0.1d0
       call chelp_grid_size_compute(rmin,rmax,delr,nr)
       nrt = nr(1)*nr(2)*nr(3)
c
c     generate uniform grid
c     ---------------------
      if(.not.ma_push_get(mt_dbl,3*nrt,'tmp h_r',h_r,i_r))
     + call errquit('qmmm: Failed to allocate memory for r',nr,
     &       MA_ERR)
      call chelp_grid_gen_uniform(rmin,rmax,nr,dbl_mb(i_r))
c
c     prune the grid if necessary
c     ---------------------------
      if (.not.rtdb_get(rtdb,"chelp:prune",
     >                   mt_log,1,oprune)) 
     >   oprune = .false.
      if(oprune) then
         if (.not.rtdb_get(rtdb,"chelp:prune:exterior",
     >                     mt_log,1,oext)) 
     >   oext = .false.
         if (.not.rtdb_get(rtdb,"chelp:rad",
     >                      mt_dbl,1,rad)) 
     >      rad=1.3d0
c        this modifies nrt,and i_r
         call chelp_grid_prune(oext,
     >      nrt,dbl_mb(i_r),nt,dbl_mb(i_c),rad)
      end if
c
c     save grid
c     --------
      call chelp_grid_data_alloc_put(3*nrt,dbl_mb(i_r))
c
c     clean up
c     -------
      if(.not.ma_chop_stack(h_c))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)

      call chelp_grid_data_print()
      return
      end
      subroutine chelp_grid_bounds(nt,c,rpad,rmin,rmax)
c     calculates bounding box
c     nt           [in] number of atoms
c     c(3,nt)      [in] coordinates atoms
c     rpad         [in] padding distance
c     rmin(3)      [in] minimum x,y,z
c     rmax(3)      [in] maximum x,y,z
      implicit none
#include "stdio.fh"
      double precision rpad
      integer nt
      double precision c(3,nt)
      double precision rmin(3), rmax(3)
c 
      integer i,k
c
c     find the extrema of the      
c     molecular geometry.                                                
c                                                                     
      do k=1,3
        rmin(k) = c(k,1)
        rmax(k) = c(k,1)
      end do
      do  i=1,nt                                                
      do k=1,3
        write(*,*) "c(k,i)=",c(k,i)
        if (c(k,i) .gt. rmax(k)) rmax(k) = c(k,i) 
        if (c(k,i) .lt. rmin(k)) rmin(k) = c(k,i)   
      end do
      end do
c      
      do k=1,3
        rmin(k) = rmin(k) - rpad
        rmax(k) = rmax(k) + rpad
      end do
c                                                                     
      write(luout,*) ' xmax = ',rmax(1),'  xmin = ',rmin(1)   
      write(luout,*) ' ymax = ',rmax(2),'  ymin = ',rmin(2)  
      write(luout,*) ' zmax = ',rmax(3),'  zmin = ',rmin(3) 
c                                                                     
      return                                                          
      end                                                             
      subroutine chelp_grid_size_compute(rmin,rmax,delr,nr)
      implicit none
#include "stdio.fh"
c     rmin(3)   [in]   box bounds
c     rmax(3)   [in]   box bounds
c     delr      [in]   grid spacing
c     nr        [out]  grid size
      double precision rmax(3),rmin(3)
      double precision delr
      integer nr(3)
c
      integer k
      double precision delr1,range
c
c     the number of grid points is at least 2 (edges)
c     ----------------------------------------------- 
      do k=1,3
        range = (rmax(k)-rmin(k))
        nr(k) = 0
c       count interior points
10      nr(k) = nr(k) +1
        delr1=range/(nr(k)+1)
        if(delr1.gt.delr) goto 10
c       add edge points
        nr(k) = nr(k)+2
      end do
c                                                                     
      write(luout,*) 'grid size', (nr(k),k=1,3)
c                                                                     
      return                                                          
      end                                                             
      subroutine chelp_grid_gen_uniform(rmin,rmax,nr,r)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
c     rmin(3)   [in] box bounds
c     rmax(3)   [in] box bounds
c     nr(3)     [in] grid size
c     r(3,nr)   [out] grid 
      double precision rmax(3),rmin(3)
      integer nr(3)
      double precision r(3,*)
c
      integer i1,i2,i3
      double precision p1,p2,p3
      integer i
      character*30 pname
      integer nrt
      double precision dx,dy,dz

      pname = "chelp_grid_gen_uniform"
      nrt=nr(1)*nr(2)*nr(3)
      dx = (rmax(1)-rmin(1))/(nr(1)-1)
      dy = (rmax(2)-rmin(2))/(nr(2)-1)
      dz = (rmax(3)-rmin(3))/(nr(3)-1)
c                                    
      i=0
      do i1=1,nr(1)
        p1=rmin(1)+dx*(i1-1)
        do i2=1,nr(2)
        p2=rmin(2)+dy*(i2-1)
          do i3=1,nr(3)
            p3=rmin(3)+dz*(i3-1)
            i = i+1
            if(i.gt.nrt)
     +        call errquit(pname//"out of bounds",i,0)
            r(1,i) = p1
            r(2,i) = p2
            r(3,i) = p3
          end do
        end do
      end do

c                                                                     
      return                                                          
      end                                                             
      subroutine chelp_grid_prune(oext,nrt,r,nt,c,rad)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
c     oext      [in] if true exterior opints are pruned
c     nrt       [inout] grid size
c     r(3,nrt)  [inout] grid 
c     nt        [in] number of atoms
c     c(3,nt)   [in] atomic coordinates
      logical oext
      integer nrt
      double precision r(3,nrt)
      integer nt
      double precision c(3,nt)
      double precision rad
c
      integer ipoint,k,ir
      double precision p1,p2,p3
      integer i
      character*30 pname
      integer h_r,i_r
      double precision dist
      logical oskip
c      
      pname = "chelp_grid_gen_uniform"
c
      if(.not.ma_push_get(mt_dbl,3*nrt,'tmp h_r',h_r,i_r))
     + call errquit(pname//'Failed to allocate memory for h_r',0,0)
c
      ipoint = 0
      do ir=1,nrt
        do i=1,nt
          dist = (r(1,ir) - c(1,i))**2 + 
     +           (r(2,ir) - c(2,i))**2 +
     +           (r(3,ir) - c(3,i))**2 
          dist = dsqrt(dist)                                           
          oskip = dist.lt.rad
          if(oext) oskip=.not.oskip
          if(oskip) goto 20
        end do
        do k=1,3
        ipoint = ipoint+1
        dbl_mb(i_r+ipoint-1) = r(k,ir)
        end do
20      continue
      end do
 
      call dcopy(ipoint,dbl_mb(i_r),1,r,1)
      nrt = ipoint/3
      if(.not.ma_chop_stack(h_r))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)
c                                                                     
      return                                                          
      end                                                             

c $Id$
