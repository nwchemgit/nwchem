      subroutine grid_initialize(rtdb,nat,c,nsp)
c
c $Id: hnd_elpmap.F 19707 2010-10-29 17:59:36Z d3y133 $
c
      implicit none
#include "prop_grid_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
c
      integer rtdb      ! [Input] rtdb        
      integer nat
      double precision c(3,nat)
      integer nsp(3)
      integer iat
      integer igrid,ix,iy,iz

c     ------------------------
c     SET UP INPUT PARAMETERS
c     ------------------------
c     limits of grid
      if(.not.rtdb_get(rtdb,'prop:pad',MT_DBL,3,pad)) then
         pad = 4.0
      end if
      if(.not.rtdb_get(rtdb,'prop:rmin',MT_DBL,3,rmin)) then
         rmin = MINVAL(c,2)-pad
      end if
      if(.not.rtdb_get(rtdb,'prop:rmax',MT_DBL,3,rmax)) then
         rmax = MAXVAL(c,2)+pad
      end if

c     number of grid points from rtdb or computed
      if(.not.rtdb_get(rtdb,'prop:ngrid',MT_INT,3,ngrid)) then
         step  = 0.2
         ngrid = CEILING((rmax-rmin)/step)+1
      end if
      npoints = PRODUCT(ngrid)
      nsp = ngrid -1
    
c     recalculate step
      step = (rmax-rmin)/nsp

      WRITE(*,*) "------PARAMETERS of gaussian cube file----------"
      WRITE(*,'(A,T12,":",3F6.3)') "step",step
      write(*,'(A,T12,":",3I3)') "ngrid",ngrid
      write(*,'(A,T12,":",I9)') "npoints",npoints
      write(*,'(A,T12,":",3I3)') "nsp",nsp
      write(*,'(A,T12,":",3(F6.3,4X))') "rmax",rmax*cau2ang
      write(*,'(A,T12,":",3(F6.3,4X))') "rmin",rmin*cau2ang
      write(*,*) "------------------------------------------------"
c     ---------------------------
c     MEMORY ALLOCATION FOR GRID
c     --------------------------
      if(.not.ma_alloc_get(MT_DBL, 3*npoints, 'grid array',
     &      h_r, i_r) ) call errquit(
     &      'grid_initialize: unable to allocate heap space',
     &      0, MA_ERR)

c     ---------------------------
c     COMPUTE GRID
c     ---------------------------
      igrid = 0
      do ix = 0, ngrid(1)-1
         do iy = 0, ngrid(2)-1
            do iz = 0, ngrid(3)-1
               igrid = igrid + 1
               dbl_mb(i_r+3*(igrid-1)  ) = rmin(1)+dble(ix)*step(1)
               dbl_mb(i_r+3*(igrid-1)+1) = rmin(2)+dble(iy)*step(2)
               dbl_mb(i_r+3*(igrid-1)+2) = rmin(3)+dble(iz)*step(3)
            end do
         end do
      end do
       
      end

      subroutine grid_write_cube(geom,n,v)
c
c
c
      implicit none
#include "prop_grid_data.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
c
      integer geom
      double precision v(n)
      integer n
      integer i,ii

      integer ix,iy,iz
      integer igrid
      integer un
      double precision charge
      character*16 atag,element
      character*2 symbol
      integer atomicno,nat
      double precision nuc_charge
      double precision c(3)
      character*30 pname

      pname = "grid_write_cube"
      if(n.ne.npoints) 
     + call errquit(pname//'incompatible grid size',911,GEOM_ERR)
      charge = 0.0d0
      un = 54

      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elpmap: geom_ncent',911,GEOM_ERR)
   
      write(un,*) "Marat"
      write(un,*) "Gaussian Cube file"
9498  format(I5,3F12.6)
      write(un,9498)
     &           nat,(rmin(i),i=1,3)
      write(un,9498)
     &           ngrid(1),step(1),0d0,0d0
      write(un,9498)
     &           ngrid(2),0.0d0,step(2),0d0
      write(un,9498)
     &           ngrid(3),0.0d0,0.0d0,step(3)
      do i=1,nat
        if(.not.geom_cent_get(geom,i,atag,c,
     &     nuc_charge)) call
     &     errquit('geom_cent_get',911,GEOM_ERR)
        if (.not. geom_tag_to_element(atag, symbol, 
     &     element, atomicno)) then
           if (symbol .ne. 'bq') call errquit
     &          ('dplotdump: center is neither atom nor bq',
     &          0,0)
        endif
        write(un,9499)
     .        atomicno,charge,c(1),c(2),c(3)
      end do
 9499 format(I5,4F12.6)
 
c      write(un,1000) (v(i),i=1,n)
      do i=1,n,ngrid(3)
        write(un,1000) (v(i+ii),ii=0,ngrid(3)-1)
      end do
1000  format(6E13.5)       
      end subroutine

      subroutine grid_get_r_ptr(nr,i_r_out)
c
c $Id: hnd_elpmap.F 19707 2010-10-29 17:59:36Z d3y133 $
c
      implicit none
#include "prop_grid_data.fh"
#include "errquit.fh"
c
      integer nr,i_r_out
      nr = PRODUCT(ngrid)
      i_r_out = i_r
       
      end

