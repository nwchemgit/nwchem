      subroutine prop_grid_initialize(rtdb,nat,c)
c
c $Id$
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
      integer iat
      integer igrid,ix,iy,iz
      integer nsp(3)
      character*256 buffer

c     ------------------------
c     SET UP INPUT PARAMETERS
c     ------------------------
c     limits of grid
      if(.not.rtdb_get(rtdb,'prop:grid:pad',MT_DBL,3,pad)) then
         pad = 4.0 * cang2au
      end if
      if(.not.rtdb_get(rtdb,'prop:grid:rmin',MT_DBL,3,rmin)) then
         rmin = MINVAL(c,2)-pad
      end if
      if(.not.rtdb_get(rtdb,'prop:grid:rmax',MT_DBL,3,rmax)) then
         rmax = MAXVAL(c,2)+pad
      end if

c     number of grid points from rtdb or computed
      if(.not.rtdb_get(rtdb,'prop:grid:ngrid',MT_INT,3,ngrid)) then
        if(.not.rtdb_get(rtdb,'prop:grid:step',MT_DBL,1,step)) then
         step  = 0.2/cau2ang
         ngrid(1) = min(CEILING((rmax(1)-rmin(1))/step(1))+1,50)  
         ngrid(2) = min(CEILING((rmax(2)-rmin(2))/step(1))+1,50)
         ngrid(3) = min(CEILING((rmax(3)-rmin(3))/step(1))+1,50)
        else
         step  = step/cau2ang
         ngrid(1) = CEILING((rmax(1)-rmin(1))/step(1))+1
         ngrid(2) = CEILING((rmax(2)-rmin(2))/step(1))+1
         ngrid(3) = CEILING((rmax(3)-rmin(3))/step(1))+1
        endif
      end if
      npoints = PRODUCT(ngrid)
      nsp = ngrid -1
    
c     recalculate step
      step = (rmax-rmin)/nsp
      if(ga_nodeid().eq.0) then
      buffer = "PARAMETERS of gaussian cube file"
      call util_print_centered(luout, trim(buffer),.true., .true.)
      WRITE(*,'(1X,A,T12,":",3F12.3)') "step (Ang)",step*cau2ang
      write(*,'(1X,A,T12,":",3I4)') "ngrid",ngrid
      write(*,'(1X,A,T12,":",I9)') "npoints",npoints
      write(*,'(1X,A,T12,":",3(F12.3,4X))') "rmax (Ang)",rmax*cau2ang
      write(*,'(1X,A,T12,":",3(F12.3,4X))') "rmin (Ang)",rmin*cau2ang
      write(*,*) "------------------------------------------------"
      end if
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

      subroutine prop_grid_destroy()
c
c $Id$
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

      if(.not.ma_free_heap(h_r) ) call errquit(
     &      'unable to deallocate heap space',
     &      0, MA_ERR)
      end

      subroutine prop_grid_writecubehead(geom,un,fname)
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
      character*(*) fname
      integer i,ii
      integer un ! [out]

      integer ix,iy,iz
      integer igrid
      double precision charge
      character*16 atag,element
      character*2 symbol
      integer atomicno,nat
      double precision nuc_charge
      double precision c(3)
      character*30 pname
      logical util_get_io_unit
      external util_get_io_unit

      pname = "grid_write_cube"
      if(.not.util_get_io_unit(un))
     >    call errquit(pname//"no free file unit",0,0)

      open(unit=un,file=fname,
     >       form='formatted',status='unknown')
      charge = 0.0d0
c      un = 54

      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elpmap: geom_ncent',911,GEOM_ERR)
   
      write(un,*) "NWChem"
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
 
      end subroutine
      subroutine prop_grid_writecubegrid(n,v,un)
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
#include "rtdb.fh"
c
      double precision v(*)
      integer n
      integer i,ii
      integer un ! [in]

      integer ix,iy,iz
      integer igrid
      character*30 pname
      logical util_get_io_unit
      external util_get_io_unit
c
      integer remaindr

      pname = "grid_writecubegrid"
      remaindr=n-(n/ngrid(3))*ngrid(3)
      if(remaindr.ne.0)  then
         write(6,*) ' n ngrid(3) ',n,ngrid(3)
      call errquit(pname//'incompatible grid size',
     N     remaindr,GEOM_ERR)
      endif
c
      do i=1,n,ngrid(3)
        write(un,1000) (v(i+ii),ii=0,ngrid(3)-1)
      end do
1000  format(6E13.5)       
      end subroutine
      subroutine prop_grid_write_cube(geom,n,v,fname)
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
      integer n
      double precision v(n)
      character*(*) fname
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
      logical util_get_io_unit
      external util_get_io_unit

      pname = "grid_write_cube"
      if(.not.util_get_io_unit(un))
     >    call errquit(pname//"no free file unit",0,0)

      open(unit=un,file=fname,
     >       form='formatted',status='unknown')
      if(n.ne.npoints) 

     + call errquit(pname//'incompatible grid size',911,GEOM_ERR)
      charge = 0.0d0
c      un = 54

      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elpmap: geom_ncent',911,GEOM_ERR)
   
      write(un,*) "NWChem"
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
      close(un)
      end subroutine

      subroutine prop_grid_get_r_ptr(nr,i_r_out,ngrid_out)
c
c $Id$
c
      implicit none
#include "prop_grid_data.fh"
#include "errquit.fh"
c
      integer nr,i_r_out
      integer ngrid_out(3)
      nr = PRODUCT(ngrid)
      i_r_out = i_r
      ngrid_out(1)=ngrid(1)
      ngrid_out(2)=ngrid(2)
      ngrid_out(3)=ngrid(3)
       
      end
c
c     Writes an unformatted CUBE file.
c

      subroutine prop_grid_write_cube_unformatted(geom,n,v,
     >i_prp_c, fname)
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
      integer n
      double precision v(n)
      character*(*) fname
      integer i,ii

      integer ix,iy,iz
      integer igrid
      integer un
      double precision charge
      character*16 atag,element
      character*2 symbol
      integer atomicno,nat
      integer i_prp_c
      double precision nuc_charge
      double precision c(3)
      double precision xp,yp,zp
      character*30 pname
      logical util_get_io_unit
      external util_get_io_unit

      pname = "grid_write_cube_unformatted"
      if(.not.util_get_io_unit(un))
     >    call errquit(pname//"no free file unit",0,0)

      open(unit=un,file=fname,
     >       form='formatted',status='unknown')
      charge = 0.0d0

      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elpmap: geom_ncent',911,GEOM_ERR)

      write(un,*) "NWChem"
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

      do i=1,n
         xp = dbl_mb(i_prp_c  +3*(i-1))
         yp = dbl_mb(i_prp_c+1+3*(i-1))
         zp = dbl_mb(i_prp_c+2+3*(i-1))
        write(un,1000) xp,yp,zp,v(i)
      end do
1000  format(4F13.6)
      close(un)
      end subroutine

