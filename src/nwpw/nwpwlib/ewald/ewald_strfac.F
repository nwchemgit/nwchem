*
* $Id: ewald_strfac.F,v 1.7 2008-06-10 19:55:11 d3p708 Exp $
*
      subroutine ewald_strfact_init()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "ewald_strfac.fh"
cccccccc arguments
cccccccc locals
      integer taskid,np
      integer nion,nx,ny,nz,nt,nloc,xtra
      logical value
cccccccc externals
      integer ion_nion
      integer ewald_grid_nx
      integer ewald_grid_ny
      integer ewald_grid_nz
      external ion_nion
      external ewald_grid_nx
      external ewald_grid_ny
      external ewald_grid_nz
cccccccc exec
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nion=ion_nion()
      nx=ewald_grid_nx()
      ny=ewald_grid_ny()
      nz=ewald_grid_nz()
      nt=nx*ny*nz
      nloc=nt/np
      xtra=nt-np*nloc
      if (taskid.lt.xtra) nloc=nloc+1   
      value=MA_alloc_get(mt_dcpl,(nx*nion),'ewx1',ewx1(2),ewx2(1))
      value=value.and.MA_alloc_get(mt_dcpl,(ny*nion),'ewx2',ewx2(2),
     > ewx2(1))
      value=value.and.MA_alloc_get(mt_dcpl,(nz*nion),'ewx3',ewx3(2),
     > ewx3(1))
      if (.not.value) 
     >  call errquit("EWALD_STRUCTFACTOR_INIT: OUT OF HEAP MEMORY",0,
     >               MA_ERR)
      return 
      end
cccccccccccccccccccccccc
      subroutine ewald_strfact_end()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "ewald_strfac.fh"
      logical value
cccccccc exec
      value=.true. 
      value=value.and.MA_free_heap(ewx1(2))
      value=value.and.MA_free_heap(ewx2(2))
      value=value.and.MA_free_heap(ewx3(2))
      if (.not.value) then
        call errquit("ewald_strfact_end: error free heap",0,MA_ERR)
      end if
      return
      end
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
      subroutine ewald_phafac()
      implicit none
#include "mafdecls.fh"
#include "ewald_strfac.fh"
cccccccc locals
      integer i,k
      integer nx,ny,nz,nxh,nyh,nzh
      integer nion   
      complex*16 cw1,cw2,cw3
      real*8 sw1,sw2,sw3,pi
ccccccccc externals
      integer ion_nion,ewald_grid_nx
      integer ewald_grid_ny,ewald_grid_nz
      real*8 lattice_unitg,ion_rion
      external ion_nion
      external ewald_grid_nz
      external ewald_grid_ny,ewald_grid_nx
      external lattice_unitg,ion_rion
ccccccc start of code       
      call nwpw_timing_start(8)
      pi=4.0d0*datan(1.0d0)
      nion=ion_nion()
      nx=ewald_grid_nx()
      ny=ewald_grid_ny()
      nz=ewald_grid_nz()
      nxh=nx/2
      nyh=ny/2
      nzh=nz/2
      do i=1,nion
        sw1=lattice_unitg(1,1)*ion_rion(1,i)+
     >       lattice_unitg(2,1)*ion_rion(2,i)+
     >       lattice_unitg(3,1)*ion_rion(3,i)+pi
        sw2=lattice_unitg(1,2)*ion_rion(1,i)+
     >       lattice_unitg(2,2)*ion_rion(2,i)+
     >       lattice_unitg(3,2)*ion_rion(3,i)+pi
        sw3=lattice_unitg(1,3)*ion_rion(1,i)+
     >       lattice_unitg(2,3)*ion_rion(2,i)+
     >       lattice_unitg(3,3)*ion_rion(3,i)+pi
        cw1=dcmplx(dcos(sw1),-dsin(sw1))
        cw2=dcmplx(dcos(sw2),-dsin(sw2))
        cw3=dcmplx(dcos(sw3),-dsin(sw3))
        dcpl_mb(ewx1(1)+(i-1)*nx)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ewx2(1)+(i-1)*ny)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ewx3(1)+(i-1)*nz)=dcmplx(1.0d0,0.0d0)

        do k=1,nxh
          dcpl_mb(ewx1(1)+k+(i-1)*nx)
     >         = dcpl_mb(ewx1(1)+k-1+(i-1)*nx)*cw1
          dcpl_mb(ewx1(1)+nx-k+(i-1)*nx)
     >         =dconjg(dcpl_mb(ewx1(1)+k+(i-1)*nx))
        end do

        do k=1,nyh
          dcpl_mb(ewx2(1)+k+(i-1)*ny)
     >         = dcpl_mb(ewx2(1)+k-1+(i-1)*ny)*cw2
          dcpl_mb(ewx2(1)+ny-k+(i-1)*ny)
     >         =dconjg(dcpl_mb(ewx2(1)+k+(i-1)*ny))
        end do

        do k=1,nzh
          dcpl_mb(ewx3(1)+k+(i-1)*nz)
     >         = dcpl_mb(ewx3(1)+k-1+(i-1)*nz)*cw3
          dcpl_mb(ewx3(1)+nz-k+(i-1)*nz)
     >         =dconjg(dcpl_mb(ewx3(1)+k+(i-1)*nz))
        end do
        dcpl_mb(ewx1(1)+nxh+(i-1)*nx)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ewx2(1)+nyh+(i-1)*ny)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ewx3(1)+nzh+(i-1)*nz)=dcmplx(0.0d0, 0.0d0)  
      end do 
      call nwpw_timing_end(8)
      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ewald_strfac(ii,exi)
      implicit none
      integer ii
      complex*16 exi(*)
#include "mafdecls.fh"
#include "ewald_strfac.fh"
ccccccc locals
      integer i,j,k,indx
      integer nx,ny,nz
      integer im1
      integer taskid,np,dutask,endx
      integer ewald_grid_nx,ewald_grid_ny,ewald_grid_nz
      logical ewald_do_point
      external ewald_grid_nx,ewald_grid_ny,ewald_grid_nz
      external ewald_do_point

      call nwpw_timing_start(8)
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nx=ewald_grid_nx()
      ny=ewald_grid_ny()
      nz=ewald_grid_nz()
      dutask=0
      im1=ii-1
      indx=1
      endx=0
      do k=0,nz-1
      do j=0,ny-1
      do i=0,nx-1
         if ((dutask.eq.taskid)) then
          if (ewald_do_point(endx)) then
           exi(indx)=dcpl_mb(ewx1(1)+i+im1*nx)*
     >           dcpl_mb(ewx2(1)+j+im1*ny)*
     >           dcpl_mb(ewx3(1)+k+im1*nz)
           indx=indx+1
          end if
          endx=endx+1
         end if
         dutask=mod(dutask+1,np)
      end do
      end do
      end do
      call nwpw_timing_end(8)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccc
      

    


  
      
