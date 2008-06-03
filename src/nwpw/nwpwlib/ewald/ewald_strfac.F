*
* $Id: ewald_strfac.F,v 1.1 2008-06-03 17:25:58 d3p708 Exp $
*
      subroutine ewald_strfact_init(nx,ny,nz,taskid,np)
      implicit none
#include "ewald.fh"
cccccccc arguments
      integer nx,ny,nz,taskid,np
cccccccc locals
      logical value
      integer i,j,k,nxh,nion
cccccccc exec
      nion=ion_nion()
      value=MA_alloc_get(mt_dcpl,(nx*nion),'ewx1',ewx1(2),ewx2(1))
      value=value.and.MA_alloc_get(mt_dcpl,(ny*nion),'ewy1',ewy1(2),
     > ewy2(1))
      value=value.and.MA_alloc_get(mt_dcpl,(nx*nion),'ewz1',ewz1(2),
     > ewz2(1))
      value=value.and.MA_alloc_get(mt_dcpl,(nbrillq*nion),'cxr',cxr(1),
     > cxr(2))
      value=value.and.MA_alloc_get(mt_int,(2*nbrillq+2),'i_indx_list',
     >    i_indx_list(2),i_indx_list(1))
      value=value.and.MA_alloc_get(mt_int,(2*nbrillq+2),'j_indx_list',
     >    j_indx_list(2),j_indx_list(1))
      value=value.and.MA_alloc_get(mt_int,(2*nbrillq+2),'k_indx_list',
     >    k_indx_list(2),k_indx_list(1))

      if (.not.value) 
     >  call errquit("EWALD_STRUCTFACTOR_INIT: OUT OF HEAP MEMORY",0,
     >               MA_ERR)
      do nb=0,nbrillq
        value=MA_alloc_get(mt_int,enloc,'i_indx0',
     >          i_indx(2,nb),i_indx(1,nb)
        value=value.and.
     >          MA_alloc_get(mt_int,enloc,'j_indx0',
     >          j_indx(2,nb),j_indx(1,nb)
        value=value.and.
     >          MA_alloc_get(mt_int,enloc,'k_indx0',
     >          k_indx(2,nb),k_indx(1,nb)
        if (.not.value)
     >       call errquit("ewald_stuctfactor:out of heap memory",nb+1,
     >                MA_ERR)       
        int_mb(i_indx_list(1)+2*nb)=i_indx(1)
        int_mb(i_indx_list(1)+2*nb+1)=i_indx(2)
        int_mb(j_indx_list(1)+2*nb)=j_indx(1)
        int_mb(j_indx_list(1)+2*nb+1)=j_indx(2)
        int_mb(k_indx_list(1)+2*nb)=k_indx(1)
        int_mb(k_indx_list(1)+2*nb+1)=k_indx(2)

        dutask=0
        cindx=0
        do k=1,nz
        do j=1,ny
        do i=1,nx
          if (dutask.eq.taskid) then
             int_mb(i_indx(1)+cindx-1)=i-1
             int_mb(j_indx(1)+cindx-1)=j-1
             int_mb(k_indx(1)+cindx-1)=k-1
             cindx=cindx+1
          end if
          dutask=mod(dutask+1,np)
        end do
        end do
        end do
      end do
      return 
      end
cccccccccccccccccccccccc
      subroutine ewald_strfact_end()
      implicit none
      integer nb
      logical value
      value=.true. 
      do nb=nbrillq,0,-1
         i_indx=int_mb(i_indx_list(1)+2*nb+1)
         j_indx=int_mb(j_indx_list(1)+2*nb+1)
         k_indx=int_mb(k_indx_list(1)+2*nb+1)
         value=value.and.MA_free_heap(k_indx)
         value=value.and.MA_free_heap(j_indx)
         value=value.and.MA_free_heap(i_indx)
      end do
      value=value.and.MA_free_heap(i_indx_list(2))
      value=value.and.MA_free_heap(j_indx_list(2))
      value=value.and.MA_free_heap(k_indx_list(2))
      value=value.and.MA_free_heap(ex1(2))
      value=value.and.MA_free_heap(ex2(2))
      value=value.and.MA_free_heap(ex3(2))
      value=value.and.MA_free_heap(cxr(2))
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
      external lattice_unitg,ion_rion       
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
        dcpl_mb(ex1(1)+(i-1)*nx)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex2(1)+(i-1)*ny)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex2(1)+(i-1)*nz)=dcmplx(1.0d0,0.0d0)

        do k=1,nxh
          dcpl_mb(ex1(1)+k+(i-1)*nx)
     >         = dcpl_mb(ex1(1)+k-1+(i-1)*nx)*cw1
          dcpl_mb(ex1(1)+nx-k+(i-1)*nx)
     >         =dconjg(dcpl_mb(ex1(1)+k+(i-1)*nx))
        end do

        do k=1,nyh
          dcpl_mb(ex2(1)+k+(i-1)*ny)
     >         = dcpl_mb(ex2(1)+k-1+(i-1)*ny)*cw2
          dcpl_mb(ex2(1)+ny-k+(i-1)*ny)
     >         =dconjg(dcpl_mb(ex2(1)+k+(i-1)*ny))
        end do

        do k=1,nzh
          dcpl_mb(ex3(1)+k+(i-1)*nz)
     >         = dcpl_mb(ex3(1)+k-1+(i-1)*nz)*cw3
          dcpl_mb(ex3(1)+nz-k+(i-1)*nz)
     >         =dconjg(dcpl_mb(ex3(1)+k+(i-1)*nz))
        end do
        dcpl_mb(ex1(1)+nxh+(i-1)*nx)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex2(1)+nyh+(i-1)*ny)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex3(1)+nzh+(i-1)*nz)=dcmplx(0.0d0, 0.0d0)  
      end do 
      call nwpw_time_end(8)
      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
      subroutine ewald_phafac_rion(rion)
      implicit none
      real*8 rion(3,*)
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
      real*8 lattice_unitg
      external ion_nion
      external lattice_unitg
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
        sw1=lattice_unitg(1,1)*rion(1,i)+
    >       lattice_unitg(2,1)*rion(2,i)+
    >       lattice_unitg(3,1)*rion(3,i)+pi
        sw2=lattice_unitg(1,2)*rion(1,i)+
    >       lattice_unitg(2,2)*rion(2,i)+
    >       lattice_unitg(3,2)*rion(3,i)+pi
        sw3=lattice_unitg(1,3)*rion(1,i)+
    >       lattice_unitg(2,3)*rion(2,i)+
    >       lattice_unitg(3,3)*rion(3,i)+pi
        cw1=dcmplx(dcos(sw1),-dsin(sw1))
        cw2=dcmplx(dcos(sw2),-dsin(sw2))
        cw3=dcmplx(dcos(sw3),-dsin(sw3))
        dcpl_mb(ex1(1)+(i-1)*nx)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex2(1)+(i-1)*ny)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex2(1)+(i-1)*nz)=dcmplx(1.0d0,0.0d0)

        do k=1,nxh
          dcpl_mb(ex1(1)+k+(i-1)*nx)
     >         = dcpl_mb(ex1(1)+k-1+(i-1)*nx)*cw1
          dcpl_mb(ex1(1)+nx-k+(i-1)*nx)
     >         =dconjg(dcpl_mb(ex1(1)+k+(i-1)*nx))
        end do

        do k=1,nyh
          dcpl_mb(ex2(1)+k+(i-1)*ny)
     >         = dcpl_mb(ex2(1)+k-1+(i-1)*ny)*cw2
          dcpl_mb(ex2(1)+ny-k+(i-1)*ny)
     >         =dconjg(dcpl_mb(ex2(1)+k+(i-1)*ny))
        end do

        do k=1,nzh
          dcpl_mb(ex3(1)+k+(i-1)*nz)
     >         = dcpl_mb(ex3(1)+k-1+(i-1)*nz)*cw3
          dcpl_mb(ex3(1)+nz-k+(i-1)*nz)
     >         =dconjg(dcpl_mb(ex3(1)+k+(i-1)*nz))
        end do
        dcpl_mb(ex1(1)+nxh+(i-1)*nx)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex2(1)+nyh+(i-1)*ny)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex3(1)+nzh+(i-1)*nz)=dcmplx(0.0d0, 0.0d0)  
      end do 
      call nwpw_time_end(8)
      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ewald_strfact(ii,exi)
      implicit none
      integer ii
      complex*16 exi(*)
#include "mafdecls.fh"
#include "ewald_strfac.fh"
ccccccc locals
      integer i,j,k,indx
      integer nx,ny,nz
      integer nxh,im1
      integer taskid,np
      
      call nwpw_timing_start(8)
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nx=ewald_grid_nx()
      ny=ewald_grid_ny()
      nz=ewald_grid_nz()
      dutask=0
      im1=ii-1
      indx=1
      do k=0,nz-1
      do j=0,ny-1
      do i=0,nx-1
         if (dutask.eq.taskid) then
           exi(indx)=dcpl_mb(ex1(1)+i+im1*nx)*
    >           dcpl_mb(ex2(1)+j+im1*ny)*
    >           dcpl_mb(ex3(1)+k+im1*nz)
           indx=indx+1
         end if
         dutask=mod(dutask+1,np)
      end do
      end do
      end do
      call nwpw_timing_end(8)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccc
      

    


  
      
