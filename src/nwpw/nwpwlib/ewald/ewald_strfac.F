*
* $Id$
*
      subroutine ewald_strfact_init()
      implicit none
#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald_strfac.fh"
cccccccc arguments
cccccccc locals
      integer nion
      logical value
cccccccc externals
      integer  ion_nion,ewald_npack
      external ion_nion,ewald_npack
      integer  ewald_grid_nx,ewald_grid_ny,ewald_grid_nz
      external ewald_grid_nx,ewald_grid_ny,ewald_grid_nz
cccccccc exec
      nion  = ion_nion()

      npack = ewald_npack()
      nx    = ewald_grid_nx()
      ny    = ewald_grid_ny()
      nz    = ewald_grid_nz()
      notzero_npack = (npack.gt.0)

      value=BA_alloc_get(mt_dcpl,(nx*nion),'ewx1',ewx1(2),ewx1(1))
      value=value.and.
     >      BA_alloc_get(mt_dcpl,(ny*nion),'ewx2',ewx2(2),ewx2(1))
      value=value.and.
     >      BA_alloc_get(mt_dcpl,(nz*nion),'ewx3',ewx3(2),ewx3(1))

      if (notzero_npack) then
         value=BA_alloc_get(mt_int,npack,'i_indxe',i_indx(2),i_indx(1))
         value=value.and.
     >         BA_alloc_get(mt_int,npack,'j_indxe',j_indx(2),j_indx(1))
         value=value.and.
     >         BA_alloc_get(mt_int,npack,'k_indxe',k_indx(2),k_indx(1))
      end if
      if (.not.value) 
     >  call errquit("EWALD_STRUCTFACTOR_INIT: OUT OF HEAP MEMORY",0,
     >               MA_ERR)
      return 
      end
cccccccccccccccccccccccc
      subroutine ewald_strfact_end()
      implicit none
#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald_strfac.fh"
      logical value
cccccccc exec
      value=          BA_free_heap(ewx1(2))
      value=value.and.BA_free_heap(ewx2(2))
      value=value.and.BA_free_heap(ewx3(2))

      if (notzero_npack) then
         value=value.and.BA_free_heap(i_indx(2))
         value=value.and.BA_free_heap(j_indx(2))
         value=value.and.BA_free_heap(k_indx(2))
      end if
      if (.not.value) then
        call errquit("ewald_strfact_end: error free heap",0,MA_ERR)
      end if
      return
      end
ccccccccccccccccccccccccccccc
      integer function ewald_strfact_i_indx()
      implicit none
#include "ewald_strfac.fh"
      ewald_strfact_i_indx = i_indx(1)
      return
      end 
      integer function ewald_strfact_j_indx()
      implicit none
#include "ewald_strfac.fh"
      ewald_strfact_j_indx = j_indx(1)
      return
      end 
      integer function ewald_strfact_k_indx()
      implicit none
#include "ewald_strfac.fh"
      ewald_strfact_k_indx = k_indx(1)
      return
      end 
ccccccccccccccccccccccccccccc
      subroutine ewald_phafac()
      implicit none
#include "bafdecls.fh"
#include "ewald_strfac.fh"
cccccccc locals
      integer i,k
      integer nxh,nyh,nzh
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

ccccccccc openmp
      integer tid,nthr
      integer Parallel_threadid, Parallel_nthreads
      external Parallel_threadid, Parallel_nthreads

ccccccc start of code       
      call nwpw_timing_start(8)

      tid = Parallel_threadid()
      nthr = Parallel_nthreads()

      pi=4.0d0*datan(1.0d0)
      nion=ion_nion()

      nxh=nx/2
      nyh=ny/2
      nzh=nz/2
c     Each thread processes independend ions (if USE_OPENMP is defined)
      do i=tid+1,nion,nthr
        sw1=lattice_unitg(1,1)*ion_rion(1,i) 
     >     +lattice_unitg(2,1)*ion_rion(2,i) 
     >     +lattice_unitg(3,1)*ion_rion(3,i)+pi

        sw2=lattice_unitg(1,2)*ion_rion(1,i)
     >     +lattice_unitg(2,2)*ion_rion(2,i)
     >     +lattice_unitg(3,2)*ion_rion(3,i)+pi

        sw3=lattice_unitg(1,3)*ion_rion(1,i)
     >     +lattice_unitg(2,3)*ion_rion(2,i)
     >     +lattice_unitg(3,3)*ion_rion(3,i)+pi

        cw1=dcmplx(dcos(sw1),-dsin(sw1))
        cw2=dcmplx(dcos(sw2),-dsin(sw2))
        cw3=dcmplx(dcos(sw3),-dsin(sw3))
        dcpl_mb(ewx1(1)+(i-1)*nx)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ewx2(1)+(i-1)*ny)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ewx3(1)+(i-1)*nz)=dcmplx(1.0d0,0.0d0)

        do k=1,nxh
          dcpl_mb(ewx1(1)+k+(i-1)*nx)
     >         =dcpl_mb(ewx1(1)+k-1+(i-1)*nx)*cw1
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
#include "bafdecls.fh"
#include "ewald_strfac.fh"
ccccccc locals

      call nwpw_timing_start(8)
      if (notzero_npack) then
         call ewald_strfac_sub(npack,
     >                     int_mb(i_indx(1)),
     >                     int_mb(j_indx(1)),
     >                     int_mb(k_indx(1)),
     >                     dcpl_mb(ewx1(1)+(ii-1)*nx),
     >                     dcpl_mb(ewx2(1)+(ii-1)*ny),
     >                     dcpl_mb(ewx3(1)+(ii-1)*nz),
     >                     exi)
      end if
      call nwpw_timing_end(8)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ewald_strfac_sub(npack,
     >                            indxi,indxj,indxk,
     >                            ex1,ex2,ex3,
     >                            exi)
      implicit none
      integer npack,indxi(*),indxj(*),indxk(*)
      complex*16 ex1(*),ex2(*),ex3(*),exi(*)
      integer indx
!$OMP DO
      do indx=1,npack
        exi(indx) = ex1(indxi(indx))*ex2(indxj(indx))*ex3(indxk(indx))
      end do
!$OMP END DO
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccc
      

    


  
      
