*     *************************************  
*     *                                   *
*     *        coulomb2_init              *
*     *                                   *
*     *************************************  
*
*       This routine initalizes coulomb
*                                                     
*       Entry:                                       
*       Uses:                                         
*
 
      subroutine coulomb2_init()      
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     *** cutoff constant ***
      real*8 EPSILON
      parameter (EPSILON=1.0d0)
      
c     ***** Greens functions common block *****
c     real*8     gk(ndfft3d)
      integer gk(2)
      common / vc_block2 / gk

*     **** Error function parameters ****
      real*8 c1,c2,c3,c4,c5,c6,yerf,xerf
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)

c     **** tempory data  ****
      integer glr(2)

      
c     **** local variables  ****  
      logical value
      real*8  domega,x,y,z
      real*8  dunita(3,3), dunitg(3,3) 
      real*8  pi,twopi,fourpi,temp,dV
      real*8  g1,g2,g3,gg,g
      
      integer nfft3d,n2ft3d
      integer i,j,k,i1,j1,k1,q,nq,nx,ny,nz
      integer index1,p,taskid

c     **** external functions ****
      logical  control_fast_erf
      external control_fast_erf
      real*8   lattice_unita,util_erf
      external lattice_unita,util_erf

      call D3dB_nfft3d(2,nfft3d)
      call D3dB_n2ft3d(2,n2ft3d)

*     **** allocate vc memory ****
      value = BA_alloc_get(mt_dbl,nfft3d,'gk',gk(2),gk(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      value = BA_push_get(mt_dbl,n2ft3d,'glr',glr(2),glr(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      call dcopy(nfft3d,0.0d0,0,dbl_mb(gk(1)),1)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(glr(1)),1)


      call Parallel2d_taskid_i(taskid)
c      call D3dB_nq(2,nq)
      call D3dB_nz(2,nz)
      call D3dB_ny(2,ny)
      call D3dB_nx(2,nx)

      pi     = 4.0d0*datan(1.0d0)
      twopi  = 2.0d0*pi
      fourpi = 4.0d0*pi

  
*     **** define lattice on expanded grid ****
      do j=1,3
        do i=1,3
           dunita(i,j) = 2.0d0*lattice_unita(i,j)
        end do
      end do
    

*     **** reciprical vectors for expanded grid ****
      dunitg(1,1) = dunita(2,2)*dunita(3,3) - dunita(3,2)*dunita(2,3)
      dunitg(2,1) = dunita(3,2)*dunita(1,3) - dunita(1,2)*dunita(3,3)
      dunitg(3,1) = dunita(1,2)*dunita(2,3) - dunita(2,2)*dunita(1,3)
      dunitg(1,2) = dunita(2,3)*dunita(3,1) - dunita(3,3)*dunita(2,1)
      dunitg(2,2) = dunita(3,3)*dunita(1,1) - dunita(1,3)*dunita(3,1)
      dunitg(3,2) = dunita(1,3)*dunita(2,1) - dunita(2,3)*dunita(1,1)
      dunitg(1,3) = dunita(2,1)*dunita(3,2) - dunita(3,1)*dunita(2,2)
      dunitg(2,3) = dunita(3,1)*dunita(1,2) - dunita(1,1)*dunita(3,2)
      dunitg(3,3) = dunita(1,1)*dunita(2,2) - dunita(2,1)*dunita(1,2)
                   
      domega = dunita(1,1)*dunitg(1,1)
     &       + dunita(2,1)*dunitg(2,1)
     &       + dunita(3,1)*dunitg(3,1)
      call dscal(9,twopi/domega,dunitg,1)


*     **** volume of a unit cell ***
      domega = dabs(domega)
      dV     = domega/dble(nx*ny*nz)


*     ***********************************
c     *** Set up the Green's function ***
*     ***********************************

*     **** short-range part of Greens function ****
      call dcopy(nfft3d,0.0d0,0,dbl_mb(gk(1)),1)
      do k = 0, nz-1
         do j = 0, ny-1
         do i = 0, (nx/2)
            i1  = i
            j1  = j
            k1  = k
           
            call D3dB_ijktoindexp(2,i+1,j+1,k+1,index1,p)
            if (p.eq.taskid) then
            if ((i1+j1+k1) .gt. 0) then
                  
               if (j1 .gt. ny/2) j1 = j1-ny
               if (k1 .gt. nz/2) k1 = k1-nz
             
               g1 = i1*dunitg(1,1) + j1*dunitg(1,2) + k1*dunitg(1,3)
               g2 = i1*dunitg(2,1) + j1*dunitg(2,2) + k1*dunitg(2,3)
               g3 = i1*dunitg(3,1) + j1*dunitg(3,2) + k1*dunitg(3,3)
               gg = g1*g1 + g2*g2 + g3*g3
               g  = dsqrt(gg)
         
               temp = (fourpi/gg)
     >              * (1.0d0 - dexp(-gg/(4.0d0*EPSILON**2)))
            else
               temp = pi/EPSILON**2
            end if

c           index1     = (i+1) 
c    >                 + j*(nx/2+1) 
c    >                 + (q-1)*(nx/2+1)*ny
            dbl_mb(gk(1)+index1-1) = temp

            end if
         end do
         end do
      end do

*     **** long-range part of Greens function ****
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(glr(1)),1)

      if (control_fast_erf()) then

      do k = 1, nz
         k1=k-1
         if (k1 .gt. nz/2) k1 = k1-nz

         do j = 1, ny
            j1=j-1
            if (j1 .gt. ny/2) j1 = j1-ny

            do i = 1, nx
               i1=i-1
               if (i1 .gt. nx/2) i1 = i1-nx

               call D3dB_ijktoindex2p(2,i,j,k,index1,p)
               if (p.eq.taskid) then
               x = i1*dunita(1,1)/dble(nx)
     >           + j1*dunita(1,2)/dble(ny)
     >           + k1*dunita(1,3)/dble(nz)
               y = i1*dunita(2,1)/dble(nx)
     >           + j1*dunita(2,2)/dble(ny)
     >           + k1*dunita(2,3)/dble(nz)
               z = i1*dunita(3,1)/dble(nx)
     >           + j1*dunita(3,2)/dble(ny)
     >           + k1*dunita(3,3)/dble(nz)

c              x = i1*dunita(1,1)/dble(nx)
c    >           + k1*dunita(1,2)/dble(ny)
c    >           + j1*dunita(1,3)/dble(nz)
c              y = i1*dunita(2,1)/dble(nx)
c    >           + k1*dunita(2,2)/dble(ny)
c    >           + j1*dunita(2,3)/dble(nz)
c              z = i1*dunita(3,1)/dble(nx)
c    >           + k1*dunita(3,2)/dble(ny)
c    >           + j1*dunita(3,3)/dble(nz)


               temp = dsqrt(x*x + y*y + z*z)
               if (temp .gt. 1.0d-15) then

                  xerf=EPSILON*temp
                  yerf = (1.0d0
     >                 + xerf*(c1 + xerf*(c2
     >                 + xerf*(c3 + xerf*(c4
     >                 + xerf*(c5 + xerf*c6))))))**4
                  yerf = (1.0d0 - 1.0d0/yerf**4)
                  temp = yerf/temp

c                  temp = util_erf(EPSILON*temp)/temp
               else
                  temp = 2.0d0*EPSILON/dsqrt(PI)
               end if

               !index1 = i + (j-1)*(nx+2) + (q-1)*(nx+2)*ny
               dbl_mb(glr(1)+index1-1) = temp*dV

               end if
            end do
         end do
      end do


      else

      do k = 1, nz
         k1=k-1
         if (k1 .gt. nz/2) k1 = k1-nz

         do j = 1, ny
            j1=j-1
            if (j1 .gt. ny/2) j1 = j1-ny

            do i = 1, nx
               i1=i-1
               if (i1 .gt. nx/2) i1 = i1-nx

               call D3dB_ijktoindex2p(2,i,j,k,index1,p)
               if (p.eq.taskid) then
               x = i1*dunita(1,1)/dble(nx)
     >           + j1*dunita(1,2)/dble(ny)
     >           + k1*dunita(1,3)/dble(nz)
               y = i1*dunita(2,1)/dble(nx)
     >           + j1*dunita(2,2)/dble(ny)
     >           + k1*dunita(2,3)/dble(nz)
               z = i1*dunita(3,1)/dble(nx)
     >           + j1*dunita(3,2)/dble(ny)
     >           + k1*dunita(3,3)/dble(nz)

c              x = i1*dunita(1,1)/dble(nx)
c    >           + k1*dunita(1,2)/dble(ny)
c    >           + j1*dunita(1,3)/dble(nz)
c              y = i1*dunita(2,1)/dble(nx)
c    >           + k1*dunita(2,2)/dble(ny)
c    >           + j1*dunita(2,3)/dble(nz)
c              z = i1*dunita(3,1)/dble(nx)
c    >           + k1*dunita(3,2)/dble(ny)
c    >           + j1*dunita(3,3)/dble(nz)

               temp = dsqrt(x*x + y*y + z*z)
               if (temp .gt. 1.0d-15) then

c                 xerf=EPSILON*temp
c                 yerf = (1.0d0
c    >                 + xerf*(c1 + xerf*(c2
c    >                 + xerf*(c3 + xerf*(c4
c    >                 + xerf*(c5 + xerf*c6))))))**4
c                 yerf = (1.0d0 - 1.0d0/yerf**4)
c                 temp = yerf/temp

                  temp = util_erf(EPSILON*temp)/temp
               else
                  temp = 2.0d0*EPSILON/dsqrt(PI)
               end if

               !index1 = i + (j-1)*(nx+2) + (q-1)*(nx+2)*ny
               dbl_mb(glr(1)+index1-1) = temp*dV

               end if
            end do
         end do
      end do

      end if

*     **** Put glr in k-space ****
      call D3dB_rc_fft3f(2,dbl_mb(glr(1)))
        

*     **** add long-range part to short-range part             ****
*     *** note that only real parts of tranformed grl are used ****
      do k=1,nfft3d
         dbl_mb(gk(1)+k-1) = dbl_mb(gk(1)+k-1) 
     >                     + dbl_mb(glr(1)+2*(k-1)) 
      end do

      if (.not.BA_pop_stack(glr(2)))
     > call errquit('coulomb2_init:error popping stack memory',0,MA_ERR)

      return
      end


      subroutine coulomb2_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      integer gk(2)
      common / vc_block2 / gk

      if (.not.BA_free_heap(gk(2)))
     > call errquit('coulomb2_end:error freeing heap memory',0,MA_ERR)

      return
      end




*     *****************************
*     *                           *
*     *        coulomb2_v          *
*     *                           *
*     *****************************
*
*  This routine calculates Poisson's equation for infinite  
*  space boundry conditions.                               
*                                                              
*       Laplacian(vh) = -4*pi*dn                      
*                                                            
*       vh(r-->infinity) = 0                
*                                                          
*                                                         
*   Entry:                                           
*         dn --- the density of the region        
*   Exit:                                         
*         vh --- the solution to Poisson's equation 

      subroutine coulomb2_v(dn, vh)
      implicit none
#include "errquit.fh"
      real*8 dn(*),vh(*)
            
#include "bafdecls.fh"

*     **** common block used for coulomb.f ****
      integer gk(2)
      common / vc_block2 / gk

*     **** local variables ****
      logical value
      integer tmpx(2)
      integer nx,ny,nz,nfft3d
      real*8 scale
                  
      call D3dB_nfft3d(2,nfft3d)
      value = BA_push_get(mt_dcpl,nfft3d,'tmpx',tmpx(2),tmpx(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(tmpx(1)),1)

      call D3dB_nx(2,nx)
      call D3dB_ny(2,ny)
      call D3dB_nz(2,nz)
      scale = 1.0d0/dble(nx*ny*nz)

c     *** Expand the density ***
      call D3dB_r_Expand(1,dn,2,dcpl_mb(tmpx(1)))

c     **** Convolution g*dn ****
      call D3dB_rc_fft3f(2,dcpl_mb(tmpx(1)))
      call D3dB_tc_Mul2(2,dbl_mb(gk(1)),dcpl_mb(tmpx(1)))
      call D3dB_cr_fft3b(2,dcpl_mb(tmpx(1)))

c     **** contract tmpx to extract vh ****
      call D3dB_r_Contract(2,dcpl_mb(tmpx(1)),1,vh)
c      call D3dB_r_SMul(1,scale,vh,vh)
      call D3dB_r_SMul1(1,scale,vh)
           
      value = BA_pop_stack(tmpx(2))
      if (.not. value) call errquit('error popping stack memory',0,
     &       MA_ERR)
      return
      end
c $Id$
