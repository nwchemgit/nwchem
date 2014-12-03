*
* $Id$
*


*     *********************************
*     *                               *
*     *      silvestrelli_overlap_old     *
*     *                               *
*     *********************************
*
*   This routine calculates the overlap
*
*          W = <psi_r(i)|exp(-i b*r)|psi_r(j)> 
*
*  for use in a Silvestrelli maximal Wannier function routine
*
*     Entry - b(3) : g vector
*             n    : number of wavefunctions
*             psi_r: wavefunctions
*
*     Exit - W(n,n): overlap matrix
*

      subroutine silvestrelli_overlap_old(b,n,psi_r,W)
      implicit none
      real*8     b(3)
      integer    n
      real*8     psi_r(*)
      complex*16 W(n,n)

#include "bafdecls.fh"

*     **** local variables ****
      logical value
      integer i,j,k,n2ft3d,nx,ny,nz
      real*8  br,Wc,Ws,scal1
      integer rgrid(2),Tc(2),Ts(2),psic,psis


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)

*     **** allocate space from stack ****
      call D3dB_n2ft3d(1,n2ft3d)
      value = BA_push_get(mt_dbl,(n2ft3d),'Tc',Tc(2),Tc(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'Ts',Ts(2),Ts(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(3*n2ft3d),'rgrid',rgrid(2),rgrid(1))
      psic = rgrid(1)
      psis = rgrid(1) + n2ft3d
      if (.not. value) 
     >  call errquit('silvestrelli_overlap:out of stack memory',0,0)


*     **** generate T = dcmplx(Tc,Ts) = exp(-i b*r) operator ****
      call lattice_r_grid(dbl_mb(rgrid(1)))
      do k=1,n2ft3d
         br = b(1)*dbl_mb(rgrid(1)+(k-1)*3)
     >      + b(2)*dbl_mb(rgrid(1)+(k-1)*3 + 1)
     >      + b(3)*dbl_mb(rgrid(1)+(k-1)*3 + 2)

         dbl_mb(Tc(1)+k-1) =  cos(br)
         dbl_mb(Ts(1)+k-1) = -sin(br)
         !dbl_mb(Ts(1)+k-1) =  sin(br)
      end do
      

*     **** generate W = <psi_r(i)|exp(-i b*r)|psi_r(j)> ****
      do j=1,n
      
        call D3dB_rr_Mul(1,dbl_mb(Tc(1)),
     >                     psi_r(1+(j-1)*n2ft3d),
     >                     dbl_mb(psic))
        call D3dB_rr_Mul(1,dbl_mb(Ts(1)),
     >                     psi_r(1+(j-1)*n2ft3d),
     >                     dbl_mb(psis))

        do i=j,n
          call D3dB_rr_dot(1,psi_r(1+(i-1)*n2ft3d),
     >                       dbl_mb(psic),
     >                       Wc)
          call D3dB_rr_dot(1,psi_r(1+(i-1)*n2ft3d),
     >                       dbl_mb(psis),
     >                       Ws)
          W(i,j) = dcmplx(Wc,Ws)*scal1
          W(j,i) = dcmplx(Wc,Ws)*scal1
        end do
      end do

*     **** pop memory ***
      value =           BA_pop_stack(rgrid(2))
      value = value.and.BA_pop_stack(Ts(2))
      value = value.and.BA_pop_stack(Tc(2))  
      if (.not. value) 
     >  call errquit('silvestrelli_overlap:popping stack memory',1,0)

      return
      end
