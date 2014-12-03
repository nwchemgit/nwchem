*
* $Id$
*


*     *********************************
*     *                               *
*     *      silvestrelli_overlap     *
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
*             psi_r2: tmp wavefunctions
*
*     Exit - W(n,n): overlap matrix
*

      subroutine silvestrelli_overlap(b,ms,ne,neq,psi_r,psi_r2,W)
      implicit none
      real*8     b(3)
      integer    ms,ne(2),neq(2)
      real*8     psi_r(*)
      real*8     psi_r2(*)
      complex*16 W(*)

#include "bafdecls.fh"

*     **** local variables ****
      logical value
      integer i,j,k,n2ft3d,nx,ny,nz
      !real*8  br,Wc,Ws,scal1
      real*8  br,scal1
      integer rgrid(2),Tc(2),Ts(2),psic,psis,Wc(2),Ws(2)

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

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
      value = value.and.Dneall_m_push_get(ms,Wc)
      value = value.and.Dneall_m_push_get(ms,Ws)
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
      do j=1,neq(ms)
        call D3dB_rr_Mul(1,dbl_mb(Tc(1)),
     >                       psi_r(1+(j-1+(ms-1)*neq(1))*n2ft3d),
     >                      psi_r2(1+(j-1+(ms-1)*neq(1))*n2ft3d))
      end do
      call Dneall_ggm_sym_Multiply(ms,psi_r,psi_r2,n2ft3d,dbl_mb(Wc(1)))
      call Dneall_m_scal(ms,scal1,dbl_mb(Wc(1)))
      do j=1,neq(ms)
        call D3dB_rr_Mul(1,dbl_mb(Ts(1)),
     >                       psi_r(1+(j-1+(ms-1)*neq(1))*n2ft3d),
     >                      psi_r2(1+(j-1+(ms-1)*neq(1))*n2ft3d))
      end do
      call Dneall_ggm_sym_Multiply(ms,psi_r,psi_r2,n2ft3d,dbl_mb(Ws(1)))
      call Dneall_m_scal(ms,scal1,dbl_mb(Ws(1)))
      call Dneall_mmtow_Cmplx(ms,dbl_mb(Wc(1)),dbl_mb(Ws(1)),W)

c      do j=1,n
c      
c        call D3dB_rr_Mul(1,dbl_mb(Tc(1)),
c     >                      psi_r(1+(j-1)*n2ft3d),
c     >                      dbl_mb(psic))
c        call D3dB_rr_Mul(1,dbl_mb(Ts(1)),
c     >                      psi_r(1+(j-1)*n2ft3d),
c     >                      dbl_mb(psis))
c
c        do i=j,n
c          call D3dB_rr_dot(1,psi_r(1+(i-1)*n2ft3d),
c     >                       dbl_mb(psic),
c     >                       Wc)
c          call D3dB_rr_dot(1,psi_r(1+(i-1)*n2ft3d),
c     >                       dbl_mb(psis),
c     >                       Ws)
c          W(i,j) = dcmplx(Wc,Ws)*scal1
c          W(j,i) = dcmplx(Wc,Ws)*scal1
c        end do
c      end do

*     **** pop memory ***
      value =           Dneall_m_pop_stack(Ws)
      value = value.and.Dneall_m_pop_stack(Wc)
      value = value.and.BA_pop_stack(rgrid(2))
      value = value.and.BA_pop_stack(Ts(2))
      value = value.and.BA_pop_stack(Tc(2))  
      if (.not. value) 
     >  call errquit('silvestrelli_overlap:popping stack memory',1,0)

      return
      end
