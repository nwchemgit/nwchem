*
*     $Id$                       
*

*     ******************************************
*     *                                        *
*     *            pspw_rmmdiis_init           *
*     *                                        *
*     ******************************************

      subroutine pspw_rmmdiis_init(max_m0)
      implicit none
      integer max_m0

#include "pspw_rmmdiis.fh"

      max_m = max_m0
      call Pack_npack(1,npack1)
      nsize = npack1
      call nwpw_list_start(2,2*max_m,2*nsize)
      return
      end

*     ******************************************
*     *                                        *
*     *            pspw_rmmdiis_end            *
*     *                                        *
*     ******************************************

      subroutine pspw_rmmdiis_end()
      implicit none
      integer max_m0

#include "pspw_rmmdiis.fh"

      call nwpw_list_end(2)
      return
      end



*     ******************************************
*     *                                        *
*     *            pspw_rmmdiis_start          *
*     *                                        *
*     ******************************************

      subroutine pspw_rmmdiis_start(lambda,residual0,psi)
      implicit none
      real*8 lambda
      complex*16 residual0(*)
      complex*16 psi(*)

#include "pspw_rmmdiis.fh"

*     **** local variables ****
      integer ngras 

      ngras = 2*npack1

      m = 1
      call nwpw_list_store(2,2*m-1,psi)
      call nwpw_list_store(2,2*m,residual0)

*     ****  psi = psi + lambda*residual0 ****
      call daxpy(ngras,lambda,residual0,1,psi,1)

      m = 2
      call nwpw_list_store(2,2*m-1,psi)

      return
      end

*     ******************************************
*     *                                        *
*     *            pspw_rmmdiis                *
*     *                                        *
*     ******************************************

      subroutine pspw_rmmdiis(lambda,residual1,psi_out)
      implicit none
      real*8     lambda
      complex*16 residual1(*)
      complex*16 psi_out(*)


#include "bafdecls.fh"
#include "pspw_rmmdiis.fh"

*     **** local variables ****
      logical value
      integer i,j,ngras,indx,indxt,info
      integer ri(2),rj(2)
      real*8  w


      ngras = 2*npack1

*     **** store residual(m) ****
      call nwpw_list_store(2,2*m,residual1)


*     **** compute B = <r(i)|r(j)> ****
      do j=1,m
      do i=1,j
        indx  = i+(j-1)*m
        indxt = j+(i-1)*m
        call nwpw_list_ptr(2,2*i,ri(1))
        call nwpw_list_ptr(2,2*j,rj(1))
        
        call Pack_cc_dot(1,dbl_mb(ri(1)),dbl_mb(rj(1)),alpha(indx))
        alpha(indxt) = alpha(indx)
      end do
      end do

*     **** compute S = <psi(i)|psi(j)> ****
      do j=1,m
      do i=1,j
        indx = i+(j-1)*m
        indx = j+(i-1)*m
        call nwpw_list_ptr(2,2*i-1,ri(1))
        call nwpw_list_ptr(2,2*j-1,rj(1))
        
        call Pack_cc_dot(1,dbl_mb(ri(1)),dbl_mb(rj(1)),S(indx))
        S(indxt) = S(indx)
      end do
      end do
  
*     **** solve B*alpha = eig*S*alpha, for alpha ****
      call dsygv(1,'V','U',m,alpha,m,S,m,eig,work,lwork,info)

c      write(*,*) "info:",info
c      do i=1,m
c        write(*,*) "eig:",i,eig(i)
c      end do

*     **** psi_out = Sum_i alpha(i)*|psi(i)> + lambda*Sum_i alpha(i)*|r(i)>  ****
      call dcopy(ngras,0.0d0,0,psi_out,1)
      do i=1,m
        call nwpw_list_ptr(2,2*i-1,ri(1))
        call nwpw_list_ptr(2,2*i,  rj(1))
c        write(*,*) "alpha:",i,alpha(i)
        
        call daxpy(ngras,       alpha(i),dbl_mb(ri(1)),1,psi_out,1)
        call daxpy(ngras,lambda*alpha(i),dbl_mb(rj(1)),1,psi_out,1)
      end do

*     **** normalize psi_out ****
      call Pack_cc_dot(1,psi_out,psi_out,w)
      w = 1.0d0/dsqrt(w)
      call dscal(ngras,w,psi_out,1)
     
      if (m.lt.max_m) then
         m = m+1
      else
         call nwpw_list_shift(2)
         call nwpw_list_shift(2)
      end if

      call nwpw_list_store(2,2*m-1,psi_out)


      return 
      end

