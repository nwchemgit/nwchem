*
*     $Id: nwpw_diis.F,v 1.1 2003-01-24 02:07:38 bylaska Exp $                       
*

      subroutine nwpw_diis_init(max_m0,nsize0,rho_in)
      implicit none
      integer max_m0,nsize0
      real*8 rho_in(*)

#include "mafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i,n2ft3d


      max_m = max_m0
      nsize = nsize0

      call nwpw_list_start(2*max_m,nsize)

      call dcopy(max_list*max_list,0.0d0,0,B,1)
      call dcopy(max_list,0.0d0,0,d,1)
      do i=2,max_list
        B(i,1) = 1.0d0
        B(1,i) = 1.0d0
      end do
      d(1) = 1.0d0

      m = 1
      call nwpw_list_store(2*m,rho_in)
      return
      end


*     *****************************************
*     *                                       *
*     *              nwpw_diis                *
*     *                                       *
*     *****************************************

      subroutine nwpw_diis(rho_in,rho_out)
      implicit none
      real*8 rho_in(*)
      real*8 rho_out(*)


#include "mafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i,j
      real*8 sum, A(max_list,max_list),scal
      integer rr_ptr,gg_ptr,ipiv(max_list),info

*     **** external functions ****
      real*8   ddot
      external ddot


*     **** rho(2*m) = rho_in - rho(2*m) ****
      call nwpw_list_ptr(2*m,rr_ptr)
      call daxpy(nsize,(-1.0d0),rho_in,1,dbl_mb(rr_ptr),1)
      call dscal(nsize,(-1.0d0),dbl_mb(rr_ptr),1)
      call nwpw_list_store(2*m-1,rho_in)

*     **** update B matrix ****
      sum = ddot(nsize,dbl_mb(rr_ptr),1,dbl_mb(rr_ptr),1)
      call D3dB_SumAll(sum)

      B(m+1,m+1) = sum
      do i=1,m-1
         call nwpw_list_ptr(2*i,gg_ptr)
         sum = ddot(nsize,dbl_mb(gg_ptr),1,dbl_mb(rr_ptr),1)
         call D3dB_SumAll(sum)
         B(i+1,m+1) = sum
         B(m+1,i+1) = sum
      end do


*     **** solve Bc = d ****
      call dcopy(max_list,d,1,c,1)
      call dcopy(max_list*max_list,B,1,A,1)
      call dgesv(m+1,1,A,max_list,ipiv,c,max_list,info)

      write(*,*) "DIIS error history"
      sum = 0.0d0
      do i=2,m+1
         sum = sum+c(i)
         write(*,*) i," error=",B(i,i),"  coef:",c(i),sum
      end do


*     **** s = sum(i) c(i) * gg(i) ****
      call dcopy(nsize,0.0d0,0,rho_out,1)
      call daxpy(nsize,c(m+1),rho_in,1,rho_out,1)

      sum = c(m+1)
      do i=1,m-1
         sum = sum + c(i+1)
         !write(*,*) i+1,c(i+1),sum
         call nwpw_list_ptr(2*i-1, gg_ptr)
         call daxpy(nsize,c(i+1),dbl_mb(gg_ptr),1,rho_out,1)
      end do

*     **** shift B ****
      if (m.lt.max_m) then
         m = m+1
      else
         call nwpw_list_shift()
         do j=2,m
         do i=2,m
           B(i,j) = B(i+1,j+1)
         end do
         end do
      end if
      call nwpw_list_store(2*m,rho_in)
     
      return 
      end

