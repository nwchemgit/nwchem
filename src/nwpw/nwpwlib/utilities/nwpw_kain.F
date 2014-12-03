*
*     $Id$                       
*

*     *****************************************
*     *                                       *
*     *            nwpw_kain_init             *
*     *                                       *
*     *****************************************

      subroutine nwpw_kain_init(max_m0,nsize0,rho_in)
      implicit none
      integer max_m0,nsize0
      real*8 rho_in(*)

#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i,n2ft3d


      max_m = max_m0
      nsize = nsize0

      call nwpw_list_start(1,2*max_m,nsize)

      call dcopy(max_list*max_list,0.0d0,0,B,1)
      call dcopy(max_list,0.0d0,0,d,1)
      do i=2,max_list
        B(i,1) = 1.0d0
        B(1,i) = 1.0d0
      end do
      d(1) = 1.0d0

      m = 1
      call nwpw_list_store(1,2*m,rho_in)
      return
      end


*     *****************************************
*     *                                       *
*     *              nwpw_kain                *
*     *                                       *
*     *****************************************

      subroutine nwpw_kain(rho_in,rho_out)
      implicit none
      real*8 rho_in(*)
      real*8 rho_out(*)


#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      logical value
      integer i,j
      real*8 sum, A(max_list,max_list),scal
      integer fm_ptr,x_ptr,f_ptr,dx_ptr,df_ptr
      integer dx_mptr,df_mptr
      integer ipiv(max_list),info

*     **** external functions ****
      real*8   ddot
      external ddot


      value = BA_push_get(mt_dbl,nsize,'dx_ptr',dx_mptr,dx_ptr)
      value = value.and.
     >        BA_push_get(mt_dbl,nsize,'df_ptr',df_mptr,df_ptr)
      if (.not. value) call errquit('nwpw_kain:push stack',0,0)


*     **** rho(2*m) = rho_in - rho(2*m) ****
      call nwpw_list_ptr(1,2*m,fm_ptr)
      call daxpy(nsize,(-1.0d0),rho_in,1,dbl_mb(fm_ptr),1)
      call dscal(nsize,(-1.0d0),dbl_mb(fm_ptr),1)
      call nwpw_list_store(1,2*m-1,rho_in)

*     **** compute B matrix and d ****
      do i=1,m-1
         call nwpw_list_ptr(1,2*i-1,x_ptr)

         !*** dx = x_i -x_m ****
         call dcopy(nsize,dbl_mb(x_ptr),1,  dbl_mb(dx_ptr),1)
         call daxpy(nsize,(-1.0d0),rho_in,1,dbl_mb(dx_ptr),1)

         !*** d_i = <dx_i|fm> ***
         sum = ddot(nsize,dbl_mb(dx_ptr),1,dbl_mb(fm_ptr),1)
         call D3dB_SumAll(sum)
         d(i) = sum

         do j=i,m-1
            call nwpw_list_ptr(1,2*j,f_ptr)
 
            !***  df_j = f_j - f_m ***
            call dcopy(nsize,dbl_mb(f_ptr),1,dbl_mb(df_ptr),1)
            call daxpy(nsize,(-1.0d0),dbl_mb(fm_ptr),1,dbl_mb(df_ptr),1)

            !*** B_ij = <dx_i|df_j> ***
            sum = ddot(nsize,dbl_mb(dx_ptr),1,dbl_mb(df_ptr),1)
            call D3dB_SumAll(sum)
            B(i,j) = sum
            B(j,i) = sum

         end do
      end do

*     **** solve Ac = d ****
 20   call dcopy(max_list,d,1,c,1)
      call dcopy(max_list*max_list,B,1,A,1)
      call dgesv((m-1),1,A,max_list,ipiv,c,max_list,info)

*     **** check for linear dependency ****
      write(*,*) "KAIN Check"
      do i=1,m-1
         write(*,*) "    coefficient:",i,c(i)
      end do

*     **** rho_out = rho_in + Sum(i=1,m-1)((dx_i-df_i)*c_i ) - f_m ****
      call dcopy(nsize,rho_in,1,rho_out,1)
      call daxpy(nsize,(-1.0d0),dbl_mb(fm_ptr),1,rho_out,1)
      do i=1,m-1
         call nwpw_list_ptr(1,2*i-1,x_ptr)
         call nwpw_list_ptr(1,2*i,  f_ptr)

         !*** dx = x_i -x_m ****
         call dcopy(nsize,dbl_mb(x_ptr),1,  dbl_mb(dx_ptr),1)
         call daxpy(nsize,(-1.0d0),rho_in,1,dbl_mb(dx_ptr),1)
 
         !***  df_i = f_i - f_m ***
         call dcopy(nsize,dbl_mb(f_ptr),1,dbl_mb(df_ptr),1)
         call daxpy(nsize,(-1.0d0),dbl_mb(fm_ptr),1,dbl_mb(df_ptr),1)

         call daxpy(nsize,(-1.0d0),dbl_mb(df_ptr),1,dbl_mb(dx_ptr),1)
         call daxpy(nsize,c(i),dbl_mb(dx_ptr),1,rho_out,1)

      end do

*     **** shift kain list ****
      if (m.lt.max_m) then
         m = m+1
      else
         call nwpw_kain_shift()
      end if
      call nwpw_list_store(1,2*m,rho_in)
     
      value =           BA_pop_stack(df_mptr)
      value = value.and.BA_pop_stack(dx_mptr)
      if (.not. value) call errquit('nwpw_kain:pop stack',0,0)
      return 
      end

*     *****************************************
*     *                                       *
*     *           nwpw_kain_shift             *
*     *                                       *
*     *****************************************

      subroutine nwpw_kain_shift()
      implicit none

#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i,j

      call nwpw_list_shift(1)
      call nwpw_list_shift(1)
      do j=1,m-2
      do i=1,m-2
        B(i,j) = B(i+1,j+1)
      end do
      end do

      return
      end
