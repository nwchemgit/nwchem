*
*     $Id$                       
*

*     *****************************************
*     *                                       *
*     *              nwpw_diis_init           *
*     *                                       *
*     *****************************************

      subroutine nwpw_diis_init(max_m0,nsize0,rho_in)
      implicit none
      integer max_m0,nsize0
      real*8 rho_in(*)

#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i

      max_m = max_m0
      nsize = nsize0

      call nwpw_list_start(1,2*max_m,nsize)
      call nwpw_diis_reset(rho_in)

      return
      end


*     *****************************************
*     *                                       *
*     *              nwpw_diis_reset          *
*     *                                       *
*     *****************************************

      subroutine nwpw_diis_reset(rho_in)
      implicit none
      real*8 rho_in(*)

#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i

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
*     *              nwpw_diis                *
*     *                                       *
*     *****************************************

      subroutine nwpw_diis(rho_in,v_in,v_out,diis_error)
      implicit none
      real*8 rho_in(*)
      real*8 v_in(*)
      real*8 v_out(*)
      real*8 diis_error


#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i,j,iwork(max_list)
      real*8 sum, scal,ferr,berr,rcond
      real*8 AP(max_list*(max_list+1)/2)
      real*8 AFP(max_list*(max_list+1)/2)
      real*8 work(3*max_list)
      real*8 A2(max_list,max_list),eig(max_list)
      integer rr_ptr,gg_ptr,ipiv(max_list),info
      integer rr1,rr2,rr3

*     **** external functions ****
      real*8   ddot
      external ddot


*     **** rho(2*m) = rho_in - rho(2*m) ****
      call nwpw_list_ptr(1,2*m,rr_ptr)
      call daxpy(nsize,(-1.0d0),rho_in,1,dbl_mb(rr_ptr),1)
      call dscal(nsize,(-1.0d0),dbl_mb(rr_ptr),1)


*     **** store v_in in rho(2*m-1) ***
      call nwpw_list_store(1,2*m-1,v_in)

*     **** update B matrix ****
      sum = ddot(nsize,dbl_mb(rr_ptr),1,dbl_mb(rr_ptr),1)
      call D3dB_SumAll(sum)

      B(m+1,m+1) = sum
      diis_error = sum

      do i=1,m-1
         call nwpw_list_ptr(1,2*i,gg_ptr)
         sum = ddot(nsize,dbl_mb(gg_ptr),1,dbl_mb(rr_ptr),1)
         call D3dB_SumAll(sum)
         B(i+1,m+1) = sum
         B(m+1,i+1) = sum
      end do

c     !*** swap vectors and residuals ***
c     if ((B(m+1,m+1).gt.B(m,m)).and.
c    >   (B(m,m)     .gt.0.0d0)) then
c        call nwpw_diis_swap(m,m-1)
c     end if


c DEBUG
c     !*** look at eigenvalues of B*B ***
c     call dcopy(max_list*max_list,0.0d0,1,A2,1)
c     call dcopy(max_list*max_list,0.0d0,1,A,1)
c     do j=1,m
c     do i=1,m
c       A(i,j) = B(i+1,j+1)
c     end do
c     end do
c     call DMMUL(max_list,m,A,A,A2)
c     call eigen_unsort(max_list,m,A2,eig,A,info)
c     write(*,*)
c     write(*,*)
c     do i=1,m
c       write(*,*) "   eigenvalue:",i,eig(i)
c     end do
c     write(*,*)
c     write(*,*)
c
c     !*** swap vectors and residuals using Bubble Sort ***
c     do i=m,2,-1
c     do j=i-1,1,-1
c     if (eig(j).gt.eig(i)) then
c        call nwpw_diis_swap(i,j)
c        sum    = eig(i)
c        eig(i) = eig(j)
c        eig(j) = sum
c     end if
c     end do
c     end do
c     write(*,*)
c     write(*,*)
c     do i=1,m
c       write(*,*) "  seigenvalue:",i,eig(i)
c     end do
c     write(*,*)
c     write(*,*)
c     
c     !*** remove small vectors ****
c     j = m-2 
c     do i=1,j
c       if ((eig(i)/eig(m)).lt.1.0d-8) then
c        call nwpw_diis_shift()
c        m = m-1
c       end if
c     end do
c DEBUG



*     **** solve Bc = d ****
 20   continue
      call dcopy(max_list,d,1,c,1)

      !call dcopy(max_list*max_list,B,1,A,1)
      !call dgesv(m+1,1,A,max_list,ipiv,c,max_list,info)
      do j=1,m+1
      do i=1,j
        AP(i+(j-1)*j/2) = B(i,j)
      end do
      end do
      call dspsvx( 'N', 'U', m+1, 1, AP, AFP,
     >             ipiv, d,  m+1, c, m+1,
     >             rcond, ferr, berr, 
     >             work,iwork,
     >             info)

*     **** check for linear dependency ****
c     sum = 0.0d0
c     do i=2,m+1
c        sum = sum+dabs(c(i))
c     end do
c     write(*,*) "DEBUG: sum=",sum,"  count=",(m+1),
c    >           "   info=",info,"   rcond=",rcond

      !if (sum.gt.1.2d0) then
      !   call nwpw_diis_gram_schmidt() !*** this routine does not work ***
      !   go to 20
      !end if

      if ((info.ne.0).or.(rcond.lt.1.0d-10)) then
         call nwpw_diis_shift()
         m = m-1
         go to 20
      end if


*     **** s = sum(i) c(i) * gg(i) ****
      call dcopy(nsize,0.0d0,0,v_out,1)
      call daxpy(nsize,c(m+1),v_in,1,v_out,1)

      sum = c(m+1)
      do i=1,m-1
         sum = sum + c(i+1)
         !write(*,*) i+1,c(i+1),sum
         call nwpw_list_ptr(1,2*i-1, gg_ptr)
         call daxpy(nsize,c(i+1),dbl_mb(gg_ptr),1,v_out,1)
      end do

*     **** shift B and diis list ****
      if (m.lt.max_m) then
         m = m+1
      else
         call nwpw_diis_shift()
      end if
      call nwpw_list_store(1,2*m,rho_in)
     
      return 
      end

*     *****************************************
*     *                                       *
*     *           nwpw_diis_shift             *
*     *                                       *
*     *****************************************

      subroutine nwpw_diis_shift()
      implicit none

#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i,j

      call nwpw_list_shift(1)
      call nwpw_list_shift(1)
      do j=2,m
      do i=2,m
        B(i,j) = B(i+1,j+1)
      end do
      end do

      return
      end

*     *****************************************
*     *                                       *
*     *           nwpw_diis_remove            *
*     *                                       *
*     *****************************************

      subroutine nwpw_diis_remove(k)
      implicit none
      integer k

#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i,j

      !*** remove k residual and density ***
      call nwpw_list_remove(1,2*k)
      call nwpw_list_remove(1,2*k-1)

      !**** remove k row ****
      do i=k,m
      do j=2,m+1
        B(i,j) = B(i+1,j)
      end do
      end do

      !**** remove k column  ****
      do j=k,m
      do i=2,m
        B(i,j) = B(i,j+1)
      end do
      end do

      m = m-1

      return
      end

*     *****************************************
*     *                                       *
*     *           nwpw_diis_swap              *
*     *                                       *
*     *****************************************

      subroutine nwpw_diis_swap(k,l)
      implicit none
      integer k,l

#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i,j
      real*8  tmp

      !*** swap  k and l residual and density ***
      call nwpw_list_swap(1,2*k,  2*l)
      call nwpw_list_swap(1,2*k-1,2*l-1)


      !**** swap k and l rows ****
      do j=2,m+1
        tmp      = B(l+1,j)
        B(l+1,j) = B(k+1,j)
        B(k+1,j) = tmp
      end do

      !**** swap k and l columns ****
      do i=2,m+1
        tmp      = B(i,l+1)
        B(i,l+1) = B(i,k+1)
        B(i,k+1) = tmp
      end do


      return
      end


*     *****************************************
*     *                                       *
*     *           nwpw_diis_gram_schmidt      *
*     *                                       *
*     *****************************************

      subroutine nwpw_diis_gram_schmidt()
      implicit none
      
#include "bafdecls.fh"
#include "nwpw_diis_common.fh"

*     **** local variables ****
      integer i,j,n
      integer ri_ptr,rj_ptr,vi_ptr,vj_ptr
      real*8  sum,sum1

*     **** external functions ****
      real*8   ddot
      external ddot

      !**** Gram-Schmidt ****
      do i=m,2,-1
        call nwpw_list_ptr(1,2*i,  ri_ptr)
        call nwpw_list_ptr(1,2*i-1,vi_ptr)
        do j= i-1, 1,-1
          call nwpw_list_ptr(1,2*j,  rj_ptr)
          call nwpw_list_ptr(1,2*j-1,vj_ptr)
          sum  = ddot(nsize,dbl_mb(rj_ptr),1,dbl_mb(ri_ptr),1)
          sum1 = ddot(nsize,dbl_mb(rj_ptr),1,dbl_mb(rj_ptr),1)
          call D3dB_SumAll(sum)
          call D3dB_SumAll(sum1)
          call daxpy(nsize,(-sum/sum1),
     >              dbl_mb(rj_ptr),1,
     >              dbl_mb(ri_ptr),1)
          call daxpy(nsize,(-sum/sum1),
     >              dbl_mb(vj_ptr),1,
     >              dbl_mb(vi_ptr),1)

        end do
      end do

      !*** remove small residuals ***
      n = m
      do i=n-1,1,-1
        call nwpw_list_ptr(1,2*i,  ri_ptr)
        sum1 = ddot(nsize,dbl_mb(ri_ptr),1,dbl_mb(ri_ptr),1)
        call D3dB_SumAll(sum1)
        write(*,*) "sum1:",i,sum1
        if (sum1.lt.1.0d-3) then
           call nwpw_list_remove(1,2*i)
           call nwpw_list_remove(1,2*i-1)
           m = m - 1
        end if
      end do

      !*** regenerate B ***
      do j=1,m
         call nwpw_list_ptr(1,2*j,rj_ptr)
         do i=j,m
           call nwpw_list_ptr(1,2*i,ri_ptr)
           sum = ddot(nsize,dbl_mb(ri_ptr),1,dbl_mb(rj_ptr),1)
           call D3dB_SumAll(sum)
           B(i+1,j+1) = sum
           B(j+1,i+1) = sum
         end do
      end do

      return
      end 
