*
* $Id: integrate.F,v 1.4 2003-02-19 01:11:22 marat Exp $
*
!*************************************************
!
!   Name    : def_integr
!
!
!   Purpose :
!
!
!   Created :
!
!*************************************************
       function c_def_integr(power_f,
     >                           f,
     >                           power_r,
     >                           r,
     >                           log_amesh,
     >                           nrange)
                    
      implicit none
      integer power_f
      complex*16  f(*)
      integer power_r
      real*8 r(*)
      real*8 log_amesh
      integer nrange
      complex*16 c_def_integr

#include "mafdecls.fh"


      !*** local variables ***
      integer i
      logical value
      integer k
      complex*16 sum_f
      complex*16 tmp_sum
      integer integrand(2)


      value = MA_push_get(mt_dcpl,nrange,
     >                    'integrand',integrand(2),integrand(1))
      if (.not.value) call errquit('def_integr: push stack',0)

      do k=1,nrange
        dcpl_mb(integrand(1)+k-1) = f(k)*r(k)**(power_r+1)
      end do

      !*** integrate from the origin to the first point ***
      sum_f = dcpl_mb(integrand(1))/(power_r+power_f+1)

      !*** the rest via trapesoidal rule ***
      tmp_sum = dcmplx(0.0,0.0)
      do i=1,nrange
       tmp_sum = tmp_sum + dcpl_mb(integrand(1)+i-1)
      end do
      
      sum_f = sum_f 
c     >      + log_amesh*zsum(nrange,dcpl_mb(integrand(1)),1) 
     >      + log_amesh*tmp_sum
     >      - 0.5d0*log_amesh
     >        *(dcpl_mb(integrand(1))
     >         +dcpl_mb(integrand(1)+nrange-1))

      value = MA_pop_stack(integrand(2))
      if (.not.value) call errquit('def_integr: pop stack',1)

      c_def_integr = sum_f

      return
      end !FUNCTION def_integr



!*************************************************
!
!   Name    : def_integr
!
!
!   Purpose :
!
!
!   Created :
!
!*************************************************
      real*8 function def_integr(power_f,
     >                           f,
     >                           power_r,
     >                           r,
     >                           log_amesh,
     >                           nrange)
                    
      implicit none
      integer power_f
      real*8  f(*)
      integer power_r
      real*8 r(*)
      real*8 log_amesh
      integer i,nrange

#include "mafdecls.fh"

      !*** local variables ***
      logical value
      integer k
      real*8 sum_f,tmp_sum
      integer integrand(2)

      !*** external functions ***
      real*8   dsum
      external dsum

      value = MA_push_get(mt_dbl,nrange,
     >                    'integrand',integrand(2),integrand(1))
      if (.not.value) call errquit('def_integr: push stack',0)

      do k=1,nrange
        dbl_mb(integrand(1)+k-1) = f(k)*r(k)**(power_r+1)
      end do

      !*** integrate from the origin to the first point ***
      sum_f = dbl_mb(integrand(1))/(power_r+power_f+1)
      
      !*** the rest via trapesoidal rule ***
      tmp_sum = 0.0d0
      do i=1,nrange
       tmp_sum = tmp_sum + dbl_mb(integrand(1)+i-1)
      end do
 
      !*** the rest via trapesoidal rule ***
      sum_f = sum_f 
c     >      + log_amesh*dsum(nrange,dbl_mb(integrand(1)),1) 
     >      + log_amesh*tmp_sum
     >      - 0.5d0*log_amesh
     >        *(dbl_mb(integrand(1))
     >         +dbl_mb(integrand(1)+nrange-1))

      value = MA_pop_stack(integrand(2))
      if (.not.value) call errquit('def_integr: pop stack',1)

      def_integr = sum_f

      return
      end !FUNCTION def_integr


!*************************************************
!
!   Name    : indef_integr
!
!
!   Purpose :
!
!
!   Created :
!
!*************************************************
      subroutine indef_integr(power_f,
     >                        f,
     >                        power_r,
     >                        r,
     >                        log_amesh,
     >                        nrange,
     >                        sum_f)
      implicit none
      integer  power_f
      real*8   f(*)
      integer  power_r
      real*8  r(*)
      real*8  log_amesh
      integer nrange
      real*8  sum_f(*)
 
#include "mafdecls.fh"

      !*** local variables ***
      logical value
      integer k
      integer integrand(2)

      value = MA_push_get(mt_dbl,nrange,
     >                    'integrand',integrand(2),integrand(1))
      if (.not.value) call errquit('indef_integr: push stack',0)

      do k=1,nrange
       dbl_mb(integrand(1)+k-1) = f(k)*r(k)**(power_r+1)
      end do

      if(nrange.le.5) then
         do k=1,nrange 
           sum_f(k) = dbl_mb(integrand(1)+k-1)/(power_r+power_f+1)
         end do
      else
         do k=1,5
           sum_f(k) = dbl_mb(integrand(1)+k-1)/(power_r+power_f+1)
         end do
         do k=6,nrange
          sum_f(k) = sum_f(k-1) 
     >             + log_amesh*0.5d0*( dbl_mb(integrand(1)+k-2)
     >                             +   dbl_mb(integrand(1)+k-1))
         end do

      end if

      value = MA_pop_stack(integrand(2))
      if (.not.value) call errquit('indef_integr: pop stack',0)


      return
      end !SUBROUTINE indef_integr


