*
* $Id$
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

#include "bafdecls.fh"
#include "errquit.fh"


      !*** local variables ***
      integer i
      logical value
      integer k
      complex*16 sum_f
      complex*16 tmp_sum
      integer integrand(2)


      value = BA_push_get(mt_dcpl,nrange,
     >                    'integrand',integrand(2),integrand(1))
      if (.not.value) call errquit('def_integr: push stack',0, MA_ERR)

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

      value = BA_pop_stack(integrand(2))
      if (.not.value) call errquit('def_integr: pop stack',1, MA_ERR)

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

#include "bafdecls.fh"
#include "errquit.fh"

      !*** local variables ***
      logical value
      integer k
      real*8 sum_f,tmp_sum
      integer integrand(2)

      !*** external functions ***
      real*8   dsum
      external dsum

      value = BA_push_get(mt_dbl,nrange,
     >                    'integrand',integrand(2),integrand(1))
      if (.not.value) call errquit('def_integr: push stack',0, MA_ERR)

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

      value = BA_pop_stack(integrand(2))
      if (.not.value) call errquit('def_integr: pop stack',1, MA_ERR)

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
 
#include "bafdecls.fh"
#include "errquit.fh"

      !*** local variables ***
      logical value
      integer k
      integer integrand(2)

      value = BA_push_get(mt_dbl,nrange,
     >                    'integrand',integrand(2),integrand(1))
      if (.not.value) call errquit('indef_integr: push stack',0, MA_ERR)

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

      value = BA_pop_stack(integrand(2))
      if (.not.value) call errquit('indef_integr: pop stack',0, MA_ERR)


      return
      end !SUBROUTINE indef_integr



!*************************************************
!
!   Name    : c_def_integr1
!
!
!   Purpose :
!
!
!   Created :
!
!*************************************************
       complex*16 function c_def_integr1(power_f,
     >                           f,
     >                           power_r,
     >                           r,
     >                           log_amesh,
     >                           nrange,istart)
                    
      implicit none
      integer power_f
      complex*16  f(*)
      integer power_r
      real*8 r(*)
      real*8 log_amesh
      integer nrange,istart

#include "bafdecls.fh"
#include "errquit.fh"


      !*** local variables ***
      integer i
      logical value
      integer k
      complex*16 sum_f
      complex*16 tmp_sum
      integer integrand(2)


      value = BA_push_get(mt_dcpl,nrange,
     >                    'integrand',integrand(2),integrand(1))
      if (.not.value) call errquit('c_def_integr1: push stack',0,MA_ERR)

      do k=1,nrange
        dcpl_mb(integrand(1)+k-1) = f(k)*r(k)**(power_r+1)
      end do

      !*** integrate from the origin to the first point ***
      if (istart.eq.1) then
        sum_f = dcpl_mb(integrand(1))/(power_r+power_f+1)
      else
        sum_f = dcmplx(0.0d0,0.0d0)
      end if

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

      value = BA_pop_stack(integrand(2))
      if (.not.value) call errquit('c_def_integr1: pop stack',1,MA_ERR)

      c_def_integr1 = sum_f

      return
      end !FUNCTION c_def_integr1



!*************************************************
!
!   Name    : def_integr1
!
!
!   Purpose :
!
!
!   Created :
!
!*************************************************
      real*8 function def_integr1(power_f,
     >                           f,
     >                           power_r,
     >                           r,
     >                           log_amesh,
     >                           nrange,istart)
                    
      implicit none
      integer power_f
      real*8  f(*)
      integer power_r
      real*8 r(*)
      real*8 log_amesh
      integer nrange,istart

#include "bafdecls.fh"
#include "errquit.fh"

      !*** local variables ***
      logical value
      integer i,k
      real*8 sum_f,tmp_sum
      integer integrand(2)

      !*** external functions ***
      real*8   dsum
      external dsum

      value = BA_push_get(mt_dbl,nrange,
     >                    'integrand',integrand(2),integrand(1))
      if (.not.value) call errquit('def_integr1: push stack',0,MA_ERR)

      do k=1,nrange
        dbl_mb(integrand(1)+k-1) = f(k)*r(k)**(power_r+1)
      end do

      !*** integrate from the origin to the first point ***
      if (istart.eq.1) then
        sum_f = dbl_mb(integrand(1))/(power_r+power_f+1)
      else
        sum_f = 0.0d0
      end if
      
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

      value = BA_pop_stack(integrand(2))
      if (.not.value) call errquit('def_integr1: pop stack',1,MA_ERR)

      def_integr1 = sum_f

      return
      end !FUNCTION def_integr1

