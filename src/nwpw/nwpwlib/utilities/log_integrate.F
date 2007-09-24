*
* $Id: log_integrate.F,v 1.1 2007-09-24 17:00:36 bylaska Exp $
*

*     ************************************************
*     *                                              *
*     *            log_integrate_def                 *
*     *                                              *
*     ************************************************
      real*8 function log_integrate_def(power_f,
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
#include "errquit.fh"

      !*** local variables ***
      logical value
      integer k
      real*8 sum_f,tmp_sum
      integer integrand(2)

      !*** external functions ***
c      real*8   dsum
c      external dsum

      value = MA_push_get(mt_dbl,nrange,
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
     >      + log_amesh*tmp_sum
     >      - 0.5d0*log_amesh
     >        *(dbl_mb(integrand(1))
     >         +dbl_mb(integrand(1)+nrange-1))

      value = MA_pop_stack(integrand(2))
      if (.not.value) call errquit('def_integr: pop stack',1, MA_ERR)

      log_integrate_def = sum_f
      return
      end 

