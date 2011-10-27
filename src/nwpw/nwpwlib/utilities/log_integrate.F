*
* $Id$
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

*     ************************************************
*     *                                              *
*     *            log_integrate_indef               *
*     *                                              *
*     ************************************************
      subroutine log_integrate_indef(power_f,f,
     >                        power_r,r,
     >                        log_amesh,
     >                        nrange,sum_f)
      implicit none
      integer  power_f
      real*8   f(*)
      integer  power_r
      real*8  r(*)
      real*8  log_amesh
      integer nrange
      real*8  sum_f(*)

#include "mafdecls.fh"
#include "errquit.fh"

      !*** local variables ***
      integer k
      integer integrand(2)

      if (.not.MA_push_get(mt_dbl,nrange,
     >                    'integrand',integrand(2),integrand(1)))
     >   call errquit('indef_integr: push stack',0, MA_ERR)

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

      if (.not.MA_pop_stack(integrand(2)))
     >   call errquit('log_integrate_indef: pop stack',0,MA_ERR)
      return
      end 



*     ************************************************
*     *                                              *
*     *            log_multipole_energy              *
*     *                                              *
*     ************************************************
      real*8 function log_multipole_energy(l,nrange,g_r,
     >                                     power_q1,q1,
     >                                     power_q2,q2,log_amesh)
      implicit none
      integer l,nrange
      real*8  g_r(*)
      integer power_q2
      real*8  q1(*) 
      integer power_q1
      real*8  q2(*)
      real*8  log_amesh

#include "mafdecls.fh"
#include "errquit.fh"

c     **** local variables ***
      logical ok
      integer power_f,ii
      integer q1_l(2)
      integer q2_l(2)
      integer f_tmp(2) 
      double precision  w,fourpi

*     **** external functions ****
      real*8   log_integrate_def
      external log_integrate_def

      fourpi = 16.0d0*datan(1.0d0)
      ok =        MA_push_get(mt_dbl,nrange,'q1_l',q1_l(2),q1_l(1))
      ok = ok.and.MA_push_get(mt_dbl,nrange,'q2_l',q2_l(2),q2_l(1))
      ok = ok.and.MA_push_get(mt_dbl,nrange,'f_tmp',f_tmp(2),f_tmp(1))
      if (.not.ok) 
     >   call errquit('log_multipole_energy:push stack',0,MA_ERR)


      call log_integrate_indef(power_q1,q1,l,g_r,
     >                         log_amesh,nrange,dbl_mb(q1_l(1)))
      call log_integrate_indef(power_q1,q2,l,g_r,
     >                         log_amesh,nrange,dbl_mb(q2_l(1)))

      do ii=1,nrange
        dbl_mb(f_tmp(1)+ii-1) = (q1(ii)*dbl_mb(q2_l(1)+ii-1) + 
     >                           dbl_mb(q1_l(1)+ii-1)*q2(ii))/
     >                           g_r(ii)**(l+1)
      end do
      power_f = power_q1 +  power_q2

      w = log_integrate_def(power_f,dbl_mb(f_tmp(1)),0,g_r,
     >                      log_amesh,nrange)
      w = w*fourpi/(2.0d0*l + 1.0d0)

      ok = MA_pop_stack(f_tmp(2))
      ok = ok.and.MA_pop_stack(q2_l(2))
      ok = ok.and.MA_pop_stack(q1_l(2))
      if (.not.ok) 
     >   call errquit('log_multipole_energy: push stack',0,MA_ERR)

      log_multipole_energy = w
      return
      end




