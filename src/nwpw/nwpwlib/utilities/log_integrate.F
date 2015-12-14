*
* $Id$
*

*     ************************************************
*     *                                              *
*     *            log_integrate_def                 *
*     *                                              *
*     ************************************************
*
*    Computes the following integral
*
*
*         /                        /
*         | f(r) r**power_r dr =   | f(i) * r(i)**(power_r+1) * log(a) di
*         /                        /
*
*    where r(i) = r0*a**i   and f(r-->0) = r**power_f

      real*8 function log_integrate_def_destroyf(power_f,
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


      !*** local variables ***
      logical value
      integer k
      real*8 sum_f,tmp_sum
      integer integrand(2)

      !*** external functions ***
c      real*8   dsum
c      external dsum


      do k=1,nrange
        f(k) = f(k)*r(k)**(power_r+1)
      end do

      !*** integrate from the origin to the first point ***
      sum_f = f(1)/(power_r+power_f+1)
      
      !*** the rest via trapesoidal rule ***
      tmp_sum = 0.0d0
      do i=1,nrange
       tmp_sum = tmp_sum + f(i)
      end do
 
      !*** the rest via trapesoidal rule ***
      sum_f = sum_f 
     >      + log_amesh*tmp_sum
     >      - 0.5d0*log_amesh
     >        *(f(1) + f(nrange))


      log_integrate_def_destroyf = sum_f
      return
      end 


*     ************************************************
*
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

#include "bafdecls.fh"
#include "errquit.fh"

      !*** local variables ***
      logical value
      integer k
      real*8 sum_f,tmp_sum
      integer integrand(2)

      !*** external functions ***
c      real*8   dsum
c      external dsum

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
     >      + log_amesh*tmp_sum
     >      - 0.5d0*log_amesh
     >        *(dbl_mb(integrand(1))
     >         +dbl_mb(integrand(1)+nrange-1))

      value = BA_pop_stack(integrand(2))
      if (.not.value) call errquit('def_integr: pop stack',1, MA_ERR)

      log_integrate_def = sum_f
      return
      end 


*     ************************************************
*     *                                              *
*     *            log_integrate_def0                *
*     *                                              *
*     ************************************************
*
*    Computes the following integral
*
*
*         /               /
*         | f(r)/r dr =   | f(i) * log(a) di
*         /               /
*
*    where r(i) = r0*a**i   and f(r-->0) = r**power_f
*
      real*8 function log_integrate_def0(power_f,
     >                           f,
     >                           r,
     >                           log_amesh,
     >                           nrange)
      implicit none
      integer power_f
      real*8  f(*)
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
c      real*8   dsum
c      external dsum

      value = BA_push_get(mt_dbl,nrange,
     >                    'integrand',integrand(2),integrand(1))
      if (.not.value) call errquit('def_integr: push stack',0, MA_ERR)

      do k=1,nrange
        dbl_mb(integrand(1)+k-1) = f(k)
      end do

      !*** integrate from the origin to the first point ***
      sum_f = dbl_mb(integrand(1))/(power_f)

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

      value = BA_pop_stack(integrand(2))
      if (.not.value) call errquit('def_integr: pop stack',1, MA_ERR)

      log_integrate_def0 = sum_f
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

#include "bafdecls.fh"
#include "errquit.fh"

      !*** local variables ***
      integer k
      integer integrand(2)

      if (.not.BA_push_get(mt_dbl,nrange,
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

      if (.not.BA_pop_stack(integrand(2)))
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

#include "bafdecls.fh"
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
      ok =        BA_push_get(mt_dbl,nrange,'q1_l',q1_l(2),q1_l(1))
      ok = ok.and.BA_push_get(mt_dbl,nrange,'q2_l',q2_l(2),q2_l(1))
      ok = ok.and.BA_push_get(mt_dbl,nrange,'f_tmp',f_tmp(2),f_tmp(1))
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
      w = w*fourpi/dble(2.0d0*l + 1.0d0)

      ok = BA_pop_stack(f_tmp(2))
      ok = ok.and.BA_pop_stack(q2_l(2))
      ok = ok.and.BA_pop_stack(q1_l(2))
      if (.not.ok) 
     >   call errquit('log_multipole_energy: push stack',0,MA_ERR)

      log_multipole_energy = w
      return
      end


*     ************************************************
*     *                                              *
*     *            log_coulomb_energy                *
*     *                                              *
*     ************************************************
      real*8 function log_coulomb_energy(rho,charge,
     >                                   rgrid,ngrid,log_amesh)
      implicit none
      integer ngrid
      real*8  rho(*),rgrid(*),log_amesh,charge

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i
      integer tmp(2),vh(2)
      real*8  E,fourpi,tt

*     **** external functions ****
      real*8   log_corrector_iF,log_r2integrate_eric
      external log_corrector_iF,log_r2integrate_eric

      fourpi = 16.0d0*datan(1.0d0)

      value =           BA_push_get(mt_dbl,ngrid,'tmp',tmp(2),tmp(1))
      value = value.and.BA_push_get(mt_dbl,ngrid,'vh',vh(2),vh(1))
      if (.not.value)
     >   call errquit('log_coulomb_energy:push stack',0,MA_ERR)

      do i=1,ngrid
         dbl_mb(tmp(1)+i-1) = fourpi*rho(i)*log_amesh*rgrid(i)**3
      end do

      dbl_mb(vh(1)+ngrid-1) = charge
      dbl_mb(vh(1)+ngrid-2) = charge
      dbl_mb(vh(1)+ngrid-3) = charge
      do i=(ngrid-3), 1, -1
         dbl_mb(vh(1)+i-1) = dbl_mb(vh(1)+i) 
     >                     + log_corrector_iF(i,dbl_mb(tmp(1)))
      end do

      do i=1,ngrid
         dbl_mb(tmp(1)+i-1)=fourpi*rho(i)*log_amesh*rgrid(i)**2
      end do

      tt = 0.0d0
      do i=(ngrid-3), 1, -1
         tt = tt + log_corrector_iF(i,dbl_mb(tmp(1)))
         dbl_mb(vh(1)+i-1) = dbl_mb(vh(1)+i-1) - rgrid(i)*tt
      end do

      do i=1,ngrid
         dbl_mb(vh(1)+i-1) = dbl_mb(vh(1)+i-1)/rgrid(i)
      end do

      do i=1,ngrid
         dbl_mb(tmp(1)+i-1) = fourpi*rho(i)*dbl_mb(vh(1)+i-1)
      end do
      E  = log_r2integrate_eric(ngrid,log_amesh,rgrid,dbl_mb(tmp(1)))

      value =           BA_pop_stack(vh(2))
      value = value.and.BA_pop_stack(tmp(2))
      if (.not.value)
     >   call errquit('log_coulomb_energy:pop stack',1,MA_ERR)

      log_coulomb_energy = E*0.5d0
      return
      end 

*     ************************************************
*     *                                              *
*     *            log_coulomb0_energy               *
*     *                                              *
*     ************************************************
      real*8 function log_coulomb0_energy(rho,charge,
     >                                    rgrid,ngrid,log_amesh,zion)
      implicit none
      integer ngrid
      real*8  rho(*),rgrid(*),log_amesh,charge,zion

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer i,tmp(2)
      real*8  E,fourpi

*     **** external functions ****
      real*8   log_r2integrate_eric
      external log_r2integrate_eric

      fourpi = 16.0d0*datan(1.0d0)
      if (.not.BA_push_get(mt_dbl,ngrid,'tmp',tmp(2),tmp(1)))
     >   call errquit('log_coulomb0_energy:push stack',0,MA_ERR)

      do i=1,ngrid
         dbl_mb(tmp(1)+i-1) = fourpi*rho(i)/rgrid(i)
      end do
      E  = log_r2integrate_eric(ngrid,log_amesh,rgrid,dbl_mb(tmp(1)))
      E  = -E*zion

      if (.not.BA_pop_stack(tmp(2)))
     >   call errquit('log_coulomb0_energy:pop stack',1,MA_ERR)

      log_coulomb0_energy = E
      return
      end


*     ************************************************
*     *                                              *
*     *            log_corrector_iF                  *
*     *                                              *
*     ************************************************
*
*    Computes y(i) <-- y(i+1) + F(f(i),f(i+1),f(i+2),f(i+3))
*   where F is a 5 point corrector
*
      real*8 function log_corrector_iF(i,f)
      implicit none
      integer i
      real*8 f(*)

*     **** local variables ****
      real*8 oneo24
      parameter (oneo24=1.0d0/24.0d0)
      log_corrector_iF = -oneo24 * (   9.0d0*f(i)
     >                              + 19.0d0*f(i+1)
     >                              -  5.0d0*f(i+2)
     >                              +  1.0d0*f(i+3))
      return
      end

*     ************************************************
*     *                                              *
*     *            log_r2integrate_eric              *
*     *                                              *
*     ************************************************
      real*8 function log_r2integrate_eric(ngrid,log_amesh,rgrid,f)
      implicit none
      integer ngrid
      real*8 log_amesh,rgrid(*),f(*)

*     **** local variables ****
      integer i
      real*8 mysum

      mysum = (   9.0d0*f(1)*(rgrid(1)**3)
     >         + 28.0d0*f(2)*(rgrid(2)**3)
     >         + 23.0d0*f(3)*(rgrid(3)**3)
     >        )/24.0d0
      do i=4,ngrid
         mysum = mysum + f(i)*(rgrid(i)**3)
      end do
      log_r2integrate_eric = log_amesh*mysum + f(1)*(rgrid(1)**3)/3.0d0
      return
      end
