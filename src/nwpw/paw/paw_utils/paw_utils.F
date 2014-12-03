c
c $Id$
c

      function find_multipole_energy(l,nrange,g_r,
     >                               power_q1,q1,power_q2,q2,log_amesh) 
      implicit none
      integer           l
      integer           nrange
      double precision  g_r(*)
      integer           power_q2
      double precision  q1(*)
      integer           power_q1
      double precision  q2(*)
      double precision  log_amesh

      double precision find_multipole_energy  !result
    
#include "paw_params.fh"
#include "integrate.fh"
#include "bafdecls.fh"

      !*** local variables ***
      logical ok
      integer   power_f,ii
      integer   q1_l(2)
      integer q2_l(2)
      integer f_tmp(2) 
c      double precision  q1_l(nrange)
c      double precision  q2_l(nrange)
c      double precision  f_tmp(nrange) 
      double precision  w

      ok = BA_push_get(mt_dbl,nrange,
     >                    'q1_l',q1_l(2),q1_l(1))
      if (.not.ok) call errquit('find_multipole_energy: push stack',0,0)
      ok = BA_push_get(mt_dbl,nrange,
     >                    'q2_l',q2_l(2),q2_l(1))
      if (.not.ok) call errquit('find_multipole_energy: push stack',0,0)
      ok = BA_push_get(mt_dbl,nrange,
     >                    'f_tmp',f_tmp(2),f_tmp(1))
      if (.not.ok) call errquit('find_multipole_energy: push stack',0,0)


      call indef_integr(power_q1,q1,l,g_r,log_amesh,
     >         nrange,dbl_mb(q1_l(1)))
      call indef_integr(power_q1,q2,l,g_r,log_amesh,
     >         nrange,dbl_mb(q2_l(1)))

      do ii=1,nrange
        dbl_mb(f_tmp(1)+ii-1) = (q1(ii)*dbl_mb(q2_l(1)+ii-1) + 
     >                           dbl_mb(q1_l(1)+ii-1)*q2(ii))/
     >                           g_r(ii)**(l+1)
      end do
      power_f = power_q1 +  power_q2 

      w = def_integr(power_f,dbl_mb(f_tmp(1)),0,g_r,log_amesh,nrange)
      w = w*FOURPI/(2.0d0*l + 1.0d0)

      find_multipole_energy = w

      ok = BA_pop_stack(f_tmp(2))
      if (.not.ok) call errquit('find_multipole_energy: push stack',0,0)
      ok = BA_pop_stack(q2_l(2))
      if (.not.ok) call errquit('find_multipole_energy: push stack',0,0)
      ok = BA_pop_stack(q1_l(2))
      if (.not.ok) call errquit('find_multipole_energy: push stack',0,0)

      return
      end 

