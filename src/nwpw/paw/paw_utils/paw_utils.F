c
c $Id: paw_utils.F,v 1.1 2003-01-31 02:01:13 bylaska Exp $
c

      function find_multipole_energy(l,nrange,g_r,
     >                               power_q1,q1,power_q2,q2,log_amesh) 
      implicit none
      integer           l
      integer           nrange
      double precision  g_r(nrange)
      integer           power_q2
      double precision  q1(nrange)
      integer           power_q1
      double precision  q2(nrange)
      double precision  log_amesh

      double precision find_multipole_energy  !result
    
#include "paw_params.fh"
#include "integrate.fh"

      !*** local variables ***
      integer   power_f,ii
      double precision  q1_l(nrange)
      double precision  q2_l(nrange)
      double precision  f_tmp(nrange) 
      double precision  w

      call indef_integr(power_q1,q1,l,g_r,log_amesh,nrange,q1_l)
      call indef_integr(power_q1,q2,l,g_r,log_amesh,nrange,q2_l)

      do ii=1,nrange
        f_tmp(ii)   = q1(ii)*q2_l(ii) + q1_l(ii)*q2(ii)
        f_tmp(ii)  = f_tmp(ii)/g_r(ii)**(l+1)
      end do
      power_f = power_q1 +  power_q2 

      w = def_integr(power_f,f_tmp,0,g_r,log_amesh,nrange)
      w = w*FOURPI/(2.0d0*l + 1.0d0)

      find_multipole_energy = w
      return
      end 

