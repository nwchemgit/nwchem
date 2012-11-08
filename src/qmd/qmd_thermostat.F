c
c     thermostat the velocities
c
      subroutine qmd_thermostat(thermostat,dt,tau,targ_temp,nat,ndof,
     &         mass,v,ekin,final_temp)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
      character*20 thermostat
      double precision dt, tau, targ_temp
      integer nat
      integer ndof
      double precision mass(nat)
      double precision v(3,nat)
      double precision ekin
      double precision final_temp
c
      integer i,j
      double precision kb,kbinv
      double precision inst_temp,ratio_temp
      double precision velscale
c
      character*30 pname
      pname='qmd_thermostat: '
c
c     preliminaries
      kb    = 0.316682942d-05   ! Boltzmann constant in atomic units
      kbinv = 3.157732432d05    ! Inverse Boltzmann constant
c
c     instantaneous temperature
      inst_temp = 2.0d0*ekin/kb/dble(ndof)
c
      if (thermostat.eq.'berendsen') then
        call qmd_berendsen(dt,tau,inst_temp,targ_temp,nat,v)
      else if (thermostat.eq.'nose-hoover') then
c       call qmd_nosehoover(dt,tau,targ_temp,nat,v)
      else
        call errquit(pname//'unknown thermostat',0,CALC_ERR)
      end if ! choice of thermostat
c
c     calculate the kinetic energy
      call qmd_kinetic(nat,mass,v,ekin)
c
c     instantaneous temperature
      final_temp = 2.0d0*ekin/kb/dble(ndof)
c
      return
      end
