c
c     thermostat the velocities
c
      subroutine qmd_thermostat(thermostat,dt,tau,targ_temp,nat,ndof,
     &         mass,v,ekin,final_temp)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
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
c     preliminaries
      kb    = 0.316682942d-05   ! Boltzmann constant in atomic units
      kbinv = 3.157732432d05    ! Inverse Boltzmann constant
c
c     instantaneous temperature
      inst_temp = 2.0d0*ekin/kb/dble(ndof)
c
      if (ga_nodeid().eq.0) 
     &    write(luout,*) "xxx thermostat: ",thermostat
c
      if (thermostat.eq.'berendsen') then
c       calculate velocity scaling factor: Berendsen scaling
        ratio_temp = targ_temp/inst_temp
        velscale = dsqrt(1.d0 + (dt/tau)*(ratio_temp - 1.d0))
        do i = 1,nat
          do j = 1,3
            v(j,i) = velscale*v(j,i)
          end do
        end do
      else
              
      end if
c
c     calculate the kinetic energy
      call qmd_kinetic(nat,mass,v,ekin)
c
c     instantaneous temperature
      final_temp = 2.0d0*ekin/kb/dble(ndof)
c
      return
      end
