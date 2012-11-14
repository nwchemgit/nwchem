c
c     QMD input handler
c
      subroutine qmd_input(rtdb)
c
c     QMD
c     END
c     
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "inp.fh"
#include "stdio.fh"
#include "case.fh"
c
      integer rtdb
c
      double precision dt
      integer nsteps
      double precision dt_elec
      integer nsteps_elec
      double precision targ_temp
      double precision tau
      double precision nh_mass
      character*32 test
      character*32 thermostat
      character*32 integrator
      character*32 integrator_elec
c
      character*30 pname
c
      pname = "qmd_input: "
c
c     ------------------
c     Set default values
c     ------------------
c
c     Nuclear MD time step: dt
      dt = 0.01d0
      if (.not.rtdb_put(rtdb,'qmd:dt',mt_dbl,1,dt))
     & call errquit(pname//'failed to write dt',0,RTDB_ERR)
c
c     Total number of nuclear time steps: nsteps
      nsteps = 1000
      if (.not.rtdb_put(rtdb,'qmd:nsteps',mt_int,1,nsteps))
     & call errquit(pname//'failed to write nsteps',0,RTDB_ERR)
c
c     Electronic dynamics time step: dt_elec
      dt_elec = 0.01d0
      if (.not.rtdb_put(rtdb,'qmd:dt_elec',mt_dbl,1,dt_elec))
     & call errquit(pname//'failed to write dt_elec',0,RTDB_ERR)
c
c     Total number of electronic time steps: nsteps_elec
      nsteps_elec = 100
      if (.not.rtdb_put(rtdb,'qmd:nsteps_elec',mt_int,1,nsteps_elec))
     & call errquit(pname//'failed to write nsteps_elec',0,RTDB_ERR)
c
c     Target nuclear MD temperature: targ_temp (Kelvin)
      targ_temp = 300.d0  ! Kelvin
      if (.not.rtdb_put(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))
     & call errquit(pname//'failed to write targ_temp',0,RTDB_ERR)
c
c     Nuclear thermostat
      thermostat='berendsen'
      if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     & call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
c
c     Berendsen tau parameter
      tau=0.1d0
      if (.not.rtdb_put(rtdb,'qmd:tau',mt_dbl,1,tau))
     & call errquit(pname//'failed to write tau',0,RTDB_ERR)
c
c     Nose-Hoover mass parameter
      nh_mass=0.1d0
      if (.not.rtdb_put(rtdb,'qmd:nh_mass',mt_dbl,1,nh_mass))
     & call errquit(pname//'failed to write nh_mass',0,RTDB_ERR)
c
c     Nuclear MD integrator
      integrator='velocity_verlet'
      if (.not.rtdb_cput(rtdb,'qmd:integrator',1,integrator))
     & call errquit(pname//'failed to write nuclear integrator',0,
     &     RTDB_ERR)
c
c     ----------
c     Read input
c     ----------
c
 10   if (.not. inp_read()) 
     1  call errquit(pname//'failed reading input',0,INPUT_ERR)
      if (.not. inp_a(test)) 
     1  call errquit(pname//'failed reading keyword',0,INPUT_ERR)
c
c     Nuclear MD time step: dt
      if (inp_compare(.false.,test,'dt')) then
       if (.not.inp_f(dt)) then
          write(luout,*) 'dt: value not found; ',
     1      'default value of 0.01d0 will be used'
          dt = 0.01d0
       endif
       if (.not.rtdb_put(rtdb,'qmd:dt',mt_dbl,1,dt))
     &  call errquit(pname//'failed to write dt',0,RTDB_ERR)
c
c     Total number of nuclear time steps: nsteps
      else if (inp_compare(.false.,test,'nsteps')) then
       if (.not.inp_i(nsteps)) then
          write(luOut,*) 'nsteps: value not found; ',
     1      'default value of 1000 will be used'
          nsteps = 1000
       endif
       if (.not.rtdb_put(rtdb,'qmd:nsteps',mt_int,1,nsteps))
     &  call errquit(pname//'failed to write nsteps',0,RTDB_ERR)
c
c     Electronic dynamics time step: dt_elec
      else if (inp_compare(.false.,test,'dt_elec')) then
       if (.not.inp_f(dt_elec)) then
          write(luOut,*) 'dt_elec: value not found; ',
     1      'default value of 0.01d0 will be used'
          dt_elec = 0.01d0
       endif
       if (.not.rtdb_put(rtdb,'qmd:dt_elec',mt_dbl,1,dt_elec))
     &  call errquit(pname//'failed to write dt_elec',0,RTDB_ERR)
c
c     Total number of electronic time steps: nsteps_elec
      else if (inp_compare(.false.,test,'nsteps_elec')) then
       if (.not.inp_i(nsteps_elec)) then
          write(luOut,*) 'nsteps_elec: value not found; ',
     1      'default value of 100 will be used'
          nsteps_elec = 100
       endif
       if (.not.rtdb_put(rtdb,'qmd:nsteps_elec',mt_int,1,nsteps_elec))
     &  call errquit(pname//'failed to write nsteps_elec',0,
     &     RTDB_ERR)
c
c     Target nuclear MD temperature: targ_temp (Kelvin)
      else if (inp_compare(.false.,test,'targ_temp')) then
       if (.not.inp_f(targ_temp)) then
          write(luOut,*) 'targ_temp: value not found; ',
     1      'default value of 300.d0 will be used'
          targ_temp = 300.d0  ! Kelvin
       endif
       if (.not.rtdb_put(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))
     &  call errquit(pname//'failed to write targ_temp',0,
     &     RTDB_ERR)
c
c     Nuclear thermostat
      else if (inp_compare(.false.,test,'thermostat')) then
       if (.not.inp_a(thermostat)) then
          write(LuOut,*) 'thermostat not found; ',
     1      'default berendsen thermostat will be used'
          thermostat='berendsen'
       endif
c
       if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     1  call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
c
c     Berendsen tau parameter 
      else if (inp_compare(.false.,test,'tau')) then
       if (.not.inp_f(tau)) then
          write(luOut,*) 'tau: value not found; ',
     1      'default value of 0.1d0 will be used'
          tau = 0.1d0
       endif
       if (.not.rtdb_put(rtdb,'qmd:tau',mt_dbl,1,tau))
     &  call errquit(pname//'failed to write tau',0,RTDB_ERR)
c
c     Nose-Hoover mass parameter
      else if (inp_compare(.false.,test,'nh_mass')) then
       if (.not.inp_f(nh_mass)) then
          write(luOut,*) 'nh_mass: value not found; ',
     1      'default value of 0.1d0 will be used'
          nh_mass = 0.1d0
       endif
       if (.not.rtdb_put(rtdb,'qmd:nh_mass',mt_dbl,1,nh_mass))
     &  call errquit(pname//'failed to write nh_mass',0,RTDB_ERR)
c
c     Nuclear MD integrator
      else if (inp_compare(.false.,test,'integrator')) then
       if (.not.inp_a(integrator)) then
          write(LuOut,*) 'nuclear integrator not found; ',
     1      'default velocity verlet will be used'
          integrator='velocity_verlet'
       endif
       if (.not.rtdb_cput(rtdb,'qmd:integrator',1,integrator))
     1  call errquit(pname//'failed to write integrator',0,RTDB_ERR)
c
c     Print
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'qmd')
c
c     END
      else if (inp_compare(.false.,test,'end')) then
        goto 20
      else
        call errquit(pname//'unknown directive',0, INPUT_ERR)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   return
      end
