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
      double precision dt_nucl
      integer nstep_nucl
      double precision dt_elec
      integer nstep_elec
      integer fstep
      double precision targ_temp
      double precision friction
      integer idum
      double precision tau
      double precision nh_mass(2)
      character*32 test
      character*32 thermostat
      character*32 integrator_nucl
      character*32 integrator_elec
      logical do_linear
c
      character*30 pname
c
      pname = "qmd_input: "
c
c     ------------------
c     Set default values
c     ------------------
c
c     Nuclear MD time step: dt_nucl
      dt_nucl = 10.d0
      if (.not.rtdb_put(rtdb,'qmd:dt_nucl',mt_dbl,1,dt_nucl))
     & call errquit(pname//'failed to write dt_nucl',0,RTDB_ERR)
c
c     Total number of nuclear time steps: nstep_nucl
      nstep_nucl = 1000
      if (.not.rtdb_put(rtdb,'qmd:nstep_nucl',mt_int,1,nstep_nucl))
     & call errquit(pname//'failed to write nstep_nucl',0,RTDB_ERR)
c
c     Electronic dynamics time step: dt_elec
      dt_elec = 0.01d0
      if (.not.rtdb_put(rtdb,'qmd:dt_elec',mt_dbl,1,dt_elec))
     & call errquit(pname//'failed to write dt_elec',0,RTDB_ERR)
c
c     Total number of electronic time steps: nstep_elec
      nstep_elec = 100
      if (.not.rtdb_put(rtdb,'qmd:nstep_elec',mt_int,1,nstep_elec))
     & call errquit(pname//'failed to write nstep_elec',0,RTDB_ERR)
c
c     Target nuclear MD temperature: targ_temp (Kelvin)
      targ_temp = 298.d0  ! Kelvin
      if (.not.rtdb_put(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))
     & call errquit(pname//'failed to write targ_temp',0,RTDB_ERR)
c
c     Nuclear thermostat
      thermostat='none'
      if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     & call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
c
c     Berendsen tau parameter, also used for BDP
      tau=1.d3
      if (.not.rtdb_put(rtdb,'qmd:tau',mt_dbl,1,tau))
     & call errquit(pname//'failed to write tau',0,RTDB_ERR)
c
c     Nose-Hoover mass parameter
      nh_mass(1)=0.1d0
      nh_mass(2)=0.1d0
      if (.not.rtdb_put(rtdb,'qmd:nh_mass',mt_dbl,2,nh_mass))
     &  call errquit(pname//'failed to write nh_mass',0,RTDB_ERR)
c
c     Nuclear MD integrator
      integrator_nucl='velocity-verlet'
      if (.not.rtdb_cput(rtdb,'qmd:integrator_nucl',1,integrator_nucl))
     & call errquit(pname//'failed to write nuclear integrator',0,
     &     RTDB_ERR)
c
c     Is linear molecule ?
      do_linear=.false.
      if (.not.rtdb_put(rtdb,'qmd:linear',mt_log,1,do_linear))
     & call errquit(pname//'failed to write do_linear',0,RTDB_ERR)
c
c     Remove translations/rotations every fsteps
      fstep=1000
      if (.not.rtdb_put(rtdb,'qmd:fstep',mt_int,1,fstep))
     & call errquit(pname//'failed to write fstep',0,RTDB_ERR)
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
c     Nuclear MD time step: dt_nucl in femtoseconds
      if (inp_compare(.false.,test,'dt_nucl')) then
       if (.not.inp_f(dt_nucl)) then
          write(luout,*) 'dt_nucl: value not found; ',
     1      'default value of 10.0 au will be used'
          dt_nucl = 10.d0
       endif
       if (.not.rtdb_put(rtdb,'qmd:dt_nucl',mt_dbl,1,dt_nucl))
     &  call errquit(pname//'failed to write dt_nucl',0,RTDB_ERR)
c
c     Total number of nuclear time steps: nstep_nucl
      else if (inp_compare(.false.,test,'nstep_nucl')) then
       if (.not.inp_i(nstep_nucl)) then
          write(luOut,*) 'nstep_nucl: value not found; ',
     1      'default value of 1000 will be used'
          nstep_nucl = 1000
       endif
       if (.not.rtdb_put(rtdb,'qmd:nstep_nucl',mt_int,1,nstep_nucl))
     &  call errquit(pname//'failed to write nstep_nucl',0,RTDB_ERR)
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
c     Total number of electronic time steps: nstep_elec
      else if (inp_compare(.false.,test,'nstep_elec')) then
       if (.not.inp_i(nstep_elec)) then
          write(luOut,*) 'nstep_elec: value not found; ',
     1      'default value of 100 will be used'
          nstep_elec = 100
       endif
       if (.not.rtdb_put(rtdb,'qmd:nstep_elec',mt_int,1,nstep_elec))
     &  call errquit(pname//'failed to write nstep_elec',0,
     &     RTDB_ERR)
c
c     Target nuclear MD temperature: targ_temp (Kelvin)
      else if (inp_compare(.false.,test,'targ_temp')) then
       if (.not.inp_f(targ_temp)) then
          write(luOut,*) 'targ_temp: value not found; ',
     1      'default value of 298.d0 will be used'
          targ_temp = 298.d0  ! Kelvin
       endif
       if (.not.rtdb_put(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))
     &  call errquit(pname//'failed to write targ_temp',0,
     &     RTDB_ERR)
c
c     Nuclear thermostat: Nose-Hoover
      else if (inp_compare(.false.,test,'nose-hoover')) then
          thermostat='nose-hoover'
       if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     1  call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
       if (.not.inp_f(nh_mass(1))) nh_mass(1) = 0.1d0
       if (.not.inp_f(nh_mass(2))) nh_mass(2) = 0.1d0
       if (.not.rtdb_put(rtdb,'qmd:nh_mass',mt_dbl,2,nh_mass))
     &  call errquit(pname//'failed to write nh_mass',0,RTDB_ERR)
c
c     Nuclear thermostat: Berendsen
      else if (inp_compare(.false.,test,'berendsen')) then
          thermostat='berendsen'
       if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     1  call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
       if (.not.inp_f(tau)) tau = 1.d3
       if (.not.rtdb_put(rtdb,'qmd:tau',mt_dbl,1,tau))
     &  call errquit(pname//'failed to write tau',0,RTDB_ERR)
c
c     Nuclear thermostat: BDP (Bussi, Donadio, Parrinello 2007)
      else if (inp_compare(.false.,test,'BDP')) then
          thermostat='BDP'
       if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     1  call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
       if (.not.inp_f(tau)) tau = 1.d2
       if (.not.rtdb_put(rtdb,'qmd:tau',mt_dbl,1,tau))
     &  call errquit(pname//'failed to write tau',0,RTDB_ERR)
c
c     Nuclear thermostat: Langevin (Bussi, Parrinello 2007))
      else if (inp_compare(.false.,test,'langevin')) then
          thermostat='langevin'
       if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     1  call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
c
c     Nuclear thermostat: Velocity rescaling
      else if (inp_compare(.false.,test,'rescale')) then
          thermostat='rescale'
       if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     1  call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
c
c     Friction
      else if (inp_compare(.false.,test,'friction')) then
       if (.not.inp_f(friction)) friction = 0.1d0
       if (.not.rtdb_put(rtdb,'qmd:friction',mt_dbl,1,friction))
     &  call errquit(pname//'failed to write friction',0,
     &     RTDB_ERR)
c
c     Random number seed
      else if (inp_compare(.false.,test,'rand_seed')) then
       if (.not.inp_i(idum)) idum = 3989
       idum = sign(idum,-1)   ! seed needs to be negative
       if (.not.rtdb_put(rtdb,'qmd:rand_seed',mt_int,1,idum))
     &  call errquit(pname//'failed to write rand_seed',0,
     &     RTDB_ERR)
c
c     Nuclear MD integrator
      else if (inp_compare(.false.,test,'velocity-verlet')) then
          integrator_nucl='velocity-verlet'
      if (.not.rtdb_cput(rtdb,'qmd:integrator_nucl',1,integrator_nucl))
     1  call errquit(pname//'failed to write integrator',0,RTDB_ERR)
c
c     Is linear molecule ?
      else if (inp_compare(.false.,test,'linear')) then
          do_linear=.true.
      if (.not.rtdb_put(rtdb,'qmd:linear',mt_log,1,do_linear))
     &  call errquit(pname//'failed to write do_linear',0,RTDB_ERR)
c
c     Remove translations and rotations every fsteps
      else if (inp_compare(.false.,test,'fstep')) then
       if (.not.inp_i(fstep)) fstep=1000
       if (.not.rtdb_put(rtdb,'qmd:fstep',mt_int,1,fstep))
     &  call errquit(pname//'failed to write fstep',0,RTDB_ERR)
c
c     Print
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'qmd')
c
c     END
      else if (inp_compare(.false.,test,'end')) then
        goto 20
c
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
c $Id$
