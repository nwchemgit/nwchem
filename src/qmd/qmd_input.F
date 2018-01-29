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
      integer nstep_elec
      integer com_step
      integer print_xyz
      integer step_property
      integer step_tddft
      integer step_freq
      integer tvals(8)
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
      logical do_cavity
      logical do_namd
      logical do_property
      logical do_tddft
      logical do_freq
      double precision cavity_axis(3)
      double precision cavity_spring
     
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
c     Target nuclear MD temperature: targ_temp (Kelvin)
      targ_temp = 298.15d0  ! Kelvin
      if (.not.rtdb_put(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))
     & call errquit(pname//'failed to write targ_temp',0,RTDB_ERR)
c
c     Nuclear thermostat
      thermostat='none'
      if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     & call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
c
c     Berendsen tau parameter, also used for svr
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
c     Remove translations/rotations every com_step steps
      com_step=100
      if (.not.rtdb_put(rtdb,'qmd:com_step',mt_int,1,com_step))
     & call errquit(pname//'failed to write com_step',0,RTDB_ERR)
c
c     Print xyz trajectory every print_xyz steps
      print_xyz=1
      if (.not.rtdb_put(rtdb,'qmd:print_xyz',mt_int,1,print_xyz))
     & call errquit(pname//'failed to write print_xyz',0,RTDB_ERR)
c
c     Do we need cavity constraints ?
      do_cavity=.false.
      if (.not.rtdb_put(rtdb,'qmd:cavity',mt_log,1,do_cavity))
     & call errquit(pname//'failed to write do_cavity',0,RTDB_ERR)
c
c     Are we doing NAMD?
      do_namd=.false.
      if (.not.rtdb_put(rtdb,'qmd:namd',mt_log,1,do_namd))
     & call errquit(pname//'failed to write do_namd',0,RTDB_ERR)
c
c     Are we calculating any properties ?
      do_property=.false.
      if (.not.rtdb_put(rtdb,'qmd:property',mt_log,1,do_property))
     & call errquit(pname//'failed to write do_property',0,RTDB_ERR)
c
c     Calculate molecular properties every step_property steps
      step_property=1
      if (.not.rtdb_put(rtdb,'qmd:step_property',mt_int,
     &  1,step_property))
     & call errquit(pname//'failed to write step_property',0,RTDB_ERR)
c
c     Are we performing any tddft calculations along the ground state trajectory ?
      do_tddft=.false.
      if (.not.rtdb_put(rtdb,'qmd:tddft',mt_log,1,do_tddft))
     & call errquit(pname//'failed to write do_tddft',0,RTDB_ERR)
c
c     Calculate tddft every step_tddft steps
      step_tddft=1
      if (.not.rtdb_put(rtdb,'qmd:step_tddft',mt_int,
     &  1,step_tddft))
     & call errquit(pname//'failed to write step_tddft',0,RTDB_ERR)
c
c     Are we performing any frequency calculations along the ground state trajectory ?
      do_freq=.false.
      if (.not.rtdb_put(rtdb,'qmd:freq',mt_log,1,do_freq))
     & call errquit(pname//'failed to write do_freq',0,RTDB_ERR)
c
c     Calculate frequencies every step_freq steps
      step_freq=1
      if (.not.rtdb_put(rtdb,'qmd:step_freq',mt_int,
     &  1,step_freq))
     & call errquit(pname//'failed to write step_freq',0,RTDB_ERR)
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
c     Target nuclear MD temperature: targ_temp (Kelvin)
      else if (inp_compare(.false.,test,'targ_temp')) then
       if (.not.inp_f(targ_temp)) then
          write(luOut,*) 'targ_temp: value not found; ',
     1      'default value of 298.15 will be used'
          targ_temp = 298.15d0  ! Kelvin
       endif
       if (.not.rtdb_put(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))
     &  call errquit(pname//'failed to write targ_temp',0,
     &     RTDB_ERR)
c
c     Nuclear thermostat
      else if (inp_compare(.false.,test,'thermostat')) then
        if (.not.inp_a(thermostat)) thermostat='none'
        if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     1   call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
c       Nose-Hoover: get masses
        if (thermostat.eq.'nose-hoover') then
          if (.not.inp_f(nh_mass(1))) nh_mass(1) = 0.1d0
          if (.not.inp_f(nh_mass(2))) nh_mass(2) = 0.1d0
          if (.not.rtdb_put(rtdb,'qmd:nh_mass',mt_dbl,2,nh_mass))
     &     call errquit(pname//'failed to write nh_mass',0,RTDB_ERR)
c       Berendsen or svr (Bussi, Donadio, Parrinello 2007): get tau
        else if (thermostat.eq.'berendsen'.or.thermostat.eq.'svr') then
          if (.not.inp_f(tau)) tau = 1.d3
          if (.not.rtdb_put(rtdb,'qmd:tau',mt_dbl,1,tau))
     &     call errquit(pname//'failed to write tau',0,RTDB_ERR)
c       Langevin (Bussi, Parrinello 2007): get friction
        else if (thermostat.eq.'langevin') then
          if (.not.inp_f(friction)) friction = 0.1d0
          if (.not.rtdb_put(rtdb,'qmd:friction',mt_dbl,1,friction))
     &     call errquit(pname//'failed to write friction',0,RTDB_ERR)
c       Rescale doesn't have any parameters
        end if
c
c     Random number seed
      else if (inp_compare(.false.,test,'rand_seed')) then
       if (.not.inp_i(idum)) then
        call errquit(pname//'failed to parse rand_seed',0,
     &      INPUT_ERR)
       else
        if (.not.rtdb_put(rtdb,'qmd:rand_seed',mt_int,1,idum))
     &   call errquit(pname//'failed to write rand_seed',0,
     &      RTDB_ERR)
       end if
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
c     Do we need cavity constraints  ?
      else if (inp_compare(.false.,test,'cavity')) then
        do_cavity=.true.
       if (.not.rtdb_put(rtdb,'qmd:cavity',mt_log,1,do_cavity))
     &  call errquit(pname//'failed to write do_cavity',0,RTDB_ERR)
       if (.not.inp_f(cavity_axis(1))) cavity_axis(1) = 10.d0
       if (.not.inp_f(cavity_axis(2))) cavity_axis(2) = 10.d0
       if (.not.inp_f(cavity_axis(3))) cavity_axis(3) = 10.d0
       if (.not.rtdb_put(rtdb,'qmd:cavity_axis',mt_dbl,3,cavity_axis))
     &  call errquit(pname//'failed to write cavity_axis',0,RTDB_ERR)
c
c     Cavity force
      else if (inp_compare(.false.,test,'cavity_spring')) then
       if (.not.inp_f(cavity_spring)) cavity_spring = 0.d0
       if (.not.rtdb_put(rtdb,'qmd:cavity_spring',mt_dbl,
     &    1,cavity_spring))
     &  call errquit(pname//'failed to write cavity_spring',0,RTDB_ERR)
c
c     Remove translations and rotations every com_step steps
      else if (inp_compare(.false.,test,'com_step')) then
       if (.not.inp_i(com_step)) com_step=100
       if (.not.rtdb_put(rtdb,'qmd:com_step',mt_int,1,com_step))
     &  call errquit(pname//'failed to write com_step',0,RTDB_ERR)
c
c     Print xyz trajectory every print_xyz steps
      else if (inp_compare(.false.,test,'print_xyz')) then
       if (.not.inp_i(print_xyz)) print_xyz=1
       if (.not.rtdb_put(rtdb,'qmd:print_xyz',mt_int,1,print_xyz))
     &  call errquit(pname//'failed to write print_xyz',0,RTDB_ERR)
c
c     NAMD input
      else if (inp_compare(.false.,test,'namd')) then
       do_namd=.true.
       if (.not.rtdb_put(rtdb,'qmd:namd',mt_log,1,do_namd))
     &  call errquit(pname//'failed to write do_namd',0,RTDB_ERR)
       call qmd_namd_input(rtdb)
c
c     Calculate properties as specified in the property block every step_property steps
      else if (inp_compare(.false.,test,'property')) then
       do_property=.true.
       if (.not.rtdb_put(rtdb,'qmd:property',mt_log,1,do_property))
     &  call errquit(pname//'failed to write do_property',0,RTDB_ERR)
       if (.not.inp_i(step_property)) step_property=1
       if (.not.rtdb_put(rtdb,'qmd:step_property',mt_int,1,
     & step_property))
     &  call errquit(pname//'failed to write step_property',0,RTDB_ERR)
c
c     Calculate tddft spectra as specified in the tddft block every step_tddft steps
      else if (inp_compare(.false.,test,'tddft')) then
       do_tddft=.true.
       if (.not.rtdb_put(rtdb,'qmd:tddft',mt_log,1,do_tddft))
     &  call errquit(pname//'failed to write do_tddft',0,RTDB_ERR)
       if (.not.inp_i(step_tddft)) step_tddft=1
       if (.not.rtdb_put(rtdb,'qmd:step_tddft',mt_int,1,
     & step_tddft))
     &  call errquit(pname//'failed to write step_tddft',0,RTDB_ERR)
c
c     Calculate freq spectra as specified in the hessian block every step_freq steps
      else if (inp_compare(.false.,test,'freq')) then
       do_freq=.true.
       if (.not.rtdb_put(rtdb,'qmd:freq',mt_log,1,do_freq))
     &  call errquit(pname//'failed to write do_freq',0,RTDB_ERR)
       if (.not.inp_i(step_freq)) step_freq=1
       if (.not.rtdb_put(rtdb,'qmd:step_freq',mt_int,1,
     & step_freq))
     &  call errquit(pname//'failed to write step_freq',0,RTDB_ERR)
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
