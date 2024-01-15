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
      integer idum, nnhos, nnhsy
      double precision tau
      double precision nh_freq
      character*32 test
      character*32 thermostat
      character*32 integrator_nucl
      character*32 integrator_elec
      character*32 ext_forces
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
      nh_freq = 0.013669d0
      nnhos = 3
      nnhsy = 3
      if (.not.rtdb_put(rtdb,'qmd:nhfreq',mt_dbl,1,nh_freq))
     & call errquit(pname//'failed to write tau',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'qmd:nnhos',mt_int,1,nnhos))
     & call errquit(pname//'failed to write tau',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'qmd:nnhsy',mt_int,1,nnhsy))
     & call errquit(pname//'failed to write tau',0,RTDB_ERR)
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
  200     if (.not.inp_a(test)) goto 10
  201     if (inp_compare(.false.,'nnhos',test)) then
            if (.not.inp_i(nnhos))
     &      call errquit(pname//'failed to read number of oscillators',
     &         0,INPUT_ERR)
            if (.not.rtdb_put(rtdb,'qmd:nnhos',mt_int,1,nnhos))
     &        call errquit(pname//'failed to write nnhos',0,RTDB_ERR)
            goto 200
          elseif(inp_compare(.false.,'frequency',test)) then
            if (.not.inp_f(nh_freq))
     &      call errquit(pname//'failed to read frequency',0,INPUT_ERR)
            if (.not.rtdb_put(rtdb,'qmd:nhfreq',mt_dbl,1,nh_freq))
     &        call errquit(pname//'failed to write nhfreq',0,RTDB_ERR)
            if (.not.inp_a(test)) goto 10
            if (inp_compare(.false.,'1/cm',test)) then
              nh_freq = nh_freq*4.5563352529120d-6
            elseif (inp_compare(.false.,'cm^-1',test)) then
              nh_freq = nh_freq*4.5563352529120d-6
            elseif (inp_compare(.false.,'ev',test)) then
              nh_freq = nh_freq*3.6749322175655d-2
            elseif (inp_compare(.false.,'au',test)) then
              nh_freq = nh_freq
            elseif (inp_compare(.false.,'ha',test)) then
              nh_freq = nh_freq
            else
              goto 201
            endif
            if (.not.rtdb_put(rtdb,'qmd:nhfreq',mt_dbl,1,nh_freq))
     &        call errquit(pname//'failed to write nnhos',0,RTDB_ERR)
            goto 200
          elseif(inp_compare(.false.,'nsy',test)) then
            if (.not.inp_i(nnhsy))
     &      call errquit(pname//'failed to number of Suzuki-Yoshida'//
     &                        ' points',0,INPUT_ERR)
            if (.not.rtdb_put(rtdb,'qmd:nnhsy',mt_int,1,nnhsy))
     &        call errquit(pname//'failed to write nnhsy',0,RTDB_ERR)
            goto 200
          else
            call errquit(pname//'unknown nose-hoover option '//
     $                   trim(test),0,INPUT_ERR)
          endif

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
       if(.not.rtdb_put(rtdb,'qmd:cavity_axis',mt_dbl,3,cavity_axis(1)))
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
c     Print xyz trajectory every print_xyz steps
      else if (inp_compare(.false.,test,'ext_forces')) then
c     external forces
c     defaults to plumed for now
        if (.not.inp_a(ext_forces)) ext_forces='plumed'
        if (.not.rtdb_cput(rtdb,'qmd:ext_forces',1,
     E       ext_forces(1:inp_strlen(ext_forces))))
     1   call errquit(pname//'failed to write ext forces',0,RTDB_ERR)
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
      subroutine qmd_rdinput(rtdb)
      implicit none
      integer rtdb ! [in]
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer nstep_nucl, nnhos, nnhsy
      double precision dt_nucl
      double precision targ_temp
      double precision nhfreq
      double precision tau
      character*32 thermostat
      character*32 integrator_nucl
      character*13 pname
c
      pname = "qmd_rdinput: "
c     load default values if not set
c     total nuclear steps
      if (.not.rtdb_get(rtdb,'qmd:nstep_nucl',mt_int,1,nstep_nucl)) then
         nstep_nucl = 1000
         if (.not.rtdb_put(rtdb,'qmd:nstep_nucl',mt_int,1,nstep_nucl))
     &      call errquit(pname//'failed to write nstep_nucl',0,RTDB_ERR)
      endif
c     Nuclear MD time step: dt_nucl
      if (.not.rtdb_get(rtdb,'qmd:dt_nucl',mt_dbl,1,dt_nucl)) then
         dt_nucl = 10d0
         if (.not.rtdb_put(rtdb,'qmd:dt_nucl',mt_dbl,1,dt_nucl))
     &        call errquit(pname//'failed to write dt_nucl',0,RTDB_ERR)
      endif
c     Target nuclear MD temperature: targ_temp (Kelvin)
      if (.not.rtdb_get(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))then
         targ_temp = 298.15d0   ! Kelvin
         if (.not.rtdb_put(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))
     &       call errquit(pname//'failed to write targ_temp',0,RTDB_ERR)
      endif
c     Nuclear thermostat
      if (.not.rtdb_cget(rtdb,'qmd:thermostat',1,thermostat))then
         thermostat='none'
         if (.not.rtdb_cput(rtdb,'qmd:thermostat',1,thermostat))
     &      call errquit(pname//'failed to write thermostat',0,RTDB_ERR)
      endif
c     Nose-Hoover mass parameter
      if (.not.rtdb_get(rtdb,'qmd:nhfreq',mt_dbl,1,nhfreq))then
         nhfreq = 0.013669d0
         if (.not.rtdb_put(rtdb,'qmd:nhfreq',mt_dbl,1,nhfreq))
     &        call errquit(pname//'failed to write nhfreq',0,RTDB_ERR)
      endif
c     Nose-Hoover number of oscillators
      if (.not.rtdb_get(rtdb,'qmd:nnhos',mt_int,1,nnhos))then
         nnhos = 3
         if (.not.rtdb_put(rtdb,'qmd:nnhos',mt_int,1,nnhos))
     &        call errquit(pname//'failed to write nnhos',0,RTDB_ERR)
      endif
c     Nose-Hoover number of Suzuki-Yoshida points
      if (.not.rtdb_get(rtdb,'qmd:nnhsy',mt_int,1,nnhsy))then
         nnhsy = 3
         if (.not.rtdb_put(rtdb,'qmd:nnhsy',mt_int,1,nnhsy))
     &        call errquit(pname//'failed to write nnhsy',0,RTDB_ERR)
      endif
c     Berendsen tau parameter, also used for svr
      if (.not.rtdb_get(rtdb,'qmd:tau',mt_dbl,1,tau)) then
         tau=1.d3
         if (.not.rtdb_put(rtdb,'qmd:tau',mt_dbl,1,tau))
     &        call errquit(pname//'failed to write tau',0,RTDB_ERR)
      endif
c     Nuclear MD integrator
      if (.not.rtdb_cget(rtdb,'qmd:integrator_nucl',1,
     I     integrator_nucl)) then
         integrator_nucl='velocity-verlet'
         if (.not.rtdb_cput(rtdb,'qmd:integrator_nucl',1,
     1        integrator_nucl))
     &        call errquit(pname//'failed to write nuclear integrator',
     &        0,RTDB_ERR)
      endif
      return
      end
