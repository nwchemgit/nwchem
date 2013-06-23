c
c     qmd_init: initialize parameters
c     
      subroutine qmd_init(rtdb,nstep_nucl,dt_nucl)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
c
#include "qmd_common.fh"
c
      integer rtdb                  ! Runtime database handle
      integer ao_bas_han            ! AO basis set handle
      integer nstep_nucl            ! nuclear steps
      integer nstep_elec            ! total electronic steps
      double precision dt_nucl      ! nuclear md time step
      double precision dt_elec      ! electronic time step
      character*32 thermostat       ! thermostat name
      character*32 integrator_nucl  ! nuclear md integrator
      character*32 integrator_elec  ! electronic dynamics integrator
      integer ichain
      double precision nh_mass(2)   ! nose-hoover mass
c
      character*32 pname
c
      pname = "qmd_init: "
c
c     print header and general information
      if (ga_nodeid().eq.0) then
       write(LuOut,*)
       write(LuOut,*)
       call util_print_centered(LuOut,'NWChem QMD Module',40,.true.)
       write(LuOut,*)
       write(LuOut,*)
c
       write(LuOut,*)
       write(LuOut,*)
       call util_print_centered
     &  (LuOut,'QMD General Information',40,.true.)
       write(LuOut,*)
       write(LuOut,*)
      endif ! ga_nodeid
c
c     boltzmann constant in atomic units
      kb=3.16679d-6
c
c     convert time in au to femtoseconds
      au2fs=0.02419d0
c
c     total nuclear steps
      if (.not.rtdb_get(rtdb,'qmd:nstep_nucl',mt_int,1,nstep_nucl))
     & call errquit(pname//'failed to read nstep_nucl',0,RTDB_ERR)
      if (ga_nodeid().eq.0) write(luout,*) "nstep_nucl: ",nstep_nucl
c
c     nuclear time step
      if (.not.rtdb_get(rtdb,'qmd:dt_nucl',mt_dbl,1,dt_nucl))
     & call errquit(pname//'failed to read dt_nucl',0,RTDB_ERR)
      if (ga_nodeid().eq.0) write(luout,*) "dt_nucl: ",dt_nucl
c
c     total electronic steps
      if(.not.rtdb_get(rtdb,'qmd:nstep_elec',mt_int,1,nstep_elec))
     & call errquit(pname//'failed to read nstep_elec',0,RTDB_ERR)
      if (ga_nodeid().eq.0) write(luout,*) "nstep_elec: ",nstep_elec
c
c     electronic time step
      if (.not.rtdb_get(rtdb,'qmd:dt_elec',mt_dbl,1,dt_elec))
     & call errquit(pname//'failed to read dt_elec',0,RTDB_ERR)
      if (ga_nodeid().eq.0) write(luout,*) "dt_elec: ",dt_elec
c
c     target temperature
      if (.not.rtdb_get(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))
     & call errquit(pname//'failed to read targ_temp',0,RTDB_ERR)
      if (ga_nodeid().eq.0) 
     &        write(luout,*) "target temperature: ",targ_temp
c
c     nose-hoover mass parameter
      if (.not.rtdb_get(rtdb,'qmd:nh_mass',mt_dbl,2,nh_mass))
     & call errquit(pname//'failed to read nh_mass',0,RTDB_ERR)
      if (ga_nodeid().eq.0) 
     &     write(luout,*) "nh_mass(1): ",nh_mass(1)
      if (ga_nodeid().eq.0) 
     &     write(luout,*) "nh_mass(2): ",nh_mass(2)
c
c     thermostat
      do_berendsen = .false.
      do_nosehoover = .false.
      if (.not. rtdb_cget(rtdb,'qmd:thermostat',1,thermostat))
     $ call errquit(pname//'failed to read thermostat',0,RTDB_ERR)
      if (ga_nodeid().eq.0) write(luout,*) "thermostat: ",thermostat
      if (thermostat.eq.'berendsen') then
         do_berendsen = .true.
      else if (thermostat.eq.'nose-hoover') then
         do_nosehoover = .true.
         do ichain = 1,maxchain
            r_nh(ichain) = 0.d0            ! position of chains
            v_nh(ichain) = 0.d0            ! velocity of chains
            m_nh(ichain) = nh_mass(ichain) ! mass of chains
            g_nh(ichain) = 0.d0            ! bath coupling
         end do
      else
         if (ga_nodeid().eq.0) 
     &        write(luout,*) "unknown thermostat"
         if (ga_nodeid().eq.0) 
     &        write(luout,*) "using default: berendsen thermostat"
         do_berendsen = .true.
      end if
c
c     nuclear md integrator
      if (.not. rtdb_cget(rtdb,'qmd:integrator_nucl',1,integrator_nucl))
     $ call errquit(pname//'failed to read integrator',0,RTDB_ERR)
      if (ga_nodeid().eq.0) 
     & write(luout,*) "integrator_nucl: ",integrator_nucl
      if (integrator_nucl.eq.'velocity-verlet') do_veloverlet = .true.
c
c     berendsen tau parameter
      if (.not.rtdb_get(rtdb,'qmd:tau',mt_dbl,1,tau))
     & call errquit(pname//'failed to read tau',0,RTDB_ERR)
      if (ga_nodeid().eq.0) write(luout,*) "tau: ",tau
c
c     check dt_nucl against the berendsen tau
      if (do_berendsen.and.(dt_nucl.ge.tau)) then
        dt_nucl = 0.75d0*tau 
        if (ga_nodeid().eq.0) then
         write(luout,*) "dt_nucl has to be lesser than tau"
         write(luout,*) "adjusting dt_nucl to accomodate"
         write(luout,*) "new dt_nucl: ", dt_nucl
        end if
      end if
c
c     is linear molecule ?
      if (.not.rtdb_get(rtdb,'qmd:linear',mt_log,1,do_linear))
     &     do_linear = .false.
c
      if (ga_nodeid().eq.0) call util_flush(LuOut)
c
      return
      end
c $Id$
