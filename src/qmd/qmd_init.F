c
c     qmd_init: initialize parameters
c     
      subroutine qmd_init(rtdb,nstep_nucl,dt_nucl,do_none,fstep)
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
      integer states                ! number of electronic states
      integer init_st               ! initially populated state
      integer fstep                 ! how often to remove trans/rot
      integer tvals(8)
      integer svals
      double precision dt_nucl      ! nuclear md time step
      double precision dt_elec      ! electronic time step
      character*32 thermostat       ! thermostat name
      character*32 integrator_nucl  ! nuclear md integrator
      character*32 integrator_elec  ! electronic dynamics integrator
      integer ichain
      double precision nh_mass(2)   ! nose-hoover mass
      logical do_none
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
     &  (LuOut,'QMD Run Parameters',40,.true.)
       write(LuOut,*)
       write(LuOut,*)
      endif ! ga_nodeid
c
c     boltzmann constant in atomic units
      kb=3.1668114d-6
c
c     convert time in au to femtoseconds
      au2fs=0.02418884326505d0
c
c     convert mass in amu to au
      amu2au=1822.8884845d0
c
c     convert length in au to Angstrom
      au2ang=0.52917721092d0
c
c     total nuclear steps
      if (.not.rtdb_get(rtdb,'qmd:nstep_nucl',mt_int,1,nstep_nucl))
     & call errquit(pname//'failed to read nstep_nucl',0,RTDB_ERR)
      if (ga_nodeid().eq.0) 
     &    write(luout,110) "No. of nuclear steps:",nstep_nucl
c
c     nuclear time step
      if (.not.rtdb_get(rtdb,'qmd:dt_nucl',mt_dbl,1,dt_nucl))
     & call errquit(pname//'failed to read dt_nucl',0,RTDB_ERR)
      if (ga_nodeid().eq.0)
     &    write(luout,120) "Nuclear time step:",dt_nucl
c
c     total electronic steps
      if(.not.rtdb_get(rtdb,'qmd:nstep_elec',mt_int,1,nstep_elec))
     & call errquit(pname//'failed to read nstep_elec',0,RTDB_ERR)
c      if (ga_nodeid().eq.0) 
c     &     write(luout,110) "No. of electronic steps:",nstep_elec
c
c     electronic time step
      if (.not.rtdb_get(rtdb,'qmd:dt_elec',mt_dbl,1,dt_elec))
     & call errquit(pname//'failed to read dt_elec',0,RTDB_ERR)
c      if (ga_nodeid().eq.0)
c     &     write(luout,*) "Electronic time step:",dt_elec
c
c     target temperature
      if (.not.rtdb_get(rtdb,'qmd:targ_temp',mt_dbl,1,targ_temp))
     & call errquit(pname//'failed to read targ_temp',0,RTDB_ERR)
      if (ga_nodeid().eq.0) 
     &        write(luout,120) "Target temp. (K):",targ_temp
c
c     thermostat
      do_berendsen = .false.
      do_nosehoover = .false.
      do_langevin = .false.
      do_BDP = .false.
      do_rescale = .false.
      do_none = .false.
      if (.not. rtdb_cget(rtdb,'qmd:thermostat',1,thermostat))
     $ call errquit(pname//'failed to read thermostat',0,RTDB_ERR)
      if (ga_nodeid().eq.0)
     $    write(luout,140) "Thermostat:",trim(thermostat)
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
       else if (thermostat.eq.'langevin') then
         do_langevin = .true.
       else if (thermostat.eq.'rescale') then
         do_rescale = .true.
       else if (thermostat.eq.'BDP') then
         do_BDP = .true.
       else if (thermostat.eq.'none') then
         do_none = .true.
      else
         if (ga_nodeid().eq.0) 
     &        write(luout,*) "unknown thermostat"
         if (ga_nodeid().eq.0) 
     &        write(luout,*) "using default: NO thermostat"
         do_none = .true.
      end if
c
c     nose-hoover mass parameter
      if (.not.rtdb_get(rtdb,'qmd:nh_mass',mt_dbl,2,nh_mass))
     & call errquit(pname//'failed to read nh_mass',0,RTDB_ERR)
      if (do_nosehoover) then
       if (ga_nodeid().eq.0) 
     &      write(luout,130) "NH mass(1):",nh_mass(1)
       if (ga_nodeid().eq.0) 
     &      write(luout,130) "NH mass(2):",nh_mass(2)
      end if
c
c     langevin friction parameter
      friction = 0.1d0
      if (.not.rtdb_get(rtdb,'qmd:friction',mt_dbl,1,friction))
     & friction = 0.1d0
      if (do_langevin) then
       if (ga_nodeid().eq.0) write(luout,130) "Friction:",friction
      end if
c
c     berendsen tau parameter, also used for BDP
      if (.not.rtdb_get(rtdb,'qmd:tau',mt_dbl,1,tau))
     & call errquit(pname//'failed to read tau',0,RTDB_ERR)
      if (do_berendsen.or.do_BDP) then
       if (ga_nodeid().eq.0) write(luout,120) "Tau:",tau
      end if
c
c     check dt_nucl against the berendsen tau
      if (do_berendsen.and.(dt_nucl.ge.tau)) then
        tau = 1.25d0*dt_nucl
        if (ga_nodeid().eq.0) then
         write(luout,*) "dt_nucl has to be smaller than tau"
         write(luout,*) "adjusting tau to accomodate"
         write(luout,*) "new tau: ", tau
        end if
      end if
c
c     random number seed
      call date_and_time(values=tvals)
c     current seconds x milliseconds + minutes
      svals = tvals(7)*tvals(8)+tvals(6)
      idum = -svals
      if (.not.rtdb_get(rtdb,'qmd:rand_seed',mt_int,1,idum))
     & idum = -svals ! default
      if (ga_nodeid().eq.0) write(luout,110) "Random seed:",idum
c
c     nuclear md integrator
      if (.not. rtdb_cget(rtdb,'qmd:integrator_nucl',1,integrator_nucl))
     $ call errquit(pname//'failed to read integrator',0,RTDB_ERR)
      if (ga_nodeid().eq.0) 
     & write(luout,140) "Nuclear integrator:",trim(integrator_nucl)
      if (integrator_nucl.eq.'velocity-verlet') do_veloverlet = .true.
c
c     is linear molecule ?
      if (.not.rtdb_get(rtdb,'qmd:linear',mt_log,1,do_linear))
     &     do_linear = .false.
c
c     remove translations and rotations every fsteps
      if (.not.rtdb_get(rtdb,'qmd:fstep',mt_int,1,fstep))
     $     fstep = 1000
      if (ga_nodeid().eq.0) call util_flush(LuOut)
c
  110 format(A25,I20)
  120 format(A25,F20.2)
  130 format(A25,F20.6)
  140 format(A25,A20)
c
      return
      end
c $Id$
