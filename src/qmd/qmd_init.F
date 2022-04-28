c
c     qmd_init: initialize parameters
c     
      subroutine qmd_init(rtdb,nstep_nucl,dt_nucl,do_none,com_step,
     $                    print_xyz,restart,do_namd,step_property,
     $                    do_property,step_tddft,do_tddft,step_freq, 
     $                    do_freq)
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
#include "msgtypesf.h"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer rtdb                  ! Runtime database handle
      integer ao_bas_han            ! AO basis set handle
      integer nstep_nucl            ! nuclear steps
      integer nstep_elec            ! total electronic steps
      integer com_step              ! how often to remove trans/rot
      integer print_xyz             ! how often to print xyz trajectory
      integer step_property
      integer step_tddft
      integer step_freq
      character*8 sdate
      character*10 stime
      character*5 szone
      integer tvals(8)
      integer svals
      double precision dt_nucl      ! nuclear md time step
      character*32 thermostat       ! thermostat name
      character*32 integrator_nucl  ! nuclear md integrator
      character*32 integrator_elec  ! electronic dynamics integrator
      integer ichain
      double precision nhfreq   ! nose-hoover mass
      logical do_none
      logical restart
      logical do_namd
      logical do_property
      logical do_tddft
      logical do_freq
      logical qmmm_init
      integer inntsize
c
      character*32 pname
c
      pname = "qmd_init: "
c
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
c
c     print header and general information
      if (ga_nodeid().eq.0) then
       write(LuOut,*)
       write(LuOut,*)
       call util_print_centered(LuOut,'NWChem QMD Module',40,.true.)
       write(LuOut,*)
       write(LuOut,*)
       write(LuOut,*)
       write(LuOut,*)
       call util_print_centered
     &  (LuOut,'QMD Run Parameters',40,.true.)
       write(LuOut,*)
       write(LuOut,*)
      endif ! ga_nodeid
c     check if qmmm is on
       qmmm_init=.false.
      if (.not. rtdb_get(rtdb,'qmmm:init',mt_log,1,qmmm_init))
     $      qmmm_init = .false.
      if(qmmm_init)
     &   call errquit(pname//'qmd module not compatible with qmmm',0,
     &                CAPMIS_ERR)
c
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
      do_svr = .false.
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
       else if (thermostat.eq.'langevin') then
         do_langevin = .true.
       else if (thermostat.eq.'rescale') then
         do_rescale = .true.
       else if (thermostat.eq.'svr') then
         do_svr = .true.
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
      if (.not.rtdb_get(rtdb,'qmd:nhfreq',mt_dbl,1,nhfreq))
     & call errquit(pname//'failed to read nh_mass',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'qmd:nnhos',mt_int,1,nnhos))
     & call errquit(pname//'failed to read nh_nchain',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'qmd:nnhsy',mt_int,1,nnhsy))
     & call errquit(pname//'failed to read nh_nsy',0,RTDB_ERR)
      if (do_nosehoover) then
        do ichain = 1,nnhos
           r_nh(ichain) = 0.d0                                          ! position of chains
           m_nh(ichain) = kb_au*targ_temp/nhfreq**2                     ! mass of chains
           v_nh(ichain) = 0.d0                     ! velocity of chains
           g_nh(ichain) = 0.d0                                          ! bath coupling
        end do

       if (nnhsy.eq.1) then
         syw_nh(1) = 1d0
       elseif (nnhsy.eq.3) then
         syw_nh(1) = 1d0/(2d0 - 2d0**(1d0/3d0))
         syw_nh(2) = 1 - 2*syw_nh(1)
         syw_nh(3) = syw_nh(1)
       elseif (nnhsy.eq.7) then
         syw_nh(1) =  0.784513610477560d0
         syw_nh(2) =  0.235573213359357d0
         syw_nh(3) = -1.17767998417887d0
         syw_nh(4) = 1d0 - 2*sum(syw_nh(1:3))
         syw_nh(5) = syw_nh(3)
         syw_nh(6) = syw_nh(2)
         syw_nh(7) = syw_nh(1)
       endif

       if (ga_nodeid().eq.0) 
     &      write(luout,130) "NH frequency:",nhfreq
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
c     berendsen tau parameter, also used for svr
      if (.not.rtdb_get(rtdb,'qmd:tau',mt_dbl,1,tau))
     & call errquit(pname//'failed to read tau',0,RTDB_ERR)
      if (do_berendsen.or.do_svr) then
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
      svals=0
      if(ga_nodeid().eq.0) then
         call date_and_time(sdate,stime,szone,tvals)
c     current seconds x milliseconds + minutes
         svals = tvals(7)*tvals(8)+tvals(6)
      endif
      call ga_brdcst(Msg_Vec_Stat+MSGINT,svals,inntsize,0)
      call ga_sync()
      if (.not.rtdb_get(rtdb,'qmd:rand_seed',mt_int,1,idum))
     & idum = svals ! default
      idum = sign(idum,-1)   ! seed needs to be negative
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
c     remove translations and rotations every com_step steps
      if (.not.rtdb_get(rtdb,'qmd:com_step',mt_int,1,com_step))
     $     com_step = 100
c
c     print xyz trajectory every print_xyz steps
      if (.not.rtdb_get(rtdb,'qmd:print_xyz',mt_int,1,print_xyz))
     $     print_xyz = 1
c
c     are we restarting calculation?
      if (.not.rtdb_get(rtdb,'qmd:restart',mt_log,1,restart))
     $     restart=.false.
c
c     are we doing NAMD?
      if (.not.rtdb_get(rtdb,'qmd:namd',mt_log,1,do_namd))
     $     do_namd=.false.
c
c     are we calculating molecular properties ?
      if (.not.rtdb_get(rtdb,'qmd:property',mt_log,1,do_property))
     $     do_property=.false.
c
c     calculate molecular properties as defined in the properties block every step_property steps
      if (.not.rtdb_get(rtdb,'qmd:step_property',mt_int,1,
     $ step_property)) step_property = 1
c
c     are we performing any tddft calculations along the ground state trajectory ?
      if (.not.rtdb_get(rtdb,'qmd:tddft',mt_log,1,do_tddft))
     $     do_tddft=.false.
c
c     calculate tddft spectra every step_tddft along the ground state trajectory
      if (.not.rtdb_get(rtdb,'qmd:step_tddft',mt_int,1, step_tddft))
     $     step_tddft = 1
c
c     are we performing any frequency calculations along the ground state trajectory ?
      if (.not.rtdb_get(rtdb,'qmd:freq',mt_log,1,do_freq))
     $     do_freq=.false.
c
c     calculate frequency calculations every step_freq along the ground state trajectory
      if (.not.rtdb_get(rtdb,'qmd:step_freq',mt_int,1, step_freq))
     $     step_freq = 1
c
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
