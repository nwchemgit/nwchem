c
c     qmd_driver: main driver
c
      logical function qmd_driver(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "stdio.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer rtdb
c
      integer geom                    ! geometry object
c
      integer l_r, k_r                ! coordinates
      double precision r(3)
      integer l_c, k_c                ! centroid of geometry
c
      integer l_v, k_v                ! velocities
      double precision v(3)
c
      integer l_g, k_g                ! gradients
      integer l_m, k_m                ! masses
c
      integer l_act, k_act            ! active atoms
c
      integer nat                     ! number of atoms
      integer nactive                 ! number of active atoms
      integer nstep_nucl              ! total nuclear steps
      integer sstep_nucl              ! first nuclear step
      integer istep_nucl              ! nuclear md step
      integer com_step                ! how often to remove trans/rot
      integer print_xyz               ! how often to print trajectory
      integer step_property           ! how often to calculate properties
      integer step_tddft              ! how often to calculate tddft energies along the trajectory
      integer step_freq               ! how often to calculate frequencies along the trajectory
      double precision dt_nucl        ! nuclear md time step
      double precision eltime         ! elapsed time (in fs)
      double precision esys           ! system energy (from task_gradient)
      double precision ekin           ! kinetic energy from the md part
      double precision etotal         ! total energy: esys + ekin
      double precision ekin_targ      ! target kinetic energy
      double precision dipole(3)      ! electric dipole moment
c
      logical do_cavity               ! include cavity constraints
      double precision cavity_axis(3) ! cavity axes
      double precision cavity_spring  ! cavity force spring constant
c
c     Surface Hopping
      integer geom_p                     ! geometry, previous step
      integer g_vecs(2), g_vecs_p(2)     ! movecs
      integer l_zden, k_zden             ! density matrix handle
      integer l_coup, k_coup             ! coupling matrix, current step
      integer l_coup_p, k_coup_p         ! coupling matrix, previous step
      integer l_en, k_en                 ! energies, current step
      integer l_en_p, k_en_p             ! energies, previous step
      integer l_so, k_so                 ! index of occupied orbitals for swaps
      integer l_sv, k_sv                 ! index of virtual orbitals for swaps
      integer nstates                    ! number of states
      integer state                      ! current state
      integer state_p                    ! previous state
      integer nroots                     ! number of TDDFT roots
      integer ipol                       ! 1=RHF, 2=UHF
      integer nocc(2)                    ! number of occupied orbitals
      integer nmo(2)                     ! number of MOs
      integer ks_spin                    ! 1=alpha, 2=beta
      integer ispin
      double precision dt_elec           ! electronic time step
      double precision kbt               ! k_b*T in au (for plumed)
      double precision wall              ! wall time
      logical do_namd                    ! NAMD flag
      logical do_property                ! property flag
      logical do_tddft                   ! tddft flag
      logical do_freq                    ! frequency flag
      logical do_tdks                    ! TDKS flag
      logical tda                        ! TDA flag
      logical deco                       ! decoherence flag
c
      logical status_qmd
      logical status_property
      logical status_tddft
      logical status_freq
      logical do_none
      logical dip_exists
      logical restart
c
      character*32 pname
      character*32 ext_forces
      integer xyzfileunit
      character*255 xyzfilename
      integer rstfileunit
      character*255 rstfilename
c
      logical geom_centroid
      external geom_centroid
c
      logical qmd_hasproperty
      external qmd_hasproperty
c
      logical tddft_energy
      external tddft_energy
c
      logical task_freq
      external task_freq
c
      logical qmd_property,qmd_linear
      external qmd_property,qmd_linear
c
c     preliminaries
      status_qmd = .true.
      status_property = .true.
      status_tddft = .true.
      status_freq = .true.
      dip_exists=.true.
      pname = "qmd_driver: "
      esys = 0.d0
      ekin = 0.d0
      etotal = 0.d0
      dipole(1) = 0.d0
      dipole(2) = 0.d0
      dipole(3) = 0.d0
c
c     check if we need cavity constraints
      if(.not.rtdb_get(rtdb,'qmd:cavity',mt_log,1,do_cavity))
     &     do_cavity = .false.
      if (do_cavity) then
      if(.not.rtdb_get(rtdb,'qmd:cavity_axis',mt_dbl,3,cavity_axis))
     &    call errquit("cannot get cavity axes",0,0)
      if(.not.rtdb_get(rtdb,'qmd:cavity_spring',mt_dbl,1,cavity_spring))
     &    call errquit("cannot get cavity spring constant",0,0)
      end if ! do_cavity
c
c     check if we want external forces
c
      ext_forces='none'
      if (.not.rtdb_cget(rtdb,'qmd:ext_forces',1,ext_forces))
     E       ext_forces='none'
      if(ga_nodeid().eq.0) write(6,*) ' Ext forces = ',
     E     ext_forces(1:inp_strlen(ext_forces))
c
c     get xyzfileunit
      if(.not.util_get_io_unit(xyzfileunit))
     &   call errquit("cannot get file unit",0,0)
c
c     get xyzfilename
      call util_file_name('xyz',.false.,.false.,xyzfilename)
c
c     get rstfileunit
      if(.not.util_get_io_unit(rstfileunit))
     &   call errquit("cannot get file unit",0,0)
c
c     get rstfilename
      call util_file_name('qmdrst',.false.,.false.,rstfilename)
c
c     check defaults
c
      call qmd_rdinput(rtdb)
c     initialize parameters
      call qmd_init(rtdb,nstep_nucl,dt_nucl,do_none,com_step,print_xyz,
     &   restart,do_namd,step_property,do_property,step_tddft,do_tddft,
     &   step_freq,do_freq)
c
c     create geometry object 
      if (.not. geom_create(geom,'geometry'))
     &   call errquit(pname//'geom_create', 911, GEOM_ERR)
c
c     load in the geometry
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     &   call errquit(pname//'geom_rtdb_load', geom, RTDB_ERR)
c
c     get number of atoms
      if (.not. geom_ncent(geom,nat))
     $   call errquit(pname//'geom_ncent',geom, GEOM_ERR)
      if (nat.eq.0) return
c
c     local memory and initialization
      if(.not.ma_push_get(mt_dbl,nat,'mass',l_m,k_m))  ! mass
     $   call errquit(pname//'ma_push_get mass', nat,MA_ERR)
      call dfill(nat, 0.0d0, dbl_mb(k_m), 1)
c
      if(.not.ma_push_get(mt_dbl,3*nat,'coords',l_r,k_r)) ! coords
     $   call errquit(pname//'ma_push_get coords', 3*nat,MA_ERR)
      call dfill(3*nat, 0.0d0, dbl_mb(k_r), 1)
c
      if(.not.ma_push_get(mt_dbl,3,'centroid',l_c,k_c))  ! centroid
     $   call errquit(pname//'ma_push_get centroid', 3,MA_ERR)
      call dfill(3, 0.0d0, dbl_mb(k_c), 1)
c
      if(.not.ma_push_get(mt_dbl,3*nat,'vels',l_v,k_v))  ! vels
     $   call errquit(pname//'ma_push_get vels', 3*nat,MA_ERR)
      call dfill(3*nat, 0.0d0, dbl_mb(k_v), 1)
c
      if(.not.ma_push_get(mt_dbl,3*nat,'grads',l_g,k_g))  ! grads
     $   call errquit(pname//'ma_push_get grads',3*nat,MA_ERR)
      call dfill(3*nat, 0.0d0, dbl_mb(k_g), 1)
c
c     get the active atoms (also need a check for linearity)
c     need to implement constraints
      if (.not. ma_push_get(mt_log,nat,'active atoms',l_act,k_act)) ! active
     $   call errquit(pname//'grad_active_atoms',nat,MA_ERR)
      call grad_active_atoms(rtdb, nat, log_mb(k_act), nactive)
c
c     set degrees of freedom depending on the type of thermostat
      ndeg = 3*nactive ! initialize ndeg for local thermostats (e.g. langevin)
      if (.not.do_langevin) then  ! for global thermostats
         if(.not.do_linear) then
c     check linearity of molecular geometry
            do_linear=qmd_linear(geom,nat,dbl_mb(k_r))
               write(6,*) ' do_linear ',do_linear
         endif
       if (do_linear) then    ! check for linear flag
        ndeg = 3*nactive-5   ! 3N-5 (linear)
       else
        ndeg = 3*nactive-6   ! 3N-6 (general)
       end if ! do_linear
      end if   ! langevin check
c
      if (ndeg.le.0)
     $   call errquit(pname//'ndeg<=0',geom,GEOM_ERR)
c
c     get the masses
      if (.not.geom_masses_get(geom,nat,dbl_mb(k_m)))  ! masses
     $   call errquit(pname//'geom_masses_get',geom,GEOM_ERR)
c
c     Are we doing surface hopping?
c      ***Surface hopping only works for DFT/TDDFT currently***
c              ***Forces DFT/TDDFT calculation***
      if (do_namd) then
        call qmd_namd_init(rtdb,geom_p,g_vecs,g_vecs_p,l_zden,k_zden,
     &                     l_coup,k_coup,l_coup_p,k_coup_p,l_en,k_en,
     &                     l_en_p,k_en_p,l_so,k_so,l_sv,k_sv,state,
     &                     state_p,nstates,dt_elec,nroots,ipol,
     &                     nocc,nmo,ks_spin,do_tdks,tda,deco)
        if (mod(nint(dt_nucl*1.d8),nint(dt_elec*1.d8)).ne.0)
     &   call errquit(pname//'dt_nucl and dt_elec not compatible',0,
     &                INPUT_ERR)
      end if

      if(ext_forces.ne.'none')
     $     call qmd_ext_init(nat,dt_nucl,kb_au*targ_temp,ext_forces)
c
c     get initial coordinates and velocities
      call qmd_coords_vels(rtdb,geom,nat,sstep_nucl,ekin,rstfilename,
     $       rstfileunit,dbl_mb(k_r),dbl_mb(k_v),dbl_mb(k_m),
     $       log_mb(k_act),nactive,restart,do_namd,state,nstates,
     $       dbl_mb(k_coup_p),dcpl_mb(k_zden))

      if (do_namd) state_p=state
c
      call qmd_kinetic(nat,dbl_mb(k_m),dbl_mb(k_v),ekin)
      if (ga_nodeid().eq.0) then
        final_temp = 2.0d0*ekin/kb_au/(1.d0*ndeg)
        write(luout,120) "Current temp. (K):",final_temp
      end if
c
c     target kinetic energy (only used if thermostat is rescale)
      ekin_targ=(1.d0*ndeg)*kb_au*targ_temp/2.d0
c
      if (do_namd.and.do_tdks) then
c       get energies and movecs
        nmo=-1
        call qmd_tdks_enmo(rtdb,geom,state,nstates,ipol,nmo,nocc,
     $                     int_mb(k_so),int_mb(k_sv),dbl_mb(k_en_p),
     $                     g_vecs_p,ks_spin)
      end if
c
c     perform energy gradient calculation on initial coordinates
      call qmd_gradient(rtdb,nat,0,
     $     dbl_mb(k_m),dbl_mb(k_g),dbl_mb(k_r),esys,ext_forces)
c
      if (do_namd.and.(.not.do_tdks)) then
c       get energies and movecs
        nmo=-1
        call qmd_tddft_enmo(rtdb,geom,state,nstates,nroots,ipol,nmo,
     $                      nocc,dbl_mb(k_en_p),g_vecs_p)
c       rename civecs file so we have it later
        call qmd_move_civecs(rtdb,ipol)
      end if
c
c     calculate the centroid of the structure
      if (.not.geom_centroid(geom,dbl_mb(k_c)))
     $      call errquit(pname//'geom_centroid',0,GEOM_ERR)
c
c     include cavity constraints, if requested
      if (do_cavity) then
       call qmd_cavity(nat,cavity_axis,cavity_spring,dbl_mb(k_c),
     &       dbl_mb(k_r),dbl_mb(k_g))
      end if ! do_cavity
c
c     if starting fresh write initial structure to xyz
      if (.not.restart) then
c
c       pull dipole from RTDB
        if(.not.rtdb_get(rtdb,'task:dipole',mt_dbl,3,dipole))
     &     dip_exists=.false.
c
c       instantaneous total energy: instantaneous kinetic energy + system energy
        etotal = ekin + esys
c
        if (ga_nodeid().eq.0) then
c
c         write trajectory
          call qmd_write_trj(xyzfilename,xyzfileunit,.false.,nat,
     &           0,etotal,geom,dipole)
        end if
      end if  ! restart
      wall = util_wallsec()
c
c     start nuclear molecular dynamics loop
      do istep_nucl = sstep_nucl,nstep_nucl
c
c       don't attempt to remove COM rot/tran if we have inactive atoms
        if (mod(istep_nucl,com_step).eq.0.and.nat.eq.nactive) then
c         remove global translations and rotations
c         kinetic energy is the same on exit as was on entry
          call qmd_fix_momentum(nat,dbl_mb(k_m),dbl_mb(k_v),dbl_mb(k_r))
        end if
c
        if (.not.do_none) then
c
c         apply thermostat to velocities
          call qmd_thermostat(nat,dt_nucl,dbl_mb(k_m),dbl_mb(k_v),
     &           ekin,ekin_targ,log_mb(k_act))
        end if ! do_none
c
c       update velocity at half step (dt/2)
        call qmd_nucl_vel(nat,0.5d0*dt_nucl,dbl_mb(k_m),dbl_mb(k_v),
     &         dbl_mb(k_g))
c
c       update position at dt
        call qmd_nucl_pos(nat,dt_nucl,dbl_mb(k_m),dbl_mb(k_r),
     &        dbl_mb(k_v))
c
c       update geom and store in rtdb
        call qmd_geom_store(rtdb,geom,nat,dbl_mb(k_r))
c
        if (do_namd.and.do_tdks) then
c         get energies and movecs
          call qmd_tdks_enmo(rtdb,geom,state,nstates,ipol,nmo,nocc,
     $                       int_mb(k_so),int_mb(k_sv),dbl_mb(k_en),
     $                       g_vecs,ks_spin)
        end if
c
c       calculate energy and gradient
        call qmd_gradient(rtdb,nat,istep_nucl,
     $       dbl_mb(k_m),dbl_mb(k_g),dbl_mb(k_r),esys,ext_forces)
c
c       include cavity reflection if requested
        if (do_cavity) then
         call qmd_cavity(nat,cavity_axis,cavity_spring,dbl_mb(k_c),
     &         dbl_mb(k_r),dbl_mb(k_g))
        end if ! do_cavity
c
c       pull dipole from RTDB
        if(.not.rtdb_get(rtdb,'task:dipole',mt_dbl,3,dipole))
     &     dip_exists=.false.
c
c       update velocity at half step (dt/2)
        call qmd_nucl_vel(nat,0.5d0*dt_nucl,dbl_mb(k_m),dbl_mb(k_v),
     &         dbl_mb(k_g))
c
        if (.not.do_none) then
c         apply thermostat to velocities
          call qmd_thermostat(nat,dt_nucl,dbl_mb(k_m),dbl_mb(k_v),
     &           ekin,ekin_targ,log_mb(k_act))
        end if ! do_none
c
c       persist updated velocity in geom
        if (.not.geom_vel_set(geom,dbl_mb(k_v)))
     $   call errquit(pname//'geom_vel_set',0,GEOM_ERR)
c
        if (do_namd) then
          if (.not.do_tdks) call qmd_tddft_enmo(rtdb,geom,state,nstates,
     $                                          nroots,ipol,nmo,nocc,
     $                                          dbl_mb(k_en),g_vecs)
          esys=dbl_mb(k_en+state-1)
c     check nroots vs nstates (tddft might have changed nroots)
c     get number of TDDFT roots calculated
        if (.not.rtdb_get(rtdb,'tddft:nroots',mt_int,1,nroots))
     &    call errquit(pname//'failed to get nroots ',0,RTDB_ERR)
        if (nroots.lt.(nstates-1)) then
c     call errquit(pname//'nroots<(nstates-1) ',0,RTDB_ERR)
           if(ga_nodeid().eq.0) write(luout,124) nstates,nroots
 124       format(9x,'Number of states decreased from ',i8,' to ',i8)
           nstates=nroots
        endif

c         calculate couplings
          call qmd_coup(rtdb,geom,geom_p,ipol,nstates,dt_nucl,nmo,
     $                  nocc,g_vecs,g_vecs_p,dbl_mb(k_coup),
     $                  int_mb(k_so),int_mb(k_sv),dbl_mb(k_en),tda,
     $                  do_tdks,ks_spin)
c
c         propagate electronic system
          call qmd_elec_prop(rtdb,nat,dcpl_mb(k_zden),dbl_mb(k_en_p),
     $                       dbl_mb(k_en),dbl_mb(k_coup_p),
     $                       dbl_mb(k_coup),nstates,dbl_mb(k_v),
     $                       dbl_mb(k_m),ekin,esys,state,dt_nucl,
     $                       istep_nucl,dt_elec,do_tdks,deco)
c
c         calculate gradient for new current state
          if (state_p.ne.state) then
            if (do_tdks) then
              call qmd_tdks_swap(rtdb,state,nstates,ks_spin,
     $                           int_mb(k_so),int_mb(k_sv))
            end if
            call qmd_gradient(rtdb,nat,istep_nucl,
     $           dbl_mb(k_m),dbl_mb(k_g),dbl_mb(k_r),esys,ext_forces)
          end if
          state_p=state
c
          if (.not.do_tdks) then
c           rename civecs file so we have it later
            call qmd_move_civecs(rtdb,ipol)
          end if
        end if !do_namd
c
c       calculate current kinetic energy
        call qmd_kinetic(nat,dbl_mb(k_m),dbl_mb(k_v),ekin)
c       instantaneous total energy: instantaneous kinetic energy + system energy
        etotal = ekin + esys
c
c       instantaneous temperature
        final_temp = 2.0d0*ekin/kb_au/(1.d0*ndeg)
c       print md run information
        if (ga_nodeid().eq.0) then
c
         write(luout,*)
         write(luout,*)
         call util_print_centered(LuOut,'QMD Run Information',20,.true.)
c         write(luout,110) "Current time step:",istep_nucl
         eltime=dble(istep_nucl)*dt_nucl*cau2fs
         write(luout,130) "Time elapsed (fs) :",eltime
         write(luout,135) "Kin. energy (a.u.):",istep_nucl,ekin
         write(luout,135) "Pot. energy (a.u.):",istep_nucl,esys
         write(luout,135) "Tot. energy (a.u.):",istep_nucl,etotal
         write(luout,125) "Target temp. (K)  :",istep_nucl,targ_temp
         write(luout,125) "Current temp. (K) :",istep_nucl,final_temp
         if (do_namd) then
c          inside code state runs from 1 to nstates
c          for input and output, state runs from 0 to nstates-1
           write(luout,145) "Current state     :",istep_nucl,state-1
         end if
         if (dip_exists) then
           write(luout,140) "Dipole (a.u.)     :",istep_nucl,dipole(:)
         else
           write(luout,'(A21)') "No dipole on RTDB"
         end if
         write(luout,125) "Wall time/ts    (s) :",istep_nucl,
     ,        util_wallsec()-wall
         wall = util_wallsec()
         write(luout,*)
         write(luout,*)
c
c        check how often to write out the tracjectory
         if (mod(istep_nucl,print_xyz).eq.0) then
c
c          write trajectory
           call qmd_write_trj(xyzfilename,xyzfileunit,.true.,nat,
     &            istep_nucl,etotal,geom,dipole)
         end if
c
c        write restart
         call qmd_write_rst(rstfilename,rstfileunit,nat,istep_nucl,
     &                      ekin,geom,state,nstates,dcpl_mb(k_zden),
     &                      dbl_mb(k_coup),do_namd)
c
        end if  ! ga_nodeid()
c
c       calculate molecular properties
        if (do_property) then
         if (mod(istep_nucl,step_property).eq.0) then
c            if (ga_nodeid().eq.0) write(luout,*) "calling qmd_property"
             status_property = qmd_property(rtdb)
         end if
        end if  ! do_property
c
c       calculate tddft spectra along ground state trajectory
        if (do_tddft) then
         if (mod(istep_nucl,step_tddft).eq.0) then
c            if (ga_nodeid().eq.0) write(luout,*) "calling tddft_energy"
             status_tddft = tddft_energy(rtdb)
         end if 
        end if  ! do_tddft
c
c       calculate frequencies along ground state trajectory
        if (do_freq) then
         if (mod(istep_nucl,step_freq).eq.0) then
c            if (ga_nodeid().eq.0) write(luout,*) "calling tddft_energy"
             status_freq = task_freq(rtdb)
         end if 
        end if  ! do_freq
c
c       restart file has been written so set restart flag in RTDB
        if (.not.rtdb_put(rtdb,'qmd:restart',mt_log,1,.true.))
     &    call errquit(pname//'failed to write restart',0,RTDB_ERR)
c
      end do  ! istep_nucl
  110 format(A20,I20)
  120 format(A25,F20.2)
  125 format(A21,I9,F20.2)
  130 format(A21,F29.6)
  135 format(A21,I9,F20.6)
  140 format(A21,I9,3ES15.6)
  145 format(A21,I9,I20)
c
c     finalized ext forces
c
      if(ext_forces.ne.'none') call qmd_ext_final(ext_forces)
c
c     clear memory
      if(.not.ma_chop_stack(l_m))
     + call errquit(pname//'ma_chop_stack m',0,MA_ERR)
      if(.not.ma_verify_allocator_stuff())
     + call errquit(pname//'ma_verify_allocator_stuff',0,MA_ERR)
      if(.not.geom_destroy(geom))
     + call errquit(pname//'geom_destroy',0,GEOM_ERR)
      if (do_namd) then
       if(.not.geom_destroy(geom_p))
     +  call errquit(pname//'geom_destroy',0,GEOM_ERR)
       do ispin=1,ipol 
        if(.not.ga_destroy(g_vecs(ispin)))
     +    call errquit(pname//'ga_destroy',g_vecs(ispin),GA_ERR)
        if(.not.ga_destroy(g_vecs_p(ispin)))
     +    call errquit(pname//'ga_destroy',g_vecs_p(ispin),GA_ERR)
       end do
      end if
c
      qmd_driver = status_qmd
c
      end
c $Id$
      logical function qmd_linear(geom,n,xyz)
      implicit none
      integer n,geom
      double precision xyz(3,*)
c
      integer i
      double precision angle
      logical  geom_calc_angle,geom_cart_coords_get
      external geom_calc_angle,geom_cart_coords_get
c      
      qmd_linear=.true.
      if(n.eq.2) return
      if(.not.geom_cart_coords_get(geom, xyz))
     c        call errquit(' geom_cart_coords fails ',0,0)
      do i=1,n-2
         if(.not.geom_calc_angle(
     1        xyz(1,i),
     2        xyz(1,i+1),
     3        xyz(1,i+2),
     a        angle))
     c        call errquit(' geom_calc_angle fails ',0,0)
         if(abs(angle-180d0).gt.1d-5) qmd_linear=.false.
      enddo
      return
      end
