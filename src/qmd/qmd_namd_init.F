c
c     qmd_namd_init: initialize NAMD parameters
c     
      subroutine qmd_namd_init(rtdb,geom_p,g_vecs,g_vecs_p,l_zden,
     &                         k_zden,l_coup,k_coup,l_coup_p,k_coup_p,
     &                         l_en,k_en,l_en_p,k_en_p,l_so,k_so,l_sv,
     &                         k_sv,state,state_p,nstates,dt_elec,
     &                         nroots,ipol,nocc,nmo,ks_spin,
     &                         do_tdks,tda,deco)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "geom.fh"
#include "inp.fh"
c
#include "qmd_common.fh"
c
      integer rtdb
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
      double precision dt_elec           ! electronic time step
      logical do_tdks                    ! TDKS flag
      logical tda                        ! TDA flag
      logical deco                       ! decoherence flag
c
      integer basis                      ! basis handle
      integer ispin
      logical triplet
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      character*32 pname
c
      pname = "qmd_namd_init "
c
c     create geometry object for previous step geometry
      if (.not. geom_create(geom_p,'geometry_prev'))
     &   call errquit(pname//'geom_create ',geom_p, GEOM_ERR)
c
c     load in the geometry
      if (.not. geom_rtdb_load(rtdb,geom_p,'geometry'))
     &   call errquit(pname//'geom_rtdb_load ', geom_p, RTDB_ERR)
c
c     electronic time step
      if (.not.rtdb_get(rtdb,'qmd:dt_elec',mt_dbl,1,dt_elec))
     & call errquit(pname//'failed to read dt_elec ',0,RTDB_ERR)
c
c     decoherence?
      if (.not.rtdb_get(rtdb,'qmd:deco',mt_log,1,deco))
     $   call errquit(pname//'deco flag ',0,RTDB_ERR)
c
c     are we doing TDKS?
      if (.not.rtdb_get(rtdb,'qmd:tdks',mt_log,1,do_tdks))
     $   call errquit(pname//'tdks flag ',0,RTDB_ERR)
c
c     get number of states
      if (.not.rtdb_get(rtdb,'qmd:nstates',mt_int,1,nstates))
     $   call errquit(pname//'# of states not given ',0,RTDB_ERR)
c
c     get initial state (will be over-ridden if restart file exists)
      if (.not.rtdb_get(rtdb,'qmd:init_state',mt_int,1,state))
     $   call errquit(pname//'no initial state ',0,RTDB_ERR)
      state_p=state
c
      if (ga_nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*)
        call util_print_centered
     &   (LuOut,'QMD NAMD Parameters',40,.true.)
        write(LuOut,*)
        write(LuOut,*)
        write(LuOut,110) "Elec. time step:", dt_elec
        write(LuOut,120) "Number of states:", nstates
c    For input/output state runs 0:nstates-1
c     Inside code states runs 1:nstates 
        write(LuOut,120) "Initial state:", state-1
        write(LuOut,130) "Decoherence:", deco
      end if
c
      if (do_tdks) then
c       1=alpha,2=beta
        if (.not.rtdb_get(rtdb,'qmd:ks_spin',mt_int,1,ks_spin))
     $     call errquit(pname//'no spin specified ',0,RTDB_ERR)
        if ((ks_spin.gt.2).or.(ks_spin.lt.1))
     $     call errquit(pname//'invalid spin specified ',0,RTDB_ERR)
        if (ga_nodeid().eq.0) then
          write(LuOut,140) "Will do TDKS"
          if (ks_spin.eq.1) then
            write(LuOut,140) "TDKS spin channel: alpha"
          else
            write(LuOut,140) "TDKS spin channel: beta"
          end if
        end if
c
c       set up arrays for orbital indices for swaps
        if (.not.ma_push_get(mt_int,nstates,'swap occ',l_so,k_so))
     $     call errquit(pname//'ma_push_get swap occ ',nstates,MA_ERR)
        if (.not.ma_push_get(mt_int,nstates,'swap vir',l_sv,k_sv))
     $     call errquit(pname//'ma_push_get swap vir ',nstates,MA_ERR)
      else
c
c       get number of TDDFT roots calculated
        if (.not.rtdb_get(rtdb,'tddft:nroots',mt_int,1,nroots))
     &    call errquit(pname//'failed to get nroots ',0,RTDB_ERR)
        if (nroots.lt.(nstates-1))
     &    call errquit(pname//'nroots<(nstates-1) ',0,RTDB_ERR)
c
c       TDA or RPA
        if (.not.rtdb_get(rtdb,'tddft:tda',mt_log,1,tda))
     &       tda = .false.
      end if
      if (ga_nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*)
        call util_flush(LuOut)
      end if
c
c     create array for density matrix
      if (.not.ma_push_get(mt_dcpl,nstates**2,'zden',l_zden,k_zden))
     $   call errquit(pname//'ma_push_get density ',nstates**2,MA_ERR)
c
c     create array for coupling matrix at current step
      if (.not.ma_push_get(mt_dbl,nstates**2,'coup',l_coup,k_coup))
     $   call errquit(pname//'ma_push_get coup ',nstates**2,MA_ERR)
c
c     create array for coupling matrix at previous step
      if (.not.ma_push_get(mt_dbl,nstates**2,'coup_p',l_coup_p,
     $                     k_coup_p))
     $   call errquit(pname//'ma_push_get coup_p ',nstates**2,MA_ERR)
c
c     create array for energies at current step
      if (.not.ma_push_get(mt_dbl,nstates,'energies',l_en,k_en))
     $   call errquit(pname//'ma_push_get energies ',nstates,MA_ERR)
c
c     create array for energies at previous step
      if (.not.ma_push_get(mt_dbl,nstates,'energies_p',l_en_p,k_en_p))
     $   call errquit(pname//'ma_push_get energies_p ',nstates,MA_ERR)
c
c     get ipol; 1=RHF, 2=UHF
      if (.not.rtdb_get(rtdb,'dft:ipol',mt_int,1,ipol))
     &  call errquit(pname//'failed to get ipol ',0,RTDB_ERR)
c
c     make sure we are not doing a trivial calculation
      if ((.not.do_tdks).and.(ipol.eq.1)) then
        if (.not.rtdb_get(rtdb,'tddft:triplet',mt_log,1,triplet))
     &    call errquit(pname//'failed to get triplet ',0,RTDB_ERR)
        if (triplet) 
     &    call errquit(pname//'triplets do not couple ',0,CALC_ERR)
      end if
c
c     currently must use unrestricted with TDKS
      if (do_tdks.and.(ipol.eq.1))
     &  call errquit(pname//'use unrestricted with TDKS ',0,RTDB_ERR)
c
c     create GAs for MOvecs
      if (.not.bas_create(basis,'ao basis'))
     &     call errquit(pname//'failed to create ',0,BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom_p,basis,'ao basis'))
     &     call errquit(pname//'failed to load ',0,BASIS_ERR)
c
      do ispin=1,ipol
        g_vecs(ispin)=ga_create_atom_blocked(geom_p,basis,'MOvecs')
        g_vecs_p(ispin)=ga_create_atom_blocked(geom_p,basis,'MOvecs_p')
      end do
c
      if (.not.bas_destroy(basis))
     $  call errquit(pname//'failed to destroy basis ',basis,BASIS_ERR)
c
c     set initial density matrix (will be replaced if restart file exists)
      call qmd_set_den(nstates,state,dcpl_mb(k_zden))
c
c     force theory
      if ((state.eq.1).or.do_tdks) then
        if (.not.rtdb_cput(rtdb,'task:theory',1,'dft'))
     &    call errquit(pname//'failed to set theory ',0,RTDB_ERR)
      else if ((state.gt.1).and.(state.le.nstates).and.
     &         (.not.do_tdks)) then
        if (.not.rtdb_cput(rtdb,'task:theory',1,'tddft'))
     &    call errquit(pname//'failed to set theory ',0,RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tddft_grad:iroots',mt_int,1,
     $                    state-1))
     $    call errquit(pname//'failed to write iroot ',0,RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tddft_grad:isinglet_roots',
     $                    mt_int,1,state-1)) 
     $    call errquit(pname//'failed to write iroot ',0,RTDB_ERR)
        if (.not. rtdb_put(rtdb,'tddft:target',mt_int,1,state-1))
     $    call errquit(pname//'failed to write target ',0,RTDB_ERR)
      else
        call errquit(pname//'initial state not valid ',0,0)
      end if
c
      call ga_sync()
c
  110 format(A25,F20.4)
  120 format(A25,I20)
  130 format(A25,L20)
  140 format(A25)
c
      return
      end
