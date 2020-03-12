C
C     rt_tddft_params_init_common.F
C
C     Computes parameters which are universal to closed shell,
C     openshell, etc and stores them in the params struct.  Also loads
C     rtdb parameters into params.
C
      subroutine rt_tddft_init_common (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t), intent(inout) ::  params


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_init_common: "


C     == Variables ==
      integer n
      logical ok
      integer me

      double precision local_d_conv

      double precision center(3)
      logical dft_converged

      integer filesize, memsize

C     (tmp dip and quad integral matricies; all nbf_ao x nbf_ao)
      integer g_dip_x, g_dip_y, g_dip_z
      integer g_quad_xx, g_quad_xy, g_quad_xz
      integer g_quad_yy, g_quad_yz, g_quad_zz

      double precision elapsed
      integer g_fock_ao_core


      me = ga_nodeid ()

      
C
C     Hardcoded parameters.
C
      params%tol_zero = 1d-8

      
C
C     Get density tolerance--dictates precision of Fock matrix calculations.
C
      if (.not. rtdb_get (params%rtdb, "dft:d_conv",
     $     mt_dbl, 1, local_d_conv)) then
         
         call errquit (pname//"Failed to read d_conv", 0, 0)
      else
         params%d_conv = local_d_conv
      endif



C
C     Check that DFT SCF converged from before.
C
      if (.not. rtdb_get (params%rtdb, "dft:converged",
     $     mt_log, 1, dft_converged)) then
         
         call rt_tddft_print_warning (
     $        "Failed to find previous SCF convergence status."//
     $        " Either you are loading vectors from a previous job"//
     $        " or database is corrupt--be careful.")
         
      else
         
         if (.not. dft_converged) call rt_tddft_print_warning (
     $        "SCF did not converge previously--it is quite likely"//
     $        " that the starting density is not in the ground state")
         
      endif


C
C     Store params passed in from SCF or in headers.
C
      params%nbf_ao = nbf_ao
      params%ipol = ipol
      params%ao_bas_han = ao_bas_han
      params%cd_bas_han = cd_bas_han


C
C     Get values from rtdb and store in params struct.
C      
      call rt_tddft_init_rtdb (params)
 


C     XXX DISABLED FOR NOW
#if 0      
C
C     Disable *all* disk I/O for DFT routines (default behavior).
C
      if (params%nodisk) then
         
C     (dont use disk for gridpts)
         if(.not.rtdb_put(params%rtdb,'dft:largenode',mt_log,1,.true.))
     $        call errquit(pname//"rtdb_put largenode failed",
     $        0, RTDB_ERR)
         

C     (direct)
         if (.not. rtdb_put(params%rtdb, 'dft:direct',mt_log,1,.true.))
     $        call errquit(pname//"rtdb_put direct failed",0,RTDB_ERR)
         
         filesize = -1
         memsize  = -1
         if (.not.rtdb_put(params%rtdb,
     $        'int2e:filesize',mt_int,1,filesize))
     $        call errquit(pname//'failed to store filesize',
     $        0, RTDB_ERR)
         
         if (.not.rtdb_put(params%rtdb,
     $        'int2e:memsize',mt_int,1,memsize))
     $        call errquit(pname//'failed to store filesize',
     $        0, RTDB_ERR)
         
         
C     (no I/O)
         if (.not. rtdb_put(params%rtdb, 'dft:noio', mt_int, 1, .true.))
     $        call errquit(pname//"rtdb_put noio failed", 0, RTDB_ERR)

         if (me.eq.0) then
            write (luout, *) "*** All disk I/O has been disabled ***"
         endif

      else
         
         call rt_tddft_print_warning (
     $        "Disk I/O using DFT options; proceed with caution")
            
      endif
#endif

      

C
C     Dummy values XXX NEVER USED?
C     
CXXX  [KAL]: add check to ensure none of these are used
      params%g_vxc(1) = -999999
      params%g_vxc(2) = -999999
      params%g_vxc(3) = -999999
      params%g_vxc(4) = -999999
      params%g_xcinv  = -999999
      params%rho_n = -99d0
      params%lzorasf = .false.
      
C
C     Nuclear repulsion energy
C
      if (.not. geom_nuc_rep_energy (geom, params%Enuc))
     $     call errquit (pname//'Failed to compute Enuc', 0, BASIS_ERR)

      
C
C     I/O logical flag (used in dft call)
C     this might have been set in the override i/o block above.
C
      if (noio.eq.1) then
         params%iolgc = .false.
      else
         params%iolgc = .true.            
      endif
 

C
C     Charge density fitting flag
C     
      if (cd_bas_han.eq.-99999)then
         params%iVcoul_opt = 0
      else
         params%iVcoul_opt = 1
      endif
     

C
C     XC fitting flag
C
      if (xc_bas_han.eq.-99999)then
         params%iVxc_opt = 0
      else
         params%iVxc_opt = 1
      endif


C     
C     Check parameters.
C     
      call rt_tddft_init_check_params (params)


C
C     Calculate the total number of time steps.  Also, if number of
C     prints, checks, snapshots, etc is greater than nt, set them to nt.
C
      if (params%nprints < 0)
     $     call errquit (pname//"nprints must be >= 0", 0, 0)

      if (params%nchecks < 0)
     $     call errquit (pname//"nchecks must be >= 0", 0, 0)

c$$$      if (params%nsnapshots < 0)
c$$$     $     call errquit (pname//"nsnapshots must be >= 0", 0, 0)

      if (params%nrestarts < 0)
     $     call errquit (pname//"nrestarts must be >= 0", 0, 0)


      params%nt = ceiling((params%tmax-params%tmin)/params%dt)
      
      if (params%nprints > params%nt) params%nprints = params%nt
      if (params%nchecks > params%nt) params%nchecks = params%nt
c$$$      if (params%nsnapshots > params%nt) params%nsnapshots = params%nt
      if (params%nrestarts > params%nt) params%nrestarts = params%nt


C     hardcoded for now:
      params%nestims = params%nt / 100  !i.e., every 100 steps



C
C     Store geom tag for "active" geometry.
C
      params%geom_active_handle = geom


C
C     Center of mass.
C
      if (.not.geom_center_of_mass(params%geom_active_handle,center))
     &     call errquit (pname//'could not get center of mass',
     $     0, GEOM_ERR)

      params%origin(1) = center(1)
      params%origin(2) = center(2)
      params%origin(3) = center(3)
      

C
C     Number of exchange correlation energies.
C
      params%nExc = idecomp + 1

      
C
C     Number of atoms.
C
      if (.not. geom_ncent(params%geom_active_handle, params%natoms))
     &     call errquit(pname//'geom_ncent failed',73, GEOM_ERR)



C
C     Initialize all integrals.
C
      call rt_tddft_init_fock2e (params)
      call rt_tddft_init_dip_quad_ints (params)
      if (cdfit) call rt_tddft_init_coulcdfit (params)


C
C     Initialize geometries.
C      
      call rt_tddft_init_geoms (params)


C
C     Compute and store 1e (time-dependent) core part of the Fock
C     matrix.  Note, the GA created here will be destroyed in
C     rt_tddft_clean().
C
       if (params%prof) call prof_start (elapsed)

C     xxx check ns_ao vs nbf_ao
      if (.not. ga_create (mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "fock_ao_core", 0, 0, g_fock_ao_core))
     $     call errquit (pname//"create fock_ao_core failed",0,GA_ERR)
      
      call ga_zero (g_fock_ao_core)
      call int_1e_ga(ao_bas_han, ao_bas_han,
     $     g_fock_ao_core, 'kinetic', oskel)
      call int_1e_ga(ao_bas_han, ao_bas_han,
     $     g_fock_ao_core, 'potential', oskel)

      params%g_fock_ao_core = g_fock_ao_core

      if (params%prof) call prof_end(elapsed,"Fock core initialization")


      
      if (params%static)
     $     call rt_tddft_print_warning ("Static calculation--"//
     $     "Fock matrix will be frozen")


C
C     Spatial complex absorbing potential (CAP).
C
      if (params%cap_active) then
         if (.not. ga_create(mt_dcpl, nbf_ao, nbf_ao ,"zCAP",
     $        0, 0, params%g_zcap))
     $        call errquit (pname//"zcap create failed", 0, GA_ERR)

         call rt_tddft_spatialcap (params, nbf_ao, params%g_zcap) !calc and store CAP
      endif
      

      end subroutine rt_tddft_init_common
c $Id$
