C
C     zfock_cs_build.F
C
C     Builds closed shell complex Fock matrix from a complex density
C     matrix.
C
C     A comment on notation.  The "z" in g_zfock means this GA is
C     complex data type.  The corresponding real and imaginary parts are
C     called g_fockre and g_fockim, which are both real data types, so
C     no "z".
C
C     All quantities are in the atomic orbital (AO) basis.
C
      
      subroutine zfock_cs_build (params, g_zdens, energies, g_zfock)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "zora.fh"
#include "case.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zdens

      
C     == Outputs ==
      type(rt_energies_t), intent(out) :: energies
      integer, intent(in)              :: g_zfock


C     == Parameters ==
      character(*), parameter :: pname = "zfock_cs_build: "

      
C     == Variables ==
      logical need_exactexch
      logical need_dftxc
      integer g_densre          !real part of dens mat
      integer g_densim          !imag part of dens mat
      integer g_fockre          !real part of fock mat
      integer g_fockim          !imag part of fock mat
      integer g_v1, g_v2        !potentials--results of each piece of F
      double precision ener1, ener2 !energies-results of each piece of F
      double precision dftExc(2)
      logical status
      double precision elapsed, elapsed_tot
      integer g_zfock_mo, g_zcapmo ! for MO damping
      
C     == External routines ==
      logical xc_gotxc
      external xc_gotxc
      external sandbox_complex


      call rt_tddft_cs_confirm (params,'zfock_cs_build.F')

      
C     == Initializations ==
      if (params%prof) call prof_start (elapsed_tot)


      if (params%nExc .ne. 1)
     $     call errquit (pname//"might not work for nExc /= 1",0,0)

      if (oskel)
     $     call errquit (pname//"not working with symmetry", 0, 0)

      if (nbf_ao .ne. params%nbf_ao) call errquit (pname//
     $     "nbf_ao in params /= nbf_ao in cdft header", 0, 0)


C
C     If asked for (usually a "static" calculation), dont build but
C     instead use stored Fock matrix.
C
      if (params%skip_fock) then
         call ga_copy (params%g_zfock_ao_init(1), g_zfock)
         call rt_tddft_print_notice ("Static CS Fock matrix")
         return
      endif

      
      call ga_zero (g_zfock)
      need_exactexch = abs(xfac(1)).gt.1d-8
      need_dftxc = xc_gotxc()


      if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao,
     $     "Re[Dens]", 0, 0, g_densre))
     $     call errquit ("couldnt create g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_densim, "Im[Dens]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_fockre, "Re[F]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_fockim, "Im[F]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_v1, "V1"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_v2, "V2"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)

      
C     == Extract real and imag parts of density matrix ==
      if (params%prof) call prof_start (elapsed)
      call ga_zero (g_densre)
      call ga_zero (g_densim)
      call convert_z2d (g_zdens, g_densre, g_densim)
      if (params%prof) call prof_end (elapsed, "Fock CS z2d")


C
C     == Compute complex Fock matrix ==
C
C     For each piece we compute the energy and potential, then
C     accumulate the result in the real or imag part of F.  Note that
C     the only piece that involves the imag part of the density matrix
C     is the exact exchange.  We also only call the DFT XC routine if we
C     need it, i.e., anything but pure Hartree-Fock.
C
      energies%core = 0d0
      energies%coul = 0d0
      energies%xc(1) = 0d0
      energies%xc(2) = 0d0

      
      call ga_zero (g_fockre)
      call ga_zero (g_fockim)


C     == Standard core (kinetic+potential) ==
      call zfock_cs_core (params, g_densre, g_v1)
      call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
C      energies%core = ga_ddot (g_densre, g_v1)


C     == Scalar ZORA core ==
      if (params%lzorasf) then
C         call zfock_cs_core_scalarzora (params, g_densre, g_v1)
         call ga_add (1d0, params%g_zora_sf(1), 1d0, g_fockre, g_fockre)
C         energies%core = energies%core + ga_ddot (g_densre, g_v1)
      endif

      energies%core = ga_ddot (g_densre, g_fockre)

      
C     == Standard coul + standard exch ==
      if ((need_exactexch).and.(.not.cam_exch).and.(.not.cdfit)) then
         call zfock_cs_coul_exchre(params, g_densre, ener1, ener2, g_v1)
         call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
         energies%coul = ener1
         energies%xc(1) = ener2
         energies%xc(2) = 0d0

         call zfock_cs_exchim (params, g_densim, ener1, g_v1)
         call ga_add (1d0, g_v1, 1d0, g_fockim, g_fockim)
         energies%xc(1) = energies%xc(1) + ener1
         

C     == Standard coul + CAM exchange ==
      elseif ((need_exactexch).and.(cam_exch).and.(.not.cdfit)) then
         call zfock_cs_coul (params, g_densre, ener1, g_v1)
         call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
         energies%coul = ener1

         call zfock_cs_exch (params, g_densre, g_densim, .true.,
     $        ener1, g_v1, g_v2)
         call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
         call ga_add (1d0, g_v2, 1d0, g_fockim, g_fockim)
         energies%xc(1) = ener1
         energies%xc(2) = 0d0


C     == CD fit coul + standard exch ==
      elseif ((need_exactexch).and.(.not. cam_exch).and.(cdfit)) then
         call zfock_cs_coulcdfit (params, g_densre, ener1, g_v1)
         call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
         energies%coul = ener1

         call zfock_cs_exch (params, g_densre, g_densim, .false.,
     $        ener1, g_v1, g_v2)
         call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
         call ga_add (1d0, g_v2, 1d0, g_fockim, g_fockim)
         energies%xc(1) = ener1
         

C     == CD fit coul + CAM exch ==
      elseif ((need_exactexch).and.(cam_exch).and.(cdfit)) then
         call zfock_cs_coulcdfit (params, g_densre, ener1, g_v1)
         call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
         energies%coul = ener1

         call zfock_cs_exch (params, g_densre, g_densim, .true.,
     $        ener1, g_v1, g_v2)
         call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
         call ga_add (1d0, g_v2, 1d0, g_fockim, g_fockim)
         energies%xc(1) = ener1
         energies%xc(2) = 0d0


C     == Standard coul + no exchange ==
      elseif ((.not.need_exactexch).and.(.not.cdfit)) then
         energies%xc(1) = 0d0
         energies%xc(2) = 0d0

         call zfock_cs_coul (params, g_densre, ener1, g_v1)
         call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
         energies%coul = ener1


C     == CD fit coul + no exchange ==
      elseif ((.not.need_exactexch).and.(cdfit)) then
         energies%xc(1) = 0d0
         energies%xc(2) = 0d0

         call zfock_cs_coulcdfit (params, g_densre, ener1, g_v1)
         call ga_add (1d0, g_v1, 1d0, g_fockre, g_fockre)
         energies%coul = ener1
         
      else
         call errquit (pname//"Unsupported combination of flags", 0, 0)
      endif


C     == DFT XC ==
      if (need_dftxc) then
         dftExc(1) = 0d0
         dftExc(2) = 0d0
         call zfock_cs_dftxc (params, g_densre, dftExc, g_v1)
         call ga_dadd (1d0, g_v1, 1d0, g_fockre, g_fockre)
         energies%xc(1) = energies%xc(1) + dftExc(1)
         energies%xc(2) = energies%xc(2) + dftExc(2)
      endif


C
C     Fudge symmetry
C
      if (oskel) call errquit (pname//"not working with oskel", 0, 0)

      

C     == Pack real and imag parts of Fock matrix back into g_zfock ==
C     No need to repack density matrix, as its value should
C     not have changed. Perhaps double check this??

      if (params%prof) call prof_start (elapsed)
      call convert_d2z (1d0, g_fockre, 1d0, g_fockim, g_zfock)
      if (params%prof) call prof_end (elapsed, "Fock CS d2z")



C
C     Apply molecular orbital complex absorbing potential (MO CAP)
C
      if (params%mocap_active) then
         if (.not. ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "Fock MO tmp" , 0, 0, g_zfock_mo))
     $        call errquit(pname//"GA alloc failed",0, GA_ERR)
         
         if (.not. ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "Fock MO tmp" , 0, 0, g_zcapmo))
     $        call errquit(pname//"GA alloc failed",0, GA_ERR)

         call canorg_trans (params,"F","AO->MO", g_zfock, g_zfock_mo)
         call rt_tddft_mocap (params, g_zfock_mo, g_zcapmo)
         call ga_zadd (z1, g_zfock_mo, z1, g_zcapmo, g_zfock_mo) 
         call canorg_trans (params,"F","MO->AO", g_zfock_mo, g_zfock)
         call rt_tddft_print_notice ("Applied closed shell MO CAP")
         
         if (.not. ga_destroy(g_zfock_mo))
     $        call errquit(pname//"GA destroy failed",0, GA_ERR)

         if (.not. ga_destroy(g_zcapmo))
     $        call errquit(pname//"GA destroy failed",0, GA_ERR)
      endif


C     == Spatial complex absorbing boundary conditions (CAP) ==
      if (params%cap_active) then
         call ga_zadd (z1, g_zfock, z1, params%g_zcap, g_zfock)
         call rt_tddft_print_notice ("Applied spatial CAP to "//
     $        "restricted Fock matrix")
      endif

      
C
C     Extra checks, disabled for speed.
C
      if (params%checklvl .ge. 2) then
         if (.not. mat_is_hermitian (g_zfock, params%tol_zero))
     $        call errquit ("F not hermitian in fock builder", 0, 0)
         if (.not. mat_is_symmetric (g_fockim, "A", params%tol_zero))
     $        call errquit ("Im[F] not antisym in fock builder", 0, 0)
         if (energies%xc(2) > 1d-8)
     $        call errquit (pname//"Exc(2) /= 0?", 0, 0)
      endif

      status = .true.
      status=status.and.ga_destroy(g_densre)
      status=status.and.ga_destroy(g_densim)
      status=status.and.ga_destroy(g_fockre)
      status=status.and.ga_destroy(g_fockim)
      status=status.and.ga_destroy(g_v1)
      status=status.and.ga_destroy(g_v2)
      
      if (.not.status)
     $     call errquit (pname//"couldnt free arrays", 0, GA_ERR)

      if (params%prof) call prof_end (elapsed_tot,
     $     "Fock CS total build")

      end subroutine

c $Id$
