C
C     zfock_os_build.F
C
C     Builds open shell complex Fock matrix from a complex density
C     matrix.
C
C     A comment on notation.  The "z" in g_zfock means this GA is
C     complex data type.  The corresponding real and imaginary parts are
C     called g_fockre and g_fockim, which are both real data types, so
C     no "z".
C
C     Note the overall structure is slightly different from the closed
C     shell case.  Here, we cant lump coul and exch calculations into
C     one (for speed like we did in closed shell), so there are less
C     "combinations" of Fock build.
C     
C
C     All quantities are in the atomic orbital (AO) basis.
C
      
      subroutine zfock_os_build (params, g_zdens, energies, g_zfock)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "case.fh"
#include "zora.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zdens(2)    !alpha, beta spins; in AO bas

      
C     == Outputs ==
      type(rt_energies_t), intent(out) :: energies
      integer, intent(in)              :: g_zfock(2)   !alpha, beta spins; in AO bas


C     == Parameters ==
      character(*), parameter :: pname = "zfock_os_build: "
      integer, parameter      :: nspin = 2

      
C     == Variables ==
      logical need_exactexch
      logical need_dftxc
      integer g_densre(2)          !real part of dens mat
      integer g_densim(2)          !imag part of dens mat
      integer g_fockre(2)          !real part of fock mat
      integer g_fockim(2)          !imag part of fock mat
      integer g_v1(2), g_v2(2)     !potentials--results of each piece of F
      double precision ener1, ener2 !energies-results of each piece of F
      double precision dftExc(2)
      logical status
      double precision elapsed, elapsed_tot
      integer is  !spin counter
      integer g_zcapmo, g_zfock_mo  ! tmp GAs for MO CAP


C     == External routines ==
      logical xc_gotxc
      external xc_gotxc
      external sandbox_complex


C      call rt_tddft_os_confirm (params,'junk')

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
         call ga_copy (params%g_zfock_ao_init(1), g_zfock(1))
         call ga_copy (params%g_zfock_ao_init(2), g_zfock(2))
         call rt_tddft_print_notice ("Static OS Fock matrix")
         return
      endif



      call ga_zero (g_zfock(1))
      call ga_zero (g_zfock(2))

      need_exactexch = abs(xfac(1)).gt.1d-8
      need_dftxc = xc_gotxc()


C
C     Allocate GAs
C
      if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao,
     $     "Re[Dens](1)", 0, 0, g_densre(1)))
     $     call errquit ("couldnt create g_densre(1)", 0, GA_ERR)

      if (.not. ga_duplicate(g_densre(1), g_densre(2), "Re[Dens](2)"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)

      do is = 1, nspin
         if (.not. ga_duplicate(g_densre(1), g_densim(is), "Im[Dens]"))
     $        call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
         if (.not. ga_duplicate(g_densre(1), g_fockre(is), "Re[F]"))
     $        call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
         if (.not. ga_duplicate(g_densre(1), g_fockim(is), "Im[F]"))
     $        call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
         if (.not. ga_duplicate(g_densre(1), g_v1(is), "V1"))
     $        call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
         if (.not. ga_duplicate(g_densre(1), g_v2(is), "V2"))
     $        call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      enddo

      

C     == Extract real and imag parts of density matrix ==
      if (params%prof) call prof_start (elapsed)
      do is = 1, nspin
         call ga_zero (g_densre(is))
         call ga_zero (g_densim(is))
         call convert_z2d (g_zdens(is), g_densre(is), g_densim(is))
      enddo
      if (params%prof) call prof_end (elapsed, "Fock OS z2d")


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
      

      do is = 1, nspin
         call ga_zero (g_fockre(is))
         call ga_zero (g_fockim(is))
      enddo

C
C     == Standard core ==
C
      call zfock_os_core (params, g_densre, g_v1)
      call ga_add (1d0, g_v1(1), 1d0, g_fockre(1), g_fockre(1))
      call ga_add (1d0, g_v1(2), 1d0, g_fockre(2), g_fockre(2))
      energies%core = ga_ddot (g_densre(1), g_v1(1)) +
     $     ga_ddot (g_densre(2), g_v1(2))

      
C     == Scalar ZORA core ==
      if (params%lzorasf) then
         call ga_add (1d0, params%g_zora_sf(1),
     $        1d0, g_fockre(1), g_fockre(1))
         call ga_add (1d0, params%g_zora_sf(2),
     $        1d0, g_fockre(2), g_fockre(2))

C         call zfock_os_core_scalarzora (params, g_densre, g_v1)
C         call ga_add (1d0, g_v1(1), 1d0, g_fockre(1), g_fockre(1))
C         call ga_add (1d0, g_v1(2), 1d0, g_fockre(2), g_fockre(2))
C         energies%core = ga_ddot (g_densre(1), g_v1(1)) +
C     $        ga_ddot (g_densre(2), g_v1(2))
      endif


      energies%core = ga_ddot (g_densre(1), g_fockre(1)) +
     $     ga_ddot (g_densre(2), g_fockre(2))

C
C     == Coulomb part ==
C
C     Depends only on the real part of the *total* alpha+beta density matrix.
C
      if (cdfit) then
         call zfock_os_coulcdfit (params, g_densre, ener1, g_v1)
      else
         call zfock_os_coul (params, g_densre, ener1, g_v1)
      endif
C      write (luout, *) "ZEROED COUL"
      call ga_add (1d0, g_v1(1), 1d0, g_fockre(1), g_fockre(1))
      call ga_add (1d0, g_v1(2), 1d0, g_fockre(2), g_fockre(2))
      energies%coul = ener1


C
C     == Exact exchange ==
C
C     cam_exch is .true. if we are doing CAM, .false. otherwise.  Real
C     part will be returned in v1; imag part in v2.
C
C     Depends on the real and imaginary parts of the alpha,beta dens
C     mats.
C

      if (need_exactexch) then
         call zfock_os_exch (params, g_densre, g_densim,
     $        cam_exch, ener1, g_v1, g_v2)

         call ga_add (1d0, g_v1(1), 1d0, g_fockre(1), g_fockre(1))
         call ga_add (1d0, g_v1(2), 1d0, g_fockre(2), g_fockre(2))
         call ga_add (1d0, g_v2(1), 1d0, g_fockim(1), g_fockim(1))
         call ga_add (1d0, g_v2(2), 1d0, g_fockim(2), g_fockim(2))
         energies%xc(1) = ener1
         energies%xc(2) = 0d0
      endif

C
C     == DFT exchange-correlation ==
C
C     Depends only on the real part of the density matrix.
C      
      if (need_dftxc) then
         dftExc(1) = 0d0
         dftExc(2) = 0d0
         call zfock_os_dftxc (params, g_densre, dftExc, g_v1)
         call ga_dadd (1d0, g_v1(1), 1d0, g_fockre(1), g_fockre(1))
         call ga_dadd (1d0, g_v1(2), 1d0, g_fockre(2), g_fockre(2))
         energies%xc(1) = energies%xc(1) + dftExc(1)
         energies%xc(2) = energies%xc(2) + dftExc(2)
      endif

      
C      
C     Pack real and imag parts of Fock matrix back into g_zfock.  No
C     need to repack density matrix, as its value should not have
C     changed. Perhaps double check this??
C
      if (params%prof) call prof_start (elapsed)
      call convert_d2z (1d0, g_fockre(1), 1d0, g_fockim(1), g_zfock(1))
      call convert_d2z (1d0, g_fockre(2), 1d0, g_fockim(2), g_zfock(2))
      if (params%prof) call prof_end (elapsed, "Fock OS d2z")



C
C     Apply molecular orbital complex absorbing potential (MO CAP).
C     Potential is different for alpha and beta spins, as it is
C     built from their respective Fock matrices.
C
      if (params%mocap_active) then
         if (.not. ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "Fock MO tmp" , 0, 0, g_zfock_mo))
     $        call errquit(pname//"GA alloc failed",0, GA_ERR)
         
         if (.not. ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "Fock MO tmp" , 0, 0, g_zcapmo))
     $        call errquit(pname//"GA alloc failed",0, GA_ERR)

C     (alpha spin)
         call canorg_trans (params,"F","AO->MO",
     $        g_zfock(1), g_zfock_mo)
         call rt_tddft_mocap (params, g_zfock_mo, g_zcapmo)
         call ga_zadd (z1, g_zfock_mo,
     $        z1, g_zcapmo, g_zfock_mo) 
         call canorg_trans (params,"F","MO->AO",
     $        g_zfock_mo, g_zfock(1))
         call rt_tddft_print_notice ("Applied MO CAP to alpha spin")

C     (beta spin)
         call canorg_trans (params,"F","AO->MO",
     $        g_zfock(2), g_zfock_mo)
         call rt_tddft_mocap (params, g_zfock_mo, g_zcapmo)
         call ga_zadd (z1, g_zfock_mo,
     $        z1, g_zcapmo, g_zfock_mo) 
         call canorg_trans (params,"F","MO->AO",
     $        g_zfock_mo, g_zfock(2))
         call rt_tddft_print_notice ("Applied MO CAP to beta spin")

         
         if (.not. ga_destroy(g_zfock_mo))
     $        call errquit(pname//"GA destroy failed",0, GA_ERR)
         if (.not. ga_destroy(g_zcapmo))
     $        call errquit(pname//"GA destroy failed",0, GA_ERR)
      endif


C     == Spatial complex absorbing boundary conditions (CAP) ==
      if (params%cap_active) then
         call ga_zadd (z1, g_zfock(1), z1, params%g_zcap, g_zfock(1))
         call rt_tddft_print_notice ("Applied spatial CAP to "//
     $        "alpha Fock matrix")
         call ga_zadd (z1, g_zfock(2), z1, params%g_zcap, g_zfock(2))
         call rt_tddft_print_notice ("Applied spatial CAP to "//
     $        "beta Fock matrix")
      endif
      

C
C     Extra checks, disabled for speed.
C
      if (params%checklvl .ge. 2) then
         if (.not. mat_is_hermitian(g_zfock(1), params%tol_zero))
     $        call errquit ("F(1) not hermitian in fock builder", 0, 0)
         if (.not. mat_is_hermitian(g_zfock(2), params%tol_zero))
     $        call errquit ("F(2) not hermitian in fock builder", 0, 0)
         
         if (.not. mat_is_symmetric (g_fockim(1),"A",params%tol_zero))
     $        call errquit ("Im[F(1)] not antisym in fock build", 0, 0)
         if (.not. mat_is_symmetric (g_fockim(2),"A",params%tol_zero))
     $        call errquit ("Im[F(2)] not antisym in fock build", 0, 0)

         if (energies%xc(2) > 1d-8)
     $        call errquit (pname//"Exc(2) /= 0?", 0, 0)
      endif


C
C     Clean up.
C
      status = .true.
      do is = 1, nspin
         status=status.and.ga_destroy(g_densre(is))
         status=status.and.ga_destroy(g_densim(is))
         status=status.and.ga_destroy(g_fockre(is))
         status=status.and.ga_destroy(g_fockim(is))
         status=status.and.ga_destroy(g_v1(is))
         status=status.and.ga_destroy(g_v2(is))
      enddo
      
      if (.not.status)
     $     call errquit (pname//"couldnt free arrays", 0, GA_ERR)

      if (params%prof) call prof_end(elapsed_tot, "Fock OS total build")

      end subroutine

c $Id$
