      subroutine et_fock(rtdb,geom,basis,
     1  g_pa,g_pb,nbf,tol2e,omega2)
c
c Try using shell_fock_build to calculate omega2
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
c
      integer i,j
      integer rtdb             ! Runtime database
      integer geom             ! Geometry handle
      integer basis            ! AO basis set handle
      integer nbf              ! Number of AO basis functions
      integer g_pa             ! GA handle for alpha density
      integer g_pb             ! GA handle for beta  density
      double precision tol2e   ! 2-e integral tolerance
      integer g_fock(4)        ! GA handles for (anti-)symmetrized fock matrix
      integer g_dens(4)        ! GA handles for (anti-)symmetrized density matrix
      double precision jfac(4) ! Prefactor for HF Coulomb for fock_2e
      double precision kfac(4) ! Prefactor for HF exchange for fock_2e
      double precision omega2, pa, pb, fJa, fJb, fKa, fKb
      integer itype
      character*3 theory
c
c ----------------------------------
c Initialize fock builder parameters
c ----------------------------------
c
      call int_init(rtdb,1,basis)
      call int_acc_std()
      call schwarz_init(geom,basis)
      call scf_get_fock_param(rtdb,tol2e)
      call fock_force_direct(rtdb)
c
c --------------------------------------------------------------------
c Allocate work arrays, set bgj theory (fock_2e won't work without it)
c --------------------------------------------------------------------
c
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory)) then
        call errquit('et_fock: unable to get theory',0)
      endif
c
      if (theory.eq.'scf') then
        itype=1
      else if (theory.eq.'hf') then
        itype=1
      else if (theory.eq.'dft') then
        itype=2
      else 
         call errquit('et_fock: theory must be scf or dft',1)
      endif
c 
      if (.not. rtdb_put(rtdb,'bgj:scf_type',mt_int,1, itype)) then
         call errquit('couldnt write bgj scftype',1)
      endif
c
      do i = 1,4
        if(.not.ga_create(mt_dbl,nbf,nbf,'fock',0,0,g_fock(i)))
     $  call errquit('et_fock: ga_create failed', 0)
        call ga_zero(g_fock(i))
c
        if(.not.ga_create(mt_dbl,nbf,nbf,'density',0,0,g_dens(i)))
     $  call errquit('et_fock: ga_create failed', 0)
        call ga_zero(g_dens(i))
      enddo
c
c SYMMETRIC CONTRIBUTION:
c
c -----------------
c Set jfac and kfac
c -----------------
c
      jfac(1) = 1.0d0
      jfac(2) = 1.0d0
      jfac(3) = 0.0d0
      jfac(4) = 0.0d0
c
      kfac(1) = 0.0d0
      kfac(2) = 0.0d0
      kfac(3) = 1.0d0
      kfac(4) = 1.0d0
c
c -----------------------------
c Copy Pa and Pb and symmetrize
c -----------------------------
c
      call ga_copy(g_pa,g_dens(1))
      call ga_copy(g_pb,g_dens(2))
      call ga_copy(g_pa,g_dens(3))
      call ga_copy(g_pb,g_dens(4))
c
      do i = 1,4
        call ga_symmetrize (g_dens(i))
      enddo
c
c ---------------------------------------------
c Call fock_2e to get symmetrized fock matrices
c ---------------------------------------------
c
      call fock_2e( geom, basis, 4, jfac, kfac,
     $     tol2e, .false., g_dens, g_fock, .false.)
c
c --------------------------------------------
c Calculate symmetrized contribution to omega2
c --------------------------------------------
c
      omega2 = 0.0d0
c
      do i = 1,nbf
      do j = 1,nbf
        call ga_get(g_fock(1), i, i, j, j,fJa,1)
        call ga_get(g_fock(2), i, i, j, j,fJb,1)
        call ga_get(g_fock(3), i, i, j, j,fKa,1)
        call ga_get(g_fock(4), i, i, j, j,fKb,1)
        call ga_get(g_dens(1), i, i, j, j,pa,1)
        call ga_get(g_dens(2), i, i, j, j,pb,1)
c
        omega2 = omega2 + pa*fJa
        omega2 = omega2 + pb*fJa
        omega2 = omega2 + pb*fJb
        omega2 = omega2 + pa*fJb
        omega2 = omega2 - pa*fKa
        omega2 = omega2 - pb*fKb
      enddo
      enddo
c
c
c
c ANTI-SYMMETRIC CONTRIBUTION:
c
c -----------------
c Set jfac and kfac
c -----------------
c
      jfac(1) = 0.0d0
      jfac(2) = 0.0d0
c
      kfac(1) = 1.0d0
      kfac(2) = 1.0d0
c
c ----------------------------------
c Copy Pa and Pb and ANTI-symmetrize
c ---------------------------------
c
      call ga_copy(g_pa,g_dens(1))
      call ga_copy(g_pb,g_dens(2))
c
      do i = 1,2
        call ga_antisymmetrize (g_dens(i))
      enddo
c
c --------------------------------------------------
c Call fock_2e then get ANTI-symmetrized fock matrices
c --------------------------------------------------
c
      call fock_2e( geom, basis, 2, jfac, kfac,
     $     tol2e, .false., g_dens, g_fock, .true.)
c
      do i = 1,2
        call ga_antisymmetrize (g_fock(i))
      enddo
c
c --------------------------------------------
c Calculate ANTI-symmetrized contribution to omega2
c --------------------------------------------
c
      do i = 1,nbf
      do j = 1,nbf
        call ga_get(g_fock(1), i, i, j, j,fKa,1)
        call ga_get(g_fock(2), i, i, j, j,fKb,1)
        call ga_get(g_dens(1), i, i, j, j,pa,1)
        call ga_get(g_dens(2), i, i, j, j,pb,1)
c
        omega2 = omega2 - pa*fKa 
        omega2 = omega2 - pb*fKb 
      enddo
      enddo
c
c
c --------------
c Destroy arrays
c --------------
c
      do i = 1,4
        if(.not.ga_destroy(g_fock(i)))
     $  call errquit('et_fock: ga_destroy failed', 0)
c
        if(.not.ga_destroy(g_dens(i)))
     $  call errquit('et_fock: ga_create failed', 0)
      enddo
c
       call fock_2e_tidy(rtdb)
       call schwarz_tidy()
c
       return
       end


