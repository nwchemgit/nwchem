      subroutine et_fock(rtdb,geom,basis,
     1  g_pa,g_pb,nbf,tol2e,omega2)
c
c This routine uses fock_2e to calculate the 2e contribution to Vab
c
      implicit none
#include "errquit.fh"
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
      integer i
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
      double precision omega2
      integer itype
      character*3 theory
c
 9993 format(1x,'Two-electron integral screening (tol2e) : ',1Pd8.2)
      if(ga_nodeid().eq.0) then
       write(LuOut,*) 
       write(LuOut,*) 'Beginning calculation of 2e contribution'
       write(LuOut,9993) tol2e
       write(LuOut,*)
      endif

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
        call errquit('et_fock: unable to get theory',0, RTDB_ERR)
      endif
c
      if (theory.eq.'scf') then
        itype=1
      else if (theory.eq.'hf') then
        itype=1
      else if (theory.eq.'dft') then
        itype=2
      else 
         call errquit('et_fock: theory must be scf or dft',1, RTDB_ERR)
      endif
c 
      if (.not. rtdb_put(rtdb,'bgj:scf_type',mt_int,1, itype)) then
         call errquit('couldnt write bgj scftype',1, RTDB_ERR)
      endif
c
      do i = 1,4
        if(.not.ga_create(mt_dbl,nbf,nbf,'fock',0,0,g_fock(i)))
     $  call errquit('et_fock: ga_create failed', 0, GA_ERR)
        call ga_zero(g_fock(i))
c
        if(.not.ga_create(mt_dbl,nbf,nbf,'density',0,0,g_dens(i)))
     $  call errquit('et_fock: ga_create failed', 0, GA_ERR)
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
      omega2 = omega2 + ga_ddot (g_fock(1),g_dens(1)) 
      omega2 = omega2 + ga_ddot (g_fock(1),g_dens(2)) 
      omega2 = omega2 + ga_ddot (g_fock(2),g_dens(1)) 
      omega2 = omega2 + ga_ddot (g_fock(2),g_dens(2)) 
      omega2 = omega2 - ga_ddot (g_fock(3),g_dens(1)) 
      omega2 = omega2 - ga_ddot (g_fock(4),g_dens(2)) 
c
c      do i = 1,nbf
c      do j = 1,nbf
c        call ga_get(g_fock(1), i, i, j, j,fJa,1)
c        call ga_get(g_fock(2), i, i, j, j,fJb,1)
c        call ga_get(g_fock(3), i, i, j, j,fKa,1)
c        call ga_get(g_fock(4), i, i, j, j,fKb,1)
c        call ga_get(g_dens(1), i, i, j, j,pa,1)
c        call ga_get(g_dens(2), i, i, j, j,pb,1)
c
c        omega2 = omega2 + pa*fJa
c        omega2 = omega2 + pb*fJa
c        omega2 = omega2 + pb*fJb
c        omega2 = omega2 + pa*fJb
c        omega2 = omega2 - pa*fKa
c        omega2 = omega2 - pb*fKb
c      enddo
c      enddo
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
      omega2 = omega2 - ga_ddot (g_fock(1),g_dens(1)) 
      omega2 = omega2 - ga_ddot (g_fock(2),g_dens(2)) 
c
c      do i = 1,nbf
c      do j = 1,nbf
c        call ga_get(g_fock(1), i, i, j, j,fKa,1)
c        call ga_get(g_fock(2), i, i, j, j,fKb,1)
c        call ga_get(g_dens(1), i, i, j, j,pa,1)
c        call ga_get(g_dens(2), i, i, j, j,pb,1)
c
c        omega2 = omega2 - pa*fKa 
c        omega2 = omega2 - pb*fKb 
c      enddo
c      enddo
c
c    
c
c --------------
c Destroy arrays
c --------------
c
      do i = 1,4
        if(.not.ga_destroy(g_fock(i)))
     $  call errquit('et_fock: ga_destroy failed', 0, GA_ERR)
c
        if(.not.ga_destroy(g_dens(i)))
     $  call errquit('et_fock: ga_create failed', 0, GA_ERR)
      enddo
c
       call fock_2e_tidy(rtdb)
       call schwarz_tidy()
c
       return
       end


c $Id$
