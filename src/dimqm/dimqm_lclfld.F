      subroutine dimqm_lclfld(g_dipel, omega, lifetime, g_dipel_i)
c
c     
c    
c     Called from aoresponse_driver.F
c
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
#include "geom.fh"
#include "crohf.fh"
#include "cscf.fh"
c
c     Input Variables
      integer g_dipel        ! Global array handle to the dipole matrix
      logical lifetime        !damping or not
      integer g_dipel_i       !currently defined in dimqm.fh, might change
      double precision omega  ! freq value
c
c     Local variables
      integer g_tmp1, g_tmp2, g_dcv, g_dim_temp
      integer dims(3), chunk(3)
      integer alo(3), ahi(3)
      integer blo(2), bhi(2)
      integer clo(3), chi(3)
      integer xend
      integer icmplx
      
      double precision dx_r, dy_r, dz_r
      double precision dx_i, dy_i, dz_i
      integer l_dimxyz, k_dimxyz
      integer l_muind, k_muind
c      integer l_muold, k_muold
      double precision dsum
      external dsum
      integer i3, ivec, n
      integer l_fld, k_fld
      integer g_dim_r(2)
      integer g_dim_i(2)
      integer g_temp(2)
      integer nvir, voff, xoff
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      integer ipm
      double precision pre_factor
      character*(1) direction(3)
      character*(1) pm(2)
      logical seed_save
      data direction /'x', 'y', 'z'/
      data pm /'+', '-'/
c
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
      i3 = nDIM * 3
      icmplx = 1
      if(omega > ZERO) icmplx = 2

c      write(luout,*) "omega:", omega
c      write(luout,*) "icmplx:", icmplx

      
c      g_tmp1 = ga_create_atom_blocked(geom, basis, 'dim_lclfld: tmp1')
c      g_tmp2 = ga_create_atom_blocked(geom, basis, 'dim_lclfld: tmp2')
c     TODO: work out a way for FD without damping to run correctly, 
c           although, the case itself makes little sense
c      if(icmplx > 1 ) then
      if(lifetime) then
        alo(1) = nbf
        alo(2) = -1
        alo(3) = -1
        ahi(1) = nbf
        ahi(2) = nbf
        ahi(3) = 3
        if (.not.nga_create(MT_DBL,3,ahi,'e-dipole-i',alo,g_dipel_i))
     $    call errquit('lclfld: nga_create failed g_dipel_i',0,GA_ERR)
        call ga_zero(g_dipel_i)
      end if

      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf

      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:fld',l_fld,k_fld))
     $  call errquit('malloc dimrsp:fld failed',1,MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:muind',
     $                                            l_muind,k_muind))
     $  call errquit('malloc dimrsp:muind failed',1,MA_ERR)
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:muold',
     $                                            l_muold,k_muold))
     $  call errquit('malloc dimrsp:muind failed',1,MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,i3,'dimrsp:xyz',l_dimxyz,k_dimxyz))
     $  call errquit('malloc dimrsp:xyz failed',1,MA_ERR)
c
      if(.not.rtdb_get(dimqm_rtdb,'dimpar:coords',mt_dbl,i3,
     $                                dbl_mb(k_dimxyz)))
     $  call errquit('get dimpar:coords failed', 1, RTDB_ERR)

c     Loop over xyz dsrections
      seed_save = dimqm_noseed
      dimqm_noseed = .true.
      do n = 1, 3
c        do ipm = 1, icmplx ! +/- direction if it's complex
c
c   =============================
c   Solve for dipoles 
c   =============================
c
c   Zero arrays
        call dfill(i3*icmplx, ZERO, dbl_mb(k_muind),   1)
        call dfill(i3*icmplx, ZERO, dbl_mb(k_muold),   1)
        call dfill(i3*icmplx, ZERO, dbl_mb(k_fld),     1)
c
c   Ones for current cartesian direction
        call dfill(nDIM,      ONE,  dbl_mb(k_fld+n-1), 3)
c   
c   Calculate induced dipoles
        call dimqm_f2d(dimqm_rtdb, dbl_mb(k_fld),
     $                 dbl_mb(k_muind), dbl_mb(k_muold),
     $                 dbl_mb(k_dimxyz), icmplx, 's', ' ',.false.)
c       
        if(icmplx.eq.2 .and. lifetime) then
          if(.not.rtdb_put(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_r'//pm(1),
     $                           mt_dbl, i3, dbl_mb(k_muind)))
     $      call errquit('put dimqm:muind_r failed',1,RTDB_ERR)
          if(.not.rtdb_put(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_i'//pm(1),
     $                           mt_dbl, i3, dbl_mb(k_muind+i3)))
     $      call errquit('put dimqm:muind_i failed',1,RTDB_ERR)
        else
          if(.not.rtdb_put(dimqm_rtdb,
     $      'dimqm:muind_'//direction(n),
     $       mt_dbl, i3, dbl_mb(k_muind)))
     $      call errquit('put dimqm:muind_r failed',1,RTDB_ERR)
        end if
c        end do
      end do ! ivec = 1, 3
c
      if (.not. ma_chop_stack(l_fld)) call errquit('lclfld: MA?', 0,
     $       MA_ERR)
c
c   ====================================================
c   Solve for DIM potential
c   ====================================================
c
      dims(1) = 3
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
c
      if((icmplx.eq.2) .and. lifetime) then
c   Real +
        if (.not. nga_create (MT_DBL, 3, dims, 'lclfld:dim_r+',chunk,
     &                        g_dim_r(1)))
     &   call errquit('lclfld: could not allocate g_dim_r+',1,GA_ERR)
        call ga_zero(g_dim_r(1))
        call fock_dim(geom, nbf, basis, 3, g_dim_r(1), 1, 1)
        call ga_symmetrize(g_dim_r(1))
c!   Real -
 !       if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_r-',chunk,
 !    &                        g_dim_r(2)))
 !    &   call errquit('addop: could not allocate g_dim_r-',1,GA_ERR)
 !       call ga_zero(g_dim_r(2))
 !       call fock_dim(geom, nbf, basis, 3, g_dim_r(2), 2, 1)
 !       call ga_antisymmetrize(g_dim_r(2))
c   Imaginary +
        if (.not. nga_create (MT_DBL, 3, dims, 'lclfld:dim_i+',chunk,
     &                        g_dim_i(1)))
     &   call errquit('lclfld: could not allocate g_dim_i+',1,GA_ERR)
        call ga_zero(g_dim_i(1))
        call fock_dim(geom, nbf, basis, 3, g_dim_i(1), 1, 2)
        call ga_symmetrize(g_dim_i(1))
c!   Imaginary -
 !       if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_i-',chunk,
 !    &                        g_dim_i(2)))
 !    &   call errquit('addop: could not allocate g_dim_i-',1,GA_ERR)
 !       call ga_zero(g_dim_i(2))
 !       call fock_dim(geom, nbf, basis, 3, g_dim_i(2), 2, 2)
 !       call ga_antisymmetrize(g_dim_i(2))
 !
      else
        if (.not. nga_create (MT_DBL, 3, dims, 'lclfld:dim_r',chunk,
     &                        g_dim_r(1)))
     &   call errquit('lclfld: could not allocate g_dim_r',1,GA_ERR)
        call ga_zero(g_dim_r) 
        call fock_dim(geom, nbf, basis, 3, g_dim_r(1),  0, 1)
        call ga_symmetrize(g_dim_r(1))
c!   ======================================
c!   Undo the symmetrization to recover +/-
c!   ======================================
c!
 !       dims(1) = nbf
 !       dims(2) = nbf
 !       chunk(1) = dims(1)
 !       chunk(2) = -1
 !
 !       if (.not.nga_create(MT_DBL,2,dims,'gtemp',chunk,
 !    &     g_temp(1))) call
 !    &     errquit('dim_addop: nga_create failed gtemp',
 !    &     0,GA_ERR)
 !       call ga_zero(g_temp(1))
 !       if (.not.nga_create(MT_DBL,2,dims,'gtemp',chunk,
 !    &     g_temp(2))) call
 !    &     errquit('dim_addop: nga_create failed gtemp',
 !    &     0,GA_ERR)
 !       call ga_zero(g_temp(2))
 !
 !       do ivec = 1, 3
 !         alo(1) = ivec
 !         ahi(1) = ivec
c!       ************
c!       Real portion
c!       ************
 !       call nga_copy_patch ('N',g_dim_r(1),alo,ahi,g_temp(1),blo,bhi)
 !       call nga_copy_patch ('N',g_dim_r(2),alo,ahi,g_temp(2),blo,bhi)
!c
!c       it might be necessary to use 0.5 here instead of 1.0
!c       (note: that turned out NOT to be the case after some testing)
!          pre_factor = 1.0d0
 !         call ga_sync()
 !         call nga_add_patch (pre_factor, g_temp(1), blo, bhi,
 !    &       pre_factor, g_temp(2), blo, bhi,
 !    &       g_dim_r(1), alo, ahi)
 !         call nga_add_patch (pre_factor, g_temp(1), blo, bhi,
 !    &       -pre_factor, g_temp(2), blo, bhi,
 !    &       g_dim_r(2), alo, ahi)
c!       *****************
c!       Imaginary portion
c!       *****************
 !       call nga_copy_patch ('N',g_dim_i(1),alo,ahi,g_temp(1),blo,bhi)
 !       call nga_copy_patch ('N',g_dim_i(2),alo,ahi,g_temp(2),blo,bhi)
!c
c!       it might be necessary to use 0.5 here instead of 1.0
c!       (note: that turned out NOT to be the case after some testing)
 !       pre_factor = 1.0d0
 !         call ga_sync()
c!         real perturbation: 
 !         call nga_add_patch (pre_factor, g_temp(1), blo, bhi,
 !    &       pre_factor, g_temp(2), blo, bhi,
 !    &       g_dim_i(1), alo, ahi)
 !         call nga_add_patch (pre_factor, g_temp(1), blo, bhi,
 !    &       -pre_factor, g_temp(2), blo, bhi,
 !    &       g_dim_i(2), alo, ahi)
 !       enddo                     ! ivec = 1,nvec
      end if
c
c   ========================================
c   Add DIM local field to the dipole matrix
c   ========================================
c start: some debug stuff -------------------
c      call ga_print(g_dim_r(1))
cc      call ga_print(g_dim_r(2))
c      call ga_print(g_dipel)
c
c      if (icmplx > 1) then
c         call ga_print(g_dim_i(1))
cc         call ga_print(g_dim_i(2))
c      endif
c      write(luout,*)'end of g_dim print'
c end: some debug stuff --------------------
    
c
      xoff = 1
      xend = nvir*nclosed
c      call ga_print(g_movecs)
      call ga_sync()
      clo(1) = 1
      chi(1) = nbf
      clo(2) = 1
      chi(2) = nbf
!jbecca START: commenting this stuff out just to test
!      do ivec = 1, 3 ! Loop over perturbations
!        alo(1) = ivec
!        ahi(1) = ivec
!        clo(3) = ivec
!        chi(3) = ivec
!        call nga_copy_patch('N',g_dim_r(1),alo,ahi,g_dcv,blo,bhi)
!        call nga_add_patch(ONE, g_dipel, clo, chi,
!     $                     ONE, g_dcv, blo,  bhi,
!     $                          g_dipel, clo, chi)
!        if(icmplx > 1) then
!          call nga_copy_patch('N',g_dim_i(1),alo,ahi,g_dcv,blo,bhi)
!          call nga_add_patch(ONE, g_dipel_i, clo, chi,
!     $                       ONE, g_dcv, blo,  bhi,
!     $                            g_dipel_i, clo, chi)
! 
!c        call ga_zero(g_dipel_i)
!        end if
!      end do !ivec = 1, 3
!jbecca END
!jbecca START: It seems that ga_create_atom_blocked does not work 
!              for this case. Making this by hand now. 
      if (.not. ga_create(MT_DBL,nbf,nbf,'dim_temp_pert',chunk(1),
     $    chunk(2),g_dim_temp)) call errquit('lclfld: 
     $      nga_create failed g_dim_temp',0,GA_ERR)

      do ivec = 1, 3    ! loop over pert
         alo(1) = ivec
         ahi(1) = ivec
         clo(3) = ivec
         chi(3) = ivec
         call ga_zero(g_dim_temp)
         call nga_copy_patch('N',g_dim_r(1),alo,ahi,g_dim_temp,blo,bhi)
         call nga_add_patch(ONE, g_dipel, clo, chi,
     $                      ONE, g_dim_temp, blo,  bhi,
     $                          g_dipel, clo, chi)

c         if (icmplx > 1 ) then
         if (lifetime) then
            call ga_zero(g_dim_temp)
            call nga_copy_patch('N',g_dim_i(1),alo,ahi,g_dim_temp,
     $                           blo,bhi)
            call nga_add_patch(ONE, g_dipel_i, clo, chi,
     $                       ONE, g_dim_temp, blo,  bhi,
     $                            g_dipel_i, clo, chi)
         endif
      enddo       !ivec

      if (.not. ga_destroy(g_dim_temp)) 
     &   call errquit('lclfld: g_dim_temp destroy?',0,GA_ERR)
!jbecca END
c      call ga_print(g_dim_r(1))
c   ========
c   Clean up
c   =======
      if (.not. ga_destroy(g_dim_r(1))) 
     &      call errquit('lclfld: g_dim_r destroy?',0, GA_ERR)
      if(icmplx > 1 .and. lifetime) then
  !    if (.not. ga_destroy(g_dim_r(2))) call errquit('addop: GA?',0,
  !   &       GA_ERR)
        if (.not. ga_destroy(g_dim_i(1))) 
     &      call errquit('lclfld:g_dim_i destroy?',0,GA_ERR)
  !      if (.not. ga_destroy(g_dim_i(2))) call errquit('addop: GA?',0,
  !   &                                            GA_ERR)
  !      if (.not. ga_destroy(g_temp(1))) call errquit('addop: GA?',0,
  !   &       GA_ERR)
  !      if (.not. ga_destroy(g_temp(2))) call errquit('addop: GA?',0,
  !   &       GA_ERR)
      end if
!jbecca START
      dimqm_noseed = seed_save
c      call ga_print(g_dipel_i)
      end subroutine dimqm_lclfld



      subroutine dimqm_addlclfld(g_dip, omega)
c-----------------------------------------------------------------------
c  This subroutine is used to create and add in the local field
c  operator to the U (or A) matrix before solving the linear equations
c  in response routines. This ensures that the QM system's response is 
c  also dependent on any local fields that exist.
c
c  Called from: nothing right now
c
c  Calls:   subroutine dimqm_f2d
c
c  Author: Jeff Becca, jbb5516@psu.edu, 2017
c-----------------------------------------------------------------------  

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
#include "geom.fh"
#include "crohf.fh"
#include "cscf.fh"

c--------------------------- Input variables ---------------------------
      integer g_dip        !H10, used to initialize rhs
      integer omega        !frequency

c--------------------------- Local variables ---------------------------
      integer icmplx, i3, ivec, n
      integer alo(3), ahi(3), dims(3), chunk(3)
      integer l_fld, k_fld, l_muind, k_muind, l_dimxyz, k_dimxyz
      logical seed_save
      character*(1) direction(3)
      character*(1) pm(2)
      data direction /'x', 'y', 'z'/
      data pm /'+', '-'/ 
c---------------------------- Routine START ----------------------------
c  TODO: Update this for the case of complex local fields. Currently
c        only caring about the real part of the local field. 
      i3    =  nDIM * 3
      icmplx = 1
      if (omega > ZERO) icmplx = 2

c  get some DIM stuff needed for setting up field

      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:fld',l_fld,k_fld))
     $  call errquit('malloc dimrsp:fld failed',1,MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:muind',
     $                                            l_muind,k_muind))
     $  call errquit('malloc dimrsp:muind failed',1,MA_ERR)
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:muold',
     $                                            l_muold,k_muold))
     $  call errquit('malloc dimrsp:muind failed',1,MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,i3,'dimrsp:xyz',l_dimxyz,k_dimxyz))
     $  call errquit('malloc dimrsp:xyz failed',1,MA_ERR)
c
      if(.not.rtdb_get(dimqm_rtdb,'dimpar:coords',mt_dbl,i3,
     $                                dbl_mb(k_dimxyz)))
     $  call errquit('get dimpar:coords failed', 1, RTDB_ERR)

      seed_save = dimqm_noseed
      dimqm_noseed = .true.

      do n = 1, 3

c  zero arrays
         call dfill(i3*icmplx, ZERO, dbl_mb(k_muind), 1)
         call dfill(i3*icmplx, ZERO, dbl_mb(k_muold), 1)
         call dfill(i3*icmplx, ZERO, dbl_mb(k_fld),   1)

c  Create unit field in current cartesian direction for field
         call dfill(nDIM,  ONE,  dbl_mb(k_fld+n-1),   3)

c  Calculate the induced dipoles from the unit field
         call dimqm_f2d(dimqm_rtdb, dbl_mb(k_fld),
     $            dbl_mb(k_muind),  dbl_mb(k_muold),
     $            dbl_mb(k_dimxyz), icmplx, 's', ' ',.false.)

c  Store dipoles
         if (icmplx.eq.1) then
            if(.not.rtdb_put(dimqm_rtdb,
     $      'dimqm:muind_'//direction(n),
     $      mt_dbl, i3, dbl_mb(k_muind)))
     $      call errquit('put dimqm:muindr_r failed',1,RTDB_ERR)
         else
          if(.not.rtdb_put(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_r'//pm(1),
     $                           mt_dbl, i3, dbl_mb(k_muind)))
     $      call errquit('put dimqm:muind_r failed',1,RTDB_ERR)
          if(.not.rtdb_put(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_i'//pm(1),
     $                           mt_dbl, i3, dbl_mb(k_muind+i3)))
     $      call errquit('put dimqm:muind_i failed',1,RTDB_ERR)
        end if
      end do ! ivec = 1, 3

      if (.not. ma_chop_stack(l_fld)) call errquit('lclfld: MA?', 0,
     $       MA_ERR)
c
c   ====================================================
c   Solve for DIM potential
c   ====================================================
c
      dims(1) = 3
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      end subroutine dimqm_addlclfld
