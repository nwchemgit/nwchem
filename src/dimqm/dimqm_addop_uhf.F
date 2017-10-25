c ----------------------------------------------------------------------
c Subroutine used for calculating and adding in the DIM operator in 
c frequency dependent UHF calculations. This routine is highly experimental,
c and is still being worked on. Adapted from dimqm_addop. This routine
c works for the FD case with no damping. (Imaginary terms not needed)
c 
c The dimqm_addop_uhf_damp subroutine is a copy to have the lifetime 
c effects taken into account.
c
c Author: Jeff Becca, jbb5516@psu.edu, 2017
c ----------------------------------------------------------------------

      Subroutine dimqm_addop_uhf(g_x_r, g_x_i, g_Ax_r, g_Ax_i,
     $                      ncomp, limag, lifetime, g_dens_r, g_dens_i)

c     Called from: file - uhf_hessv2_ext.F
c                 subroutine - uhf_hessv_2e3
c
c     Subroutines called from: dimqm_EqmE.F, dimqm_f2d.F dim_fock_xc.F
c   
c     Calculates and adds the frequency-dependent DIM potential to the 
c     response Fock matricies (real and imaginary).  Requires knowledge
c     of both the real and imaginary vectors simultaneously which is
c     why this has to be located here, unlike the static routine.
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
      integer g_x_r(2)     ! A matrix handle (real)      [IN]
      integer g_x_i(2)     ! A matrix handle (imaginary) [IN]
      integer g_Ax_r(2)    ! F matrix handle (real)      [IN/OUT]
      integer g_Ax_i(2)    ! F matrix handle (imaginary) [IN/OUT]
      integer ncomp        ! num of components (+/-)     [IN]
      logical limag        ! Imaginary perturbation?     [IN]
      logical lifetime     ! Damping or no damping
      integer g_dens_r(2)  ! Perturbed pmat              [IN]
      integer g_dens_i(2)  ! perturbed pmat IMAG         [IN]
c
c     Local variables
      integer g_tmp1, g_tmp2, g_dcv
c      integer l_dimxyz, k_dimxyz
      double precision dimxyz(3, nDIM)
c      integer l_muind, k_muind
      double precision muind(3, nDIM, 2)
      integer dims(3), chunk(3)
      character*(255) cstemp
      integer g_pmats(2), g_pmata(2), g_h1mat(2)
      integer g_tmpwork(2)
      integer ipm
c      integer g_dens_r(2)
c      integer g_dens_i(2)
      integer alo(3), ahi(3)
      integer blo(2), bhi(2)
      integer g_dens_comp_r
      integer g_dens_comp_i
      integer xend
      double precision pre_factor
      double precision muold(3, nDIM, 2)
      
      double precision dx_r, dy_r, dz_r
      double precision dx_i, dy_i, dz_i
      double precision dsum, rmax
      external dsum
      integer i3, ivec, n
c      integer l_fld, k_fld
      double precision fld(3, nDIM, 2)
      integer g_dim_r(2)
      integer g_dim_i(2)
      integer nvir, voff, xoff
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      character*(1) direction(3)
      character*(1) pm(2)
      data direction /'x', 'y', 'z'/
      data pm /'+', '-'/
      logical firsttime
c      double precision screen(nDIM)
      double precision dimErr(3,2,2)
      double precision calcErr
      external calcErr
      integer id
c
      id = ga_nodeid()
      if (ldebug .and. id .eq. 0) then
        write(luout,*) "Start dimqm_addop_uhf"
      end if
      nvir = nmo - nclosed - nopen
      i3 = nDIM * 3
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'dim_addop: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'dim_addop: tmp2')

      dims(1) = nbf
      dims(2) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1

c
c   Allocate new arrays
c      if(.not.ma_push_get(mt_dbl,i3*2,'dimrsp:fld',l_fld,k_fld))
c     $  call errquit('malloc dimrsp:fld failed',1,MA_ERR)
c
c      if(.not.ma_push_get(mt_dbl,i3*2,'dimrsp:muind',
c     $                                            l_muind,k_muind))
c     $  call errquit('malloc dimrsp:muind failed',1,MA_ERR)
c
c      if(.not.ma_push_get(mt_dbl,i3,'dimrsp:xyz',l_dimxyz,k_dimxyz))
c     $  call errquit('malloc dimrsp:xyz failed',1,MA_ERR)
c
      if(.not.rtdb_get(dimqm_rtdb,'dimpar:coords', mt_dbl, i3, dimxyz))
     $  call errquit('get dimpar:coords failed', 1, RTDB_ERR)
c
      g_dens_comp_r = ga_create_atom_blocked(geom,basis,
     $                          'real density matrix comp')
      if (lifetime) then
        g_dens_comp_i = ga_create_atom_blocked(geom,basis,
     $                          'imag density matrix comp')
      end if
c
c
c      call dimqm_screening(dimqm_rtdb, geom, basis, dbl_mb(k_dimxyz),
c     $                       screen)
c      screen = ONE
c
c   =============================
c   Solve for induced dipoles +/-
c   =============================
c
c     Investigate the shape of g_dens_r
c      write(luout,*)'g_dens_r+'
c      call ga_print(g_dens_r(1))
c      write(luout,*)'g_dens_r-'
c      call ga_print(g_dens_r(2))
c     set dimension variables for patching
      alo(2)   =  1
      ahi(2)   =  nbf
      alo(3)   =  1
      ahi(3)   =  nbf
      blo(1)   =  1
      bhi(1)   =  nbf
      blo(2)   =  1
      bhi(2)   =  nbf
c     Loop over perturbations
      do n = 1, 3
        do ipm = 1, 2
          call ga_zero(g_dens_comp_r)
          if (lifetime) call ga_zero(g_dens_comp_i)
          alo(1) = n
          ahi(1) = n
c
c       Copy current perturbation into g_dens_comp
          call nga_copy_patch('N',g_dens_r(ipm), alo, ahi,
     $                          g_dens_comp_r, blo, bhi)
          if (lifetime) then
            call nga_copy_patch('N',g_dens_i(ipm), alo, ahi,
     $                              g_dens_comp_i, blo, bhi)
          end if
          muind = ZERO
          fld = ZERO
          firsttime = .false.
          if(.not.rtdb_get(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_r'//pm(ipm),
     $                              mt_dbl, i3, muold(:,:,1))) then
            if(id.eq.0 .and. ldebug) 
     $         write(luout,*) "First cycle, no old dipoles!"
            muold = ZERO
            firsttime = .true.
            dimqm_seeded = .false.
c            xyz_seeded(3*(n-1)+ipm) = .false.
            if(dimtol0 < 1.0d-4 .and. .not. dimqm_noseed) then
              dimtolxyz(ipm*3 - 1 + n) = 1.0d-4
              if(id.eq.0 .and. ldebug) then
                write(luout,*) "Requested tolerance below 1.0d-4"
                write(luout,*) "Setting "//direction(n)//pm(ipm)//
     $                         " dir tolerance to 1.0d-4 to start"
              end if
            end if
          else
            if(.not.rtdb_get(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_i'//pm(ipm),
     $                           mt_dbl, i3, muold(:,:,2)))
     $          call errquit('get dimqm:muold failed',1,RTDB_ERR)
          end if
c         Set convergence tolerance
c          dimtol = dimtolxyz(ipm*3 - 1 + n)
c          dimqm_seeded = xyz_seeded(ipm*3 - 1 + n)
c          dimtol = 1.0d-7
c          dimqm_noseed = .true.
c          call dfill(i3*2, ZERO, dbl_mb(k_muind), 1)
c          call dfill(i3*2, ZERO, dbl_mb(k_fld), 1)
c
c       Real portion of E-Field
c        write(luout,*) "REAL"
c        call ga_print(g_dens_comp_r)
          call dimqm_EqmE(dimqm_rtdb, g_dens_comp_r, geom, basis,
     $               fld(:,:,1), dimxyz)
c
c       Imaginary portion of E-Field
c        write(luout,*) "IMAG"
c        call ga_print(g_dens_comp_i)
          if (lifetime) then
            call dimqm_EqmE(dimqm_rtdb, g_dens_comp_i, geom, basis,
     $                      fld(:,:,2), dimxyz)
          end if
c
c       Solve for induced dipoles
          call dimqm_f2d(dimqm_rtdb, fld, muind, muold, dimxyz, 2,
     $                   direction(n), pm(ipm),.false.)
c
c         Write induced dipoles to RTDB
          dx_r = SUM(muind(1,:,1))
          dy_r = SUM(muind(2,:,1))
          dz_r = SUM(muind(3,:,1))
          dx_i = SUM(muind(1,:,2))
          dy_i = SUM(muind(2,:,2))
          dz_i = SUM(muind(3,:,2))
          if(id.eq.0.and.ldebug) then
            write(luout,*) "Total induced dipole moment for "//
     $                  direction(n)//pm(ipm)//" perturbation"
            write(luout,*) "X:", dx_r, dx_i
            write(luout,*) "Y:", dy_r, dy_i
            write(luout,*) "Z:", dz_r, dz_i
            write(luout,*) ''
          end if
          dimErr(n, ipm, 1) = calcErr(i3, muold(:,:,1), muind(:,:,1))
          dimErr(n, ipm, 2) = calcErr(i3, muold(:,:,2), muind(:,:,2))
          if(id.eq.0.and.ldebug) then
            write(luout,*) "Max error in real dipoles:",
     $                       dimErr(n, ipm, 1)
            write(luout,*) "Max error in imag dipoles:",
     $                       dimErr(n, ipm, 2)
          end if
c          if(dimErr(n, ipm, 1)/dimtol < HUNDRED
c     $              .and. dimErr(n, ipm, 2)/dimtol < HUNDRED
c     $              .and. .not. xyz_seeded(ipm*3 - 1 + n)
c     $              .and. .not. firsttime) then
c            xyz_seeded(ipm*3 - 1 + n) = .true.
c            write(luout,*) "Error within 10^2 of", dimtol, "for "//
c     $                     direction(n)//pm(ipm)//" dir"
c            write(luout,*) "Setting current "//direction(n)//pm(ipm)//
c     $                     " dir as seed"
c            write(luout,*)"Reverting tolerance back to", dimtol0
c            dimtolxyz(ipm*3 - 1 + n) = dimtol0
c          end if
          if(.not.rtdb_put(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_r'//pm(ipm),
     $                              mt_dbl, i3, muind(:,:,1)))
     $        call errquit('put dimqm:muind_p failed',1,RTDB_ERR)
          if(.not.rtdb_put(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_i'//pm(ipm),
     $                           mt_dbl, i3, muind(:,:,2)))
     $        call errquit('put dimqm:muind_p failed',1,RTDB_ERR)
        end do ! ipm = 1, 2
      end do ! ivec = 1, 3
c      if(MAXVAL(dimErr) <= 1.0d-4) then
c        write(luout,*) "Dipole error below 1d-4"
c        write(luout,*) "Shutting down DIM"
c        dimqm_on = .false.
c      end if
c
c   Destroy GAs we don't need anymore
      if (.not. ga_destroy(g_dens_comp_r)) call errquit
     $    ('addop: dens_comp_r GA?',0, GA_ERR)
      if (lifetime) then
      if (.not. ga_destroy(g_dens_comp_i)) call errquit
     $    ('addop: dens_comp_i GA?',0, GA_ERR)
      end if
c      do ipm = 1,2
c        if (.not. ga_destroy(g_dens_r(ipm))) call errquit
c     $     ('addop: dens_r GA?',0, GA_ERR)
c        if (lifetime) then
c         if (.not. ga_destroy(g_dens_i(ipm))) call errquit
c     $     ('addop: dens_i GA?',0, GA_ERR)
c        endif
c      end do
c
c   Deallocate l_fld, l_muind, l_dimxyz
c      if (.not. ma_chop_stack(l_fld)) call errquit
c     $   ('addop: fld MA?', 0, MA_ERR)
c
c   ====================================================
c   Solve for DIM potential, both real and imaginary S/A
c   ====================================================
c
      dims(1) = 3
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
c
c   Real +
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_r+',chunk,
     &                        g_dim_r(1)))
     &   call errquit('addop: could not allocate g_dim_r+',1,GA_ERR)
      call ga_zero(g_dim_r(1)) 
      call fock_dim(geom, nbf, basis, 3, g_dim_r(1), 1, 1)
      call ga_symmetrize(g_dim_r(1))
c
c   Real -
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_r-',chunk,
     &                        g_dim_r(2)))
     &   call errquit('addop: could not allocate g_dim_r-',1,GA_ERR)
      call ga_zero(g_dim_r(2))
      call fock_dim(geom, nbf, basis, 3, g_dim_r(2), 2, 1)
      call ga_antisymmetrize(g_dim_r(2))
      if (lifetime) then
c
c   Imaginary +
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_i+',chunk,
     &                        g_dim_i(1)))
     &   call errquit('addop: could not allocate g_dim_i+',1,GA_ERR)
      call ga_zero(g_dim_i(1))
      call fock_dim(geom, nbf, basis, 3, g_dim_i(1), 1, 2)
      call ga_symmetrize(g_dim_i(1))
c
c   Imaginary -
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_i-',chunk,
     &                        g_dim_i(2)))
     &   call errquit('addop: could not allocate g_dim_i-',1,GA_ERR)
      call ga_zero(g_dim_i(2))
      call fock_dim(geom, nbf, basis, 3, g_dim_i(2), 2, 2)
      call ga_antisymmetrize(g_dim_i(2))
      end if
c
c   ======================================
c   Undo the symmetrization to recover +/-
c   ======================================
      blo(1)   =  nbf
      blo(2)   =  nbf
      chunk(1) =  blo(1)
      chunk(2) =  -1

      do ipm   =  1, ncomp
         write(cstemp, '(a,i1)') 'g_tmpwork_',ipm
         if (.not.nga_create(MT_DBL,2,blo,cstemp(1:11),chunk,
     $         g_tmpwork(ipm))) call errquit('dim_addop_ufh: 
     $         nga_create failed '//cstemp(1:11),0,GA_ERR)
         call ga_zero(g_tmpwork(ipm))
      enddo
c     reset blo and bhi for future use
      blo(1)   =  1
      bhi(1)   =  nbf
      blo(2)   =  1
      bhi(2)   =  nbf
c
      do ivec = 1, 3
        alo(1) = ivec
        ahi(1) = ivec
c       ************
c       Real portion
c       ************
c TODO: I think the g_pmats here are being used just as temp arrays, 
c        but I really do not know for sure. 
c  all uses of g_pmats are being switched to g_tmpwork and the arrays
c  are allocated here.
        call nga_copy_patch('N',g_dim_r(1),alo,ahi,g_tmpwork(1),blo,bhi)
        call nga_copy_patch('N',g_dim_r(2),alo,ahi,g_tmpwork(2),blo,bhi)
      
c
c       it might be necessary to use 0.5 here instead of 1.0
c       (note: that turned out NOT to be the case after some testing)
        pre_factor = 1.0d0
        call ga_sync()
        if (.not.limag) then
c         real perturbation: 
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_r(1), alo, ahi)
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       -pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_r(2), alo, ahi)
        else
c         imaginary perturbation:
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_r(1), alo, ahi)
          call nga_add_patch (-pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_r(2), alo, ahi)
        end if  ! if .not.limag
        if (lifetime) then
c       *****************
c       Imaginary portion
c       *****************
        call nga_copy_patch('N',g_dim_i(1),alo,ahi,g_tmpwork(1),blo,bhi)
        call nga_copy_patch('N',g_dim_i(2),alo,ahi,g_tmpwork(2),blo,bhi)
c
c       it might be necessary to use 0.5 here instead of 1.0
c       (note: that turned out NOT to be the case after some testing)
        pre_factor = 1.0d0
        call ga_sync()
        if (.not.limag) then
c         real perturbation: 
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_i(1), alo, ahi)
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       -pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_i(2), alo, ahi)
        else
c         imaginary perturbation:
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_i(1), alo, ahi)
          call nga_add_patch (-pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_i(2), alo, ahi)
        end if  ! if .not.limag
        end if ! lifetime
      enddo                     ! ivec = 1,nvec

c Deallocate arrays no longer needed
      do ipm = 1, ncomp
         if (.not.ga_destroy(g_tmpwork(ipm))) call errquit('dim_addop
     $      _uhf: g_tmpwork GA?', 0, GA_ERR)
      enddo

100   continue
c
c   ====================================
c   Add DIM potential to the Fock matrix
c   ====================================
c
c      call ga_print(g_movecs)
      
      g_dcv = ga_create_atom_blocked(geom, basis, 'rohf_h2e3: dcv')
      xoff = 1
      voff = nclosed + nopen + 1
      xend = nvir*nclosed
      do ivec = 1, 3 ! Loop over perturbations
        alo(1) = ivec
        ahi(1) = ivec
        do ipm = 1, ncomp! Loop over +/-
c         We only add the + direction of the DIM potential to both +/- of the Fock matrix
c   Real Portion
          call nga_copy_patch('N',g_dim_r(ipm),alo,ahi,g_dcv,blo,bhi)
          call ga_scale(g_dcv, four)
          call ga_matmul_patch('n', 'n', two, zero,
     $                           g_dcv,   1, nbf, 1, nbf,
     $                           g_movecs, 1, nbf, 1, nclosed,
     $                           g_tmp1,  1, nbf, 1, nclosed)
          call ga_sync()
          call ga_matmul_patch('t', 'n', one, zero,
     $                           g_movecs, voff, nmo, 1, nbf,
     $                           g_tmp1, 1, nbf,  1, nclosed,
     $                           g_tmp2, 1, nvir, 1, nclosed)
          call ga_sync()
          call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_Ax_r(ipm),
     $                         xoff, ivec, four, '+')
c
c   Imaginary Portion
          if (lifetime) then
          call nga_copy_patch('N',g_dim_i(ipm),alo,ahi,g_dcv,blo,bhi)
          call ga_scale(g_dcv, two)
          call ga_matmul_patch('n', 'n', two, zero,
     $                           g_dcv,   1, nbf, 1, nbf,
     $                           g_movecs, 1, nbf, 1, nclosed,
     $                           g_tmp1,  1, nbf, 1, nclosed)
          call ga_sync()
          call ga_matmul_patch('t', 'n', one, zero,
     $                           g_movecs, voff, nmo, 1, nbf,
     $                           g_tmp1, 1, nbf,  1, nclosed,
     $                           g_tmp2, 1, nvir, 1, nclosed)
          call ga_sync()
c
c   ******NOTE JEM******
c   We can remove the - sign if we apply the sign change discussed in aoresponse and rohf_hessv_xx3
c   ********************
c
          call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_Ax_i(ipm),
     $                         xoff, ivec, four, '+')
          end if !lifetime
        end do !ipm = 1, 2
      end do !ivec = 1, 3
c   ========
c   Clean up
c   =======
      lfirst = .false.
      do ipm = 1,2
c        if (.not. ga_destroy(g_pmats(ipm))) call errquit
c     $     ('addop: pmats GA?', 0, GA_ERR)
c        if (.not. ga_destroy(g_pmata(ipm))) call errquit
c     $     ('addop: pmata GA?', 0, GA_ERR)
c        if (.not. ga_destroy(g_h1mat(ipm))) call errquit
c     $     ('addop: h1mat GA?', 0, GA_ERR)
        if (.not. ga_destroy(g_dim_r(ipm))) call errquit
     $     ('addop: dim_r GA?', 0, GA_ERR)
        if (lifetime) then
        if (.not. ga_destroy(g_dim_i(ipm))) call errquit
     $     ('addop: dim_i GA?', 0, GA_ERR)
        end if
      enddo                     ! ipm = 1,2
c
      if (.not. ga_destroy(g_tmp1)) call errquit
     $   ('addop: tmp1 GA?', 0, GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit
     $   ('addop: tmp2 GA?', 0, GA_ERR)
      if (.not. ga_destroy(g_dcv)) call errquit
     $   ('addop: dcv GA?',0, GA_ERR)
      write(luout,*)'end of dimqm_addop_uhf'
c
      end subroutine dimqm_addop_uhf


c ----------------------------------------------------------------------
c Subroutine used for calculating and adding in the DIM operator in 
c frequency dependent UHF calculations. This routine is highly experimental,
c and is still being worked on. Adapted from dimqm_addop. This routine
c works for the FD case with damping. 
c 
c Modeled after dimqm_addop_uhf
c
c Author: Jeff Becca, jbb5516@psu.edu, 2017
c currently not used
c ----------------------------------------------------------------------


      Subroutine dimqm_addop_uhf_damp(g_Ax_r, g_Ax_i,
     $                      ncomp, limag, lifetime, g_dens_r, g_dens_i)

c     Called from: Nothing currently, but will be uhf_hessv2_ext.
c
c     Subroutines called from: dimqm_EqmE.F
c                              dimqm_f2d.F dim_fock_xc.F
c   
c     Calculates and adds the frequency-dependent DIM potential to the 
c     response Fock matricies (real and imaginary).  Requires knowledge
c     of both the real and imaginary vectors simultaneously which is
c     why this has to be located here, unlike the static routine.
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
      integer g_Ax_r(2)    ! F matrix handle (real)      [IN/OUT]
      integer g_Ax_i(2)    ! F matrix handle (imaginary) [IN/OUT]
      integer ncomp        ! num of components (+/-)     [IN]
      logical limag        ! Imaginary perturbation?     [IN]
      logical lifetime     ! Damping or no damping
      integer g_dens_r(2)  ! Perturbed pmat              [IN]
      integer g_dens_i(2)  ! perturbed pmat IMAG         [IN]
c
c     Local variables
      integer g_tmp1, g_tmp2, g_dcv
c      integer l_dimxyz, k_dimxyz
      double precision dimxyz(3, nDIM)
c      integer l_muind, k_muind
      double precision muind(3, nDIM, 2)
      integer dims(3), chunk(3)
      character*(255) cstemp
c      integer g_pmats(2), g_pmata(2), g_h1mat(2)
      integer g_tmpwork(2)
      integer ipm
c      integer g_dens_r(2)
c      integer g_dens_i(2)
      integer alo(3), ahi(3)
      integer blo(2), bhi(2)
      integer g_dens_comp_r
      integer g_dens_comp_i
      integer xend
      double precision pre_factor
      double precision muold(3, nDIM, 2)
      
      double precision dx_r, dy_r, dz_r
      double precision dx_i, dy_i, dz_i
      double precision dsum, rmax
      external dsum
      integer i3, ivec, n
c      integer l_fld, k_fld
      double precision fld(3, nDIM, 2)
      integer g_dim_r(2)
      integer g_dim_i(2)
      integer nvir, voff, xoff
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      character*(1) direction(3)
      character*(1) pm(2)
      data direction /'x', 'y', 'z'/
      data pm /'+', '-'/
      logical firsttime
c      double precision screen(nDIM)
      double precision dimErr(3,2,2)
      double precision calcErr
      external calcErr
      integer id
c
      id = ga_nodeid()
      if (ldebug .and. id .eq. 0) then
        write(luout,*) "Start dimqm_addop_uhf"
      end if
      nvir = nmo - nclosed - nopen
      i3 = nDIM * 3
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'dim_addop: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'dim_addop: tmp2')

      dims(1) = nbf
      dims(2) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1

c
c   Allocate new arrays
c      if(.not.ma_push_get(mt_dbl,i3*2,'dimrsp:fld',l_fld,k_fld))
c     $  call errquit('malloc dimrsp:fld failed',1,MA_ERR)
c
c      if(.not.ma_push_get(mt_dbl,i3*2,'dimrsp:muind',
c     $                                            l_muind,k_muind))
c     $  call errquit('malloc dimrsp:muind failed',1,MA_ERR)
c
c      if(.not.ma_push_get(mt_dbl,i3,'dimrsp:xyz',l_dimxyz,k_dimxyz))
c     $  call errquit('malloc dimrsp:xyz failed',1,MA_ERR)
c
      if(.not.rtdb_get(dimqm_rtdb,'dimpar:coords', mt_dbl, i3, dimxyz))
     $  call errquit('get dimpar:coords failed', 1, RTDB_ERR)
c
      g_dens_comp_r = ga_create_atom_blocked(geom,basis,
     $                          'real density matrix comp')
      if (lifetime) then
        g_dens_comp_i = ga_create_atom_blocked(geom,basis,
     $                          'imag density matrix comp')
      end if
c
c
c      call dimqm_screening(dimqm_rtdb, geom, basis, dbl_mb(k_dimxyz),
c     $                       screen)
c      screen = ONE
c
c   =============================
c   Solve for induced dipoles +/-
c   =============================
c
c     Investigate the shape of g_dens_r
      write(luout,*)'g_dens_r+'
c      call ga_print(g_dens_r(1))
c      write(luout,*)'g_dens_r-'
c      call ga_print(g_dens_r(2))
c     set dimension variables for patching
      alo(2)   =  1
      ahi(2)   =  nbf
      alo(3)   =  1
      ahi(3)   =  nbf
      blo(1)   =  1
      bhi(1)   =  nbf
      blo(2)   =  1
      bhi(2)   =  nbf
c     Loop over perturbations
      do n = 1, 3
        do ipm = 1, 2
          call ga_zero(g_dens_comp_r)
          if (lifetime) call ga_zero(g_dens_comp_i)
          alo(1) = n
          ahi(1) = n
c
          write(luout,*)'debug 1'
c       Copy current perturbation into g_dens_comp
          call nga_copy_patch('N',g_dens_r(ipm), alo, ahi,
     $                          g_dens_comp_r, blo, bhi)
          write(luout,*)'debug 2'
          if (lifetime) then
            call nga_copy_patch('N',g_dens_i(ipm), alo, ahi,
     $                              g_dens_comp_i, blo, bhi)
          write(luout,*)'debug 3'
          end if
          muind = ZERO
          fld = ZERO
          firsttime = .false.
          if(.not.rtdb_get(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_r'//pm(ipm),
     $                              mt_dbl, i3, muold(:,:,1))) then
            if(id.eq.0 .and. ldebug) 
     $         write(luout,*) "First cycle, no old dipoles!"
            muold = ZERO
            firsttime = .true.
            dimqm_seeded = .false.
c            xyz_seeded(3*(n-1)+ipm) = .false.
            if(dimtol0 < 1.0d-4 .and. .not. dimqm_noseed) then
              dimtolxyz(ipm*3 - 1 + n) = 1.0d-4
              if(id.eq.0) then
                write(luout,*) "Requested tolerance below 1.0d-4"
                write(luout,*) "Setting "//direction(n)//pm(ipm)//
     $                         " dir tolerance to 1.0d-4 to start"
              end if
            end if
          else
            if(.not.rtdb_get(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_i'//pm(ipm),
     $                           mt_dbl, i3, muold(:,:,2)))
     $          call errquit('get dimqm:muold failed',1,RTDB_ERR)
          end if
          write(luout,*)'debug 4'
c         Set convergence tolerance
c          dimtol = dimtolxyz(ipm*3 - 1 + n)
c          dimqm_seeded = xyz_seeded(ipm*3 - 1 + n)
c          dimtol = 1.0d-7
c          dimqm_noseed = .true.
c          call dfill(i3*2, ZERO, dbl_mb(k_muind), 1)
c          call dfill(i3*2, ZERO, dbl_mb(k_fld), 1)
c
c       Real portion of E-Field
c        write(luout,*) "REAL"
c        call ga_print(g_dens_comp_r)
          call dimqm_EqmE(dimqm_rtdb, g_dens_comp_r, geom, basis,
     $               fld(:,:,1), dimxyz)
c
c       Imaginary portion of E-Field
c        write(luout,*) "IMAG"
c        call ga_print(g_dens_comp_i)
          if (lifetime) then
            call dimqm_EqmE(dimqm_rtdb, g_dens_comp_i, geom, basis,
     $                      fld(:,:,2), dimxyz)
          end if
c
c       Solve for induced dipoles
c  TODO: this is always called in with 2, which makes it go to complex
c        solver. Should this happen for FD not damped case?
          call dimqm_f2d(dimqm_rtdb, fld, muind, muold, dimxyz, 2,
     $                   direction(n), pm(ipm),.false.)
c
c         Write induced dipoles to RTDB
          dx_r = SUM(muind(1,:,1))
          dy_r = SUM(muind(2,:,1))
          dz_r = SUM(muind(3,:,1))
          dx_i = SUM(muind(1,:,2))
          dy_i = SUM(muind(2,:,2))
          dz_i = SUM(muind(3,:,2))
          if(id.eq.0) then
            write(luout,*) "Total induced dipole moment for "//
     $                  direction(n)//pm(ipm)//" perturbation"
            write(luout,*) "X:", dx_r, dx_i
            write(luout,*) "Y:", dy_r, dy_i
            write(luout,*) "Z:", dz_r, dz_i
            write(luout,*) ''
          end if
          dimErr(n, ipm, 1) = calcErr(i3, muold(:,:,1), muind(:,:,1))
          dimErr(n, ipm, 2) = calcErr(i3, muold(:,:,2), muind(:,:,2))
          if(id.eq.0) then
            write(luout,*) "Max error in real dipoles:",
     $                       dimErr(n, ipm, 1)
            write(luout,*) "Max error in imag dipoles:",
     $                       dimErr(n, ipm, 2)
          end if
c          if(dimErr(n, ipm, 1)/dimtol < HUNDRED
c     $              .and. dimErr(n, ipm, 2)/dimtol < HUNDRED
c     $              .and. .not. xyz_seeded(ipm*3 - 1 + n)
c     $              .and. .not. firsttime) then
c            xyz_seeded(ipm*3 - 1 + n) = .true.
c            write(luout,*) "Error within 10^2 of", dimtol, "for "//
c     $                     direction(n)//pm(ipm)//" dir"
c            write(luout,*) "Setting current "//direction(n)//pm(ipm)//
c     $                     " dir as seed"
c            write(luout,*)"Reverting tolerance back to", dimtol0
c            dimtolxyz(ipm*3 - 1 + n) = dimtol0
c          end if
          if(.not.rtdb_put(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_r'//pm(ipm),
     $                              mt_dbl, i3, muind(:,:,1)))
     $        call errquit('put dimqm:muind_p failed',1,RTDB_ERR)
          if(.not.rtdb_put(dimqm_rtdb,
     $                'dimqm:muind_'//direction(n)//'_i'//pm(ipm),
     $                           mt_dbl, i3, muind(:,:,2)))
     $        call errquit('put dimqm:muind_p failed',1,RTDB_ERR)
        end do ! ipm = 1, 2
      end do ! ivec = 1, 3
c      if(MAXVAL(dimErr) <= 1.0d-4) then
c        write(luout,*) "Dipole error below 1d-4"
c        write(luout,*) "Shutting down DIM"
c        dimqm_on = .false.
c      end if
c
c   Destroy GAs we don't need anymore
      if (.not. ga_destroy(g_dens_comp_r)) call errquit
     $    ('addop: dens_comp_r GA?',0, GA_ERR)
      if (lifetime) then
      if (.not. ga_destroy(g_dens_comp_i)) call errquit
     $    ('addop: dens_comp_i GA?',0, GA_ERR)
      end if
c      do ipm = 1,2
c        if (.not. ga_destroy(g_dens_r(ipm))) call errquit
c     $     ('addop: dens_r GA?',0, GA_ERR)
c        if (lifetime) then
c         if (.not. ga_destroy(g_dens_i(ipm))) call errquit
c     $     ('addop: dens_i GA?',0, GA_ERR)
c        endif
c      end do
c
c   Deallocate l_fld, l_muind, l_dimxyz
c      if (.not. ma_chop_stack(l_fld)) call errquit
c     $   ('addop: fld MA?', 0, MA_ERR)
c
c   ====================================================
c   Solve for DIM potential, both real and imaginary S/A
c   ====================================================
c
      dims(1) = 3
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
c
c   Real +
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_r+',chunk,
     &                        g_dim_r(1)))
     &   call errquit('addop: could not allocate g_dim_r+',1,GA_ERR)
      call ga_zero(g_dim_r(1)) 
      call fock_dim(geom, nbf, basis, 3, g_dim_r(1), 1, 1)
      call ga_symmetrize(g_dim_r(1))
c
c   Real -
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_r-',chunk,
     &                        g_dim_r(2)))
     &   call errquit('addop: could not allocate g_dim_r-',1,GA_ERR)
      call ga_zero(g_dim_r(2))
      call fock_dim(geom, nbf, basis, 3, g_dim_r(2), 2, 1)
      call ga_antisymmetrize(g_dim_r(2))
      if (lifetime) then
c
c   Imaginary +
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_i+',chunk,
     &                        g_dim_i(1)))
     &   call errquit('addop: could not allocate g_dim_i+',1,GA_ERR)
      call ga_zero(g_dim_i(1))
      call fock_dim(geom, nbf, basis, 3, g_dim_i(1), 1, 2)
      call ga_symmetrize(g_dim_i(1))
c
c   Imaginary -
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_i-',chunk,
     &                        g_dim_i(2)))
     &   call errquit('addop: could not allocate g_dim_i-',1,GA_ERR)
      call ga_zero(g_dim_i(2))
      call fock_dim(geom, nbf, basis, 3, g_dim_i(2), 2, 2)
      call ga_antisymmetrize(g_dim_i(2))
      end if
c
c   ======================================
c   Undo the symmetrization to recover +/-
c   ======================================
      blo(1)   =  nbf
      blo(2)   =  nbf
      chunk(1) =  blo(1)
      chunk(2) =  -1

      do ipm   =  1, ncomp
         write(cstemp, '(a,i1)') 'g_tmpwork_',ipm
         if (.not.nga_create(MT_DBL,2,blo,cstemp(1:11),chunk,
     $         g_tmpwork(ipm))) call errquit('dim_addop_ufh: 
     $         nga_create failed '//cstemp(1:11),0,GA_ERR)
         call ga_zero(g_tmpwork(ipm))
      enddo
      write(luout,*)'debug 5'
c     reset blo and bhi for future use
      blo(1)   =  1
      bhi(1)   =  nbf
      blo(2)   =  1
      bhi(2)   =  nbf
c
      do ivec = 1, 3
        alo(1) = ivec
        ahi(1) = ivec
c       ************
c       Real portion
c       ************
c TODO: I think the g_pmats here are being used just as temp arrays, 
c        but I really do not know for sure. 
c  all uses of g_pmats are being switched to g_tmpwork and the arrays
c  are allocated here.
        call nga_copy_patch('N',g_dim_r(1),alo,ahi,g_tmpwork(1),blo,bhi)
        call nga_copy_patch('N',g_dim_r(2),alo,ahi,g_tmpwork(2),blo,bhi)
      
        write(luout,*)'debug 6'
c
c       it might be necessary to use 0.5 here instead of 1.0
c       (note: that turned out NOT to be the case after some testing)
        pre_factor = 1.0d0
        call ga_sync()
        if (.not.limag) then
c         real perturbation: 
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_r(1), alo, ahi)
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       -pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_r(2), alo, ahi)
        write(luout,*)'debug 7'
        else
c         imaginary perturbation:
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_r(1), alo, ahi)
          call nga_add_patch (-pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_r(2), alo, ahi)
        end if  ! if .not.limag
        if (lifetime) then
c       *****************
c       Imaginary portion
c       *****************
        call nga_copy_patch('N',g_dim_i(1),alo,ahi,g_tmpwork(1),blo,bhi)
        call nga_copy_patch('N',g_dim_i(2),alo,ahi,g_tmpwork(2),blo,bhi)
c
c       it might be necessary to use 0.5 here instead of 1.0
c       (note: that turned out NOT to be the case after some testing)
        pre_factor = 1.0d0
        call ga_sync()
        if (.not.limag) then
c         real perturbation: 
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_i(1), alo, ahi)
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       -pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_i(2), alo, ahi)
        write(luout,*)'debug 8'
        else
c         imaginary perturbation:
          call nga_add_patch (pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_i(1), alo, ahi)
          call nga_add_patch (-pre_factor, g_tmpwork(1), blo, bhi,
     &       pre_factor, g_tmpwork(2), blo, bhi,
     &       g_dim_i(2), alo, ahi)
        end if  ! if .not.limag
        end if ! lifetime
      enddo                     ! ivec = 1,nvec

        write(luout,*)'debug 9'
c Deallocate arrays no longer needed
      do ipm = 1, ncomp
         if (.not.ga_destroy(g_tmpwork(ipm))) call errquit('dim_addop
     $      _uhf: g_tmpwork GA?', 0, GA_ERR)
      enddo

100   continue
c
c   ====================================
c   Add DIM potential to the Fock matrix
c   ====================================
c
c      call ga_print(g_movecs)
c        write(luout,*)'debug 10'
      
      g_dcv = ga_create_atom_blocked(geom, basis, 'rohf_h2e3: dcv')
      xoff = 1
      voff = nclosed + nopen + 1
      xend = nvir*nclosed
      do ivec = 1, 3 ! Loop over perturbations
        alo(1) = ivec
        ahi(1) = ivec
        do ipm = 1, ncomp! Loop over +/-
c         We only add the + direction of the DIM potential to both +/- of the Fock matrix
c   Real Portion
          call nga_copy_patch('N',g_dim_r(ipm),alo,ahi,g_dcv,blo,bhi)
          call ga_scale(g_dcv, four)
          call ga_matmul_patch('n', 'n', two, zero,
     $                           g_dcv,   1, nbf, 1, nbf,
     $                           g_movecs, 1, nbf, 1, nclosed,
     $                           g_tmp1,  1, nbf, 1, nclosed)
          call ga_sync()
          call ga_matmul_patch('t', 'n', one, zero,
     $                           g_movecs, voff, nmo, 1, nbf,
     $                           g_tmp1, 1, nbf,  1, nclosed,
     $                           g_tmp2, 1, nvir, 1, nclosed)
          call ga_sync()
          write(luout,*)'debug 10.5'
          call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_Ax_r(ipm),
     $                         xoff, ivec, four, '+')
c
        write(luout,*)'debug 11'
c   Imaginary Portion
          if (lifetime) then
          call nga_copy_patch('N',g_dim_i(ipm),alo,ahi,g_dcv,blo,bhi)
          call ga_scale(g_dcv, four)
          call ga_matmul_patch('n', 'n', two, zero,
     $                           g_dcv,   1, nbf, 1, nbf,
     $                           g_movecs, 1, nbf, 1, nclosed,
     $                           g_tmp1,  1, nbf, 1, nclosed)
          call ga_sync()
          call ga_matmul_patch('t', 'n', one, zero,
     $                           g_movecs, voff, nmo, 1, nbf,
     $                           g_tmp1, 1, nbf,  1, nclosed,
     $                           g_tmp2, 1, nvir, 1, nclosed)
          call ga_sync()
c
c   ******NOTE JEM******
c   We can remove the - sign if we apply the sign change discussed in aoresponse and rohf_hessv_xx3
c   ********************
c
          call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_Ax_i(ipm),
     $                         xoff, ivec, four, '+')
          end if !lifetime
        end do !ipm = 1, 2
      end do !ivec = 1, 3
c   ========
c   Clean up
c   =======
      lfirst = .false.
      do ipm = 1,2
c        if (.not. ga_destroy(g_pmats(ipm))) call errquit
c     $     ('addop: pmats GA?', 0, GA_ERR)
c        if (.not. ga_destroy(g_pmata(ipm))) call errquit
c     $     ('addop: pmata GA?', 0, GA_ERR)
c        if (.not. ga_destroy(g_h1mat(ipm))) call errquit
c     $     ('addop: h1mat GA?', 0, GA_ERR)
        if (.not. ga_destroy(g_dim_r(ipm))) call errquit
     $     ('addop: dim_r GA?', 0, GA_ERR)
        if (lifetime) then
        if (.not. ga_destroy(g_dim_i(ipm))) call errquit
     $     ('addop: dim_i GA?', 0, GA_ERR)
        end if
      enddo                     ! ipm = 1,2
c
      if (.not. ga_destroy(g_tmp1)) call errquit
     $   ('addop: tmp1 GA?', 0, GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit
     $   ('addop: tmp2 GA?', 0, GA_ERR)
      if (.not. ga_destroy(g_dcv)) call errquit
     $   ('addop: dcv GA?',0, GA_ERR)
      write(luout,*)'end of dimqm_addop_uhf'
c
      end subroutine dimqm_addop_uhf_damp
