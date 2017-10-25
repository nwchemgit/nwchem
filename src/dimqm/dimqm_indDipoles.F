      subroutine dimqm_indDipoles(rtdb, nbf, geom, basis,
     $                              g_dens, luhf_fd)
c
c     Determines the induced dipoles for the given density matrix 
c    
c     Called from: dimqm_rohf.F
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
c
c     Input Variables
      integer rtdb                ! RTDB handle
      integer g_dens              ! perturbed density matrix handle
      integer nbf                 ! number of basis functions
      integer geom                ! geometry handle
      integer basis               ! basis set handle
c
c     Local variables
      double precision dx, dy, dz
c      double precision dsum
c      external dsum
      integer i3, n
c      integer icmplx
c      integer l_fld, k_fld
      double precision fld(3, nDIM)
c      integer nNuc
c      integer l_qmxyz, k_qmxyz
c      integer l_dimxyz, k_dimxyz
      double precision dimxyz(3, nDIM)
c      integer l_muind, k_muind
      double precision muind(3, nDIM)
      double precision muold(3, nDIM)
      integer g_dens_comp
      integer alo(3), ahi(3)
      integer blo(2), bhi(2)
      double precision dimErr(3)
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      character*(1) direction(3)
      data direction /'x', 'y', 'z'/
      double precision maxErr, calcErr
      external calcErr
      logical firsttime
      integer id
c jbecca START
      logical luhf_fd   ! variable for using this routine for FD Open Shell
c  jbecca END
c      double precision screen(nDIM)

c
c     Common variables used
c
c     logical lcmplx

      id = ga_nodeid()
      i3 = nDIM * 3
c      icmplx = REALMATH
c      if(lfreq .and. lcmplx) icmplx = COMPLEXMATH
c      call ga_print(g_dens)
c
c     Allocate memory for arrays
c      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:fld',l_fld,k_fld))
c     $  call errquit('malloc dimrsp:fld failed',1,MA_ERR)
c
c
c      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:muind',
c     $                                            l_muind,k_muind))
c     $  call errquit('malloc dimrsp:muind failed',1,MA_ERR)
c
c
c      if(.not.ma_push_get(mt_dbl,i3,'dimrsp:xyz',l_dimxyz,k_dimxyz))
c     $  call errquit('malloc dimrsp:xyz failed',1,MA_ERR)
c
      if(.not.rtdb_get(rtdb, 'dimpar:coords', mt_dbl, i3, dimxyz))
     $  call errquit('get dimpar:coords failed', 1, RTDB_ERR)
c
      g_dens_comp = ga_create_atom_blocked(geom,basis,
     $                          'density matrix comp')

c      call dimqm_screening(rtdb, geom, basis, dbl_mb(k_dimxyz),
c     $                     screen)
c
      dimErr = ZERO
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
c
      do n = 1, 3 ! Loop over perturbations
        call ga_zero(g_dens_comp)
        alo(1) = n
        ahi(1) = n 
c
c       Copy current perturbation into g_dens_comp
        call nga_copy_patch('N',g_dens, alo, ahi,
     $                          g_dens_comp, blo, bhi)
c
        call ga_sync() 
c
c       Initizalize arrays
        muind = ZERO
        fld = ZERO
        firsttime = .false.
        if(.not.rtdb_get(rtdb, 'dimqm:muind_'//direction(n), mt_dbl,
     $                   i3, muold)) then
          if(id.eq.0 .and. ldebug) write(luout,*)
     $                  "First cycle, no old dipoles! Setting to zero"
          muold = ZERO
          firsttime = .true.
          xyz_seeded(n) = .false.
          if(dimtol0 < 1.0d-4 .and. .not. dimqm_noseed) then
            dimtolxyz(n) = 1.0d-4
            if(id.eq.0 .and. ldebug ) then
              write(luout,*) "Requested tolerance below 1.0d-4"
              write(luout,*) "Setting "//direction(n)//" dir tolerance"
     $                                // " to 1.0d-4 to start"
            end if
          end if
        end if
c       Set covergence tolerance
        dimtol = dimtolxyz(n)
        dimqm_seeded = xyz_seeded(n)
c        call dfill(i3*icmplx, ZERO, dbl_mb(k_muind), 1)
c        call dfill(i3*icmplx, ZERO, dbl_mb(k_fld), 1)
c
c       Calculate the QM E-field from the QM density
c  jbecca START
c        write(luout,*)'cpks dens'
c        call ga_print(g_dens_comp)
c  jbecca END
        call dimqm_EqmE(rtdb, g_dens_comp, geom, basis, fld, dimxyz)
c
c       Calculate the dipoles induced by the QM density
c  jbecca START
        if (luhf_fd) then
            call dimqm_f2d(rtdb, fld, muind, muold, dimxyz, 2,
     $                 direction(n), ' ',.false.)
        else
            call dimqm_f2d(rtdb, fld, muind, muold, dimxyz, 1,
     $                 direction(n), ' ',.false.)
        endif     !luhf_fd
c jbecca END
c
c       Write induced dipoles to RTDB
c        dx = dsum(nDIM, dbl_mb(k_muind), 3)
c        dy = dsum(nDIM, dbl_mb(k_muind+1), 3)
c        dz = dsum(nDIM, dbl_mb(k_muind+2), 3)
        dx = SUM(muind(:, 1))
        dy = SUM(muind(:, 2))
        dz = SUM(muind(:, 3))
        if(id.eq.0.and.ldebug) then
          write(luout,*) "Total induced dipole moment for "//
     $                  direction(n)//" perturbation"
          write(luout,*) dx, dy, dz
          write(luout,*) ''
        end if
        dimErr(n) = calcErr(i3, muold, muind)
        if(id.eq.0.and.ldebug) 
     $      write(luout,*) "Max error in dipoles:", dimErr(n)
        if(id.eq.0.and.ldebug) write(luout,*) ""
c
        if(dimErr(n) < 2.0d-1 .and. .not. xyz_seeded(n)
     $                        .and. .not. firsttime
     $                        .and. .not. dimqm_noseed) then
          xyz_seeded(n) = .true.
          if(id.eq.0.and.ldebug) then
            write(luout,*) "Error below 2.0d-1 for "//direction(n)
     $                                             //" dir"
            write(luout,*)"Setting current "//direction(n)
     $                                     //" dir as seed"
            write(luout,*)"Reverting tolerance back to", dimtol0
          end if
          dimtolxyz(n) = dimtol0
        end if        
        if(.not.rtdb_put(rtdb, 'dimqm:muind_'//direction(n), mt_dbl,
     $                                 i3, muind))
     $    call errquit('put dimqm:muind_p failed',1,RTDB_ERR)
      end do ! n = 1, 3
c      if(MAXVAL(dimErr) <= 1.0d-4) then
c        write(luout,*) "Dipole error below 1d-4"
c        write(luout,*) "Shutting down DIM"
c        dimqm_on = .false.
c      end if
c
c     Clean up arrays
      if(.not.ga_destroy(g_dens_comp))
     $  call errquit('ga_destroy failed g_dens_comp',1,GA_ERR)
c      if(.not.ma_chop_stack(l_fld))
c     $  call errquit('ma chop at l_fld failed',1,MA_ERR)
c      write(luout,*) "End Dens2Dipoles"
c
      end subroutine dimqm_indDipoles
