      subroutine grad_dens_so
     $     (geom, basis, g_dens, g_wdens, g_eigen_diag, 
     $     evals, occ, ndens, nbf, nbf2, nopen, nclosed, scftype, 
     $     movecs, omp2, pdm2, coeff)

c     $Id: grad_dens_so.F,v 1.2 2000-06-13 23:34:32 zzhang Exp $

C     calculate the density and energy-weighted density matrices from the
C     movecs (and Lagrangian in ROHF case) on disk

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "eaf.fh"

      integer ndens,geom, basis, g_dens(ndens), g_wdens, g_eigen_diag,
     $     nbf, nbf2, nopen, nclosed, g_densso(2)
      double precision evals(*), occ(*), pdm2(*), coeff(*)
      character*(*) scftype
      character*(*) movecs
      double precision offset   ! EAF file offset
      character*255 pdmfilename ! EAF file name
      integer pdmfile           ! EAF file handle
      integer l_buf, k_buf

      logical omp2

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read_so, file_read_ga
      external movecs_read_so, file_read_ga

      character*255 fname       ! filename for lagrangian (ROHF) 
                                ! or second-order contributions (MP2)

C     local variables
      integer nocc, l_ind, k_ind, l_ind2, k_ind2, g_vecs(2), i
      double precision zero, one, mone, two, mtwo
      logical status

      parameter ( mtwo = -2.D0 )
      parameter ( mone = -1.D0 )
      parameter ( zero =  0.D0 )
      parameter ( one  =  1.D0 )
      parameter ( two  =  2.D0 )

C     get MO vectors from file (taken from Roberts scf_movecs_read)

      status = ga_create(mt_dbl, nbf2, nbf2,'Movecs Re', 0, 0, 
     &                   g_vecs(1))
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' grad_dens_so: ', 'Movecs Re'
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('grad_dens_so: ga_create ', 0)
      endif

      status = ga_create(mt_dbl, nbf2, nbf2,'Movecs Im', 0, 0, 
     &         g_vecs(2))
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' grad_dens_so: ', 'Movecs Im'
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('grad_dens_so: ga_create ', 0)
      endif

      if (.not. movecs_read_so(movecs, occ, evals, g_vecs)) 
     $     call errquit('grad_dens_so: could not read mo vectors', 110)
c
      if (scftype .eq. 'UHF') then
         
C     energy-weighted density
         
C     create eigenvalue diagonal matrix
cso   nocc = nclosed + nopen ! == nalpha
         nocc = 2*nclosed + nopen 
         if (ga_nodeid() .eq. 0) then
            status = MA_alloc_get(MT_INT, nocc, 'indices', l_ind, k_ind)
            status = MA_alloc_get(MT_INT, nocc, 'indices',l_ind2,k_ind2)
            do i=1, nocc
               int_mb( k_ind+i-1 ) = i
               int_mb( k_ind2+i-1 ) = i
            end do
            call ga_scatter ( g_eigen_diag, evals, int_mb(k_ind),
     $           int_mb(k_ind2), nocc )
         end if
         status = ga_create(mt_dbl,nbf2,nbf2,'density matrix', 0, 0, 
     &                   g_densso(1))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' grad_dens_so: ', 'DenMx Re'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('grad_dens_so: ga_create ', 0)
         endif
         
         
         status = ga_create(mt_dbl,nbf2,nbf2,'density matrix', 0, 0, 
     &        g_densso(2))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' grad_dens_so: ', 'DenMx Im'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('grad_dens_so: ga_create ', 0)
         endif
         
C     multiply to weighted density
         call dft_densm_sowght(g_densso,g_vecs,g_eigen_diag,nbf,nocc)

cso   get total weighted density matrix 
         call ga_dens_sf(g_dens, g_densso, nbf)
         call ga_dadd ( one, g_dens(1), one, g_dens(2), g_wdens )

C     density matrix
         call dft_densm_so(g_densso,g_vecs,nbf,nocc)
         call ga_dens_sf(g_dens, g_densso, nbf)
         call ga_dens_so(g_dens(3), g_densso, nbf, 'z')
         call ga_dens_so(g_dens(4), g_densso, nbf, 'y')
         call ga_dens_so(g_dens(5), g_densso, nbf, 'x')

         if ( .not. omp2 ) then
C     D+ and D- (Formalism from Dupuis & King, JCP 68(9), p4000
            call ga_dadd ( one, g_dens(1),one, g_dens(2),g_dens(1))
            call ga_scale(g_dens(2),2.0d0)
            call ga_dadd ( one, g_dens(1),mone,g_dens(2),g_dens(2))
            call ga_symmetrize(g_dens(1))
            call ga_symmetrize(g_dens(2))
         end if                  !  UMP2

c
         call ga_sync()
c
         call ga_symmetrize(g_wdens)
         call ga_symmetrize(g_dens(1))
      end if
         
c     

      status = ga_destroy ( g_vecs(1) )
      if ( scftype .eq. 'UHF' ) then
         status = ga_destroy ( g_vecs(2) )
      end if
      status = ga_destroy ( g_densso(1) )
      status = ga_destroy ( g_densso(2) )

      return
      end

