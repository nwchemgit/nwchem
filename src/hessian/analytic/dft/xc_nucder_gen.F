      subroutine xc_hessian(geom,rtdb, nbf,basis,
     ,     g_dens, nactive,oactive,  hess, scftype)
c
c     Wrapper routine for XC hessian without fitting
c
c     BGJ - 8/98
c
      implicit none
c
      integer geom,rtdb,nbf,basis         ! [in]
      character *4 scftype      ! [in]
      integer g_dens(*)         ! [input] array of handles to densities
      integer nactive           ! [input] number of active atoms
      logical oactive(nactive)  ! [input] key to active atoms
      double precision hess(*)  ! [output] XC hessian matrix
c
      integer ipol
      external xc_gotxc
      logical xc_gotxc
c
      if(.not.xc_gotxc()) return
      ipol=1
      if(scftype.eq.'UHF') ipol=2
c
      call xc_deriv_gen(geom,rtdb,ipol,nbf,basis,
     ,     g_dens, nactive, oactive, hess, 0, 2)
c
      return
      end
c
      subroutine xc_cpks_rhs(geom,rtdb,nbf,basis,
     ,     g_dens, nactive, oactive, g_rhs,scftype)
c
c     Wrapper routine for XC CPKS RHS without fitting
c
c     BGJ - 8/98
c
      implicit none
c
      integer geom,rtdb,nbf,basis
      character *(*) scftype    ! [in]
      integer g_dens(*)         ! [input] array of handles to densities
      integer nactive           ! [input] number of active atoms
      logical oactive(nactive)  ! [input] key to active atoms
      integer g_rhs(*)          ! [output] XC CPKS RHS matrices
c
      integer ipol
c
      ipol=1
      if(scftype.eq.'UHF') ipol=2
c
      call xc_deriv_gen(geom,rtdb,ipol,nbf,basis,
     ,     g_dens,  nactive, oactive, 0d0, g_rhs, 3)
c
      return
      end
c
      subroutine xc_deriv_gen(geom,rtdb,ipol,nbf,basis,
     ,     g_dens, nactive, oactive,  hess, g_rhs,
     &                        calc_type)
c
c     General routine for XC hessian and CPKS RHS without fitting
c     Should not be called directly by client, but through wrapper routines
c
c     BGJ - 8/98
c
      implicit none
#include "errquit.fh"
c
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c!!! BGJ
#include "bgj.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
c!!! BGJ
c
      integer geom,rtdb,nbf,basis
      integer ipol              ! [input] =1 for RKS and =2 for UKS
      integer g_dens(*)         ! [input] array of handles to densities
      integer nactive           ! [input] number of active atoms
      logical oactive(nactive)  ! [input] key to active atoms
      double precision hess(*)  ! [output] XC hessian matrix
      integer g_rhs(3,nactive,ipol) ! [output] XC CPKS RHS matrices
      integer calc_type         ! [input] calculation type
                                !   1 = XC gradient (*** not ready yet ***)
                                !   2 = XC hessian
                                !   3 = XC CPKS RHS
c
c     Local declarations
c
      logical do_hess, do_cpks_r,grid_written,grid_reopen
      external grid_reopen
      integer natoms
      integer l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &        l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom, i, j, i1
c !!! BGJ test
c
      integer l_ti_hess, k_ti_hess, g_ti(3)
      integer alo(3),ahi(3),blo(3),bhi(3)
c
      external xc_gotxc
      logical xc_gotxc
      logical oprint
c
      if(.not.xc_gotxc()) return
      oprint= util_print('xc_hessian',print_debug)
      do_hess = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      if (.not. (do_hess .or. do_cpks_r))
     &     call errquit('xc_deriv_gen: illegal calculation type',0,
     &       INPUT_ERR)
      if (oprint.and.ga_nodeid().eq.0) then
          write(luout,*)'--------Entered xc_deriv_gen-------------'
          write(luout,*)'logicals ', do_hess, do_cpks_r
          write(luout,*)'rtdb geom ipol',rtdb,geom,ipol
       endif
      if (.not. geom_ncent(geom, natoms))
     &   call errquit('xc_deriv_gen: geom_ncent failed',1, GEOM_ERR)
c      write(*,*)'natoms',natoms
#if 0
  Broken because of 3-d g_rhs
c     !!! Commented out so that routines will accumulate, allowing
c     !!! more efficient use of storage
c
c     Initialize output quantities
c
      if (do_hess) then
         call dfill(9*natoms*natoms, 0.0d0, hess, 1)
      else if (do_cpks_r) then
c
c        Collapse GA handles for RHS since the XC code processes them
c        most efficiently that way
c
         i1 = 0
         do i = 1, natoms
            if (oactive(i)) then
               i1 = i1 + 1
               do j = 1, ipol*3
                  call ga_zero(g_rhs(j,i))
                  if (i1 .lt. i) g_rhs(j,i1) = g_rhs(j,i)
               enddo
            endif
         enddo
      endif
#endif
c
c     Set up indexing arrays needed by dftgh_quadv0
c
      call xc_init_index(geom,ipol,nbf,basis,
     ,     g_dens, l_cntoce, k_cntoce, l_cntobfr,
     &                   k_cntobfr, l_cetobfr, k_cetobfr,
     &                   l_rdens_atom, k_rdens_atom)
c
c
c     check if we can get the grid from a file
c         
      grid_written=grid_reopen(geom)
c
      call dftgh_gridv0(rtdb,geom,basis, ipol,nbf,
     .        g_dens, 0d0, hess, g_rhs, calc_type,
     &        natoms, oactive, nactive, dbl_mb(k_rdens_atom),
     &     int_mb(k_cetobfr), 0, .false.,
     &     0d0, .false.)
c
c     Clean up indexing arrays
c
      call xc_exit_index(l_cntoce, l_cntobfr, l_cetobfr, l_rdens_atom)
c
      if (do_cpks_r) then
c
c        Re-expand GA handles for RHS
c
         i1 = nactive
         do i = natoms, 1, -1
#if 0
            if (oactive(i)) then
               if (i1 .lt. i) then
                  do j = 1, ipol*3
                     g_rhs(j,i) = g_rhs(j,i1)
                  enddo
               else
                  call errquit('Fix collapse/expand in xc_nucder_gen',0,
     &       UNKNOWN_ERR)
               endif
               i1 = i1 - 1
            endif
#else
            if (.not.oactive(i)) then
               call errquit('Fix collapse/expand in xc_nucder_gen',0,
     &       UNKNOWN_ERR)
            endif
#endif
         enddo
      endif
c
      if (oprint) then
         if (do_hess) then
c
            if(ga_nodeid().eq.0) then
            write(luout,*) 'xc_deriv_gen: xc hessian'
            call output(hess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
            endif
c
c     TI check
c
            if (.not. MA_push_get(MT_DBL, natoms*natoms*9, 'ti_hess',
     &           l_ti_hess, k_ti_hess))
     &           call errquit('j_nucder_gen: could not alloc ti_hess',1,
     &       MA_ERR)
            call dcopy(natoms*natoms*9, hess, 1, dbl_mb(k_ti_hess), 1)
            call hessti(dbl_mb(k_ti_hess),natoms)
            if (.not. MA_pop_stack(l_ti_hess))
     &           call errquit('j_nucder_gen: could not pop ti_hess',1,
     &       MA_ERR)
c
         else if (do_cpks_r) then
c
!            do i1 = 1, ipol
!            do i = 1, natoms
!               if (oactive(i)) then
!                  do j = 1, 3
!                     write(*,*)'*** xc cpks rhs matrix',j,i,i1
                     write(*,*)'*** xc cpks rhs matrix'
                     call ga_print(g_rhs)
!                     call ga_print(g_rhs(j,i,i1))
!                  enddo
!               endif
!            enddo
!            enddo
c
c     TI check
c
            do i = 1, 3
              if (.not. ga_duplicate(g_dens(1), g_ti(i), 'RHS TI'))
     &          call errquit('xc_deriv_gen: could not create temp GA',0,
     &       GA_ERR)
            enddo
            do i1 = 1, ipol
            do j = 1, 3
              call ga_zero(g_ti(j))
            enddo
            blo(1)=1
            bhi(1)=nbf
            blo(2)=1
            bhi(2)=nbf
            alo(2)=1
            ahi(2)=nbf
            alo(3)=1
            ahi(3)=nbf
            do i = 1, natoms
               if (oactive(i)) then
                  do j = 1, 3
                     ahi(1)=j+(i-1)*3+(i1-1)*3*natoms
                     alo(1)=ahi(1)
                     call nga_add_patch(1d0,g_rhs,alo,ahi,
     &                    1d0,g_ti(j),blo,bhi,
     &                    g_ti(j),blo,bhi)
                  enddo
               endif
            enddo
            do j = 1, 3
               write(*,*)'*** TI check',j, i1
               call ga_print(g_ti(j))
            enddo
            enddo
            do j = 1, 3
               if (.not. ga_destroy(g_ti(j)))
     &         call errquit('xc_deriv_gen: could not destroy temp GA',0,
     &       GA_ERR)
            enddo
         endif
         write(luout,*)'--------Leaving xc_deriv_gen-------------'
c
      endif
      return
      end
c !!! BGJ test
      subroutine hessti(hess, natoms)
c
c     Tests hessian TI by summing rows
c
      implicit none
c
      integer natoms
      double precision hess(3,natoms,3*natoms)
      integer i, j
      do i = 2, natoms
         do j = 1, 3*natoms
            hess(1,1,j) = hess(1,1,j) + hess(1,i,j)
            hess(2,1,j) = hess(2,1,j) + hess(2,i,j)
            hess(3,1,j) = hess(3,1,j) + hess(3,i,j)
         enddo
      enddo
      write(*,*)'*** TI check'
c      call output(hess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
      call output(hess,1,3,1,3*natoms,3*natoms,3*natoms,1)
      return
      end
c !!! BGJ test
c $Id$
