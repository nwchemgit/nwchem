      subroutine xc_hessian(geom,rtdb, nbf,basis,
     ,     g_dens, nactive,oactive,  hess, scftype)
c
c     $Id: xc_nucder_gen.F,v 1.8 2000-08-25 20:10:15 edo Exp $
c
c     Wrapper routine for XC hessian without fitting
c
c     BGJ - 8/98
c
      implicit none
c
      integer geom,rtdb,nbf,basis         ! [in]
      character *4 scftype    ! [in]
      integer g_dens(*)         ! [input] array of handles to densities
      integer nactive           ! [input] number of active atoms
      logical oactive(nactive)  ! [input] key to active atoms
      double precision hess(*)  ! [output] XC hessian matrix
c
      integer ipol
c
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
c     $Id: xc_nucder_gen.F,v 1.8 2000-08-25 20:10:15 edo Exp $
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
c     $Id: xc_nucder_gen.F,v 1.8 2000-08-25 20:10:15 edo Exp $
c
c     General routine for XC hessian and CPKS RHS without fitting
c     Should not be called directly by client, but through wrapper routines
c
c     BGJ - 8/98
c
      implicit none
c
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c!!! BGJ
#include "bgj.fh"
#include "global.fh"
c!!! BGJ
c
      integer geom,rtdb,nbf,basis
      integer ipol              ! [in]
      integer g_dens(*)         ! [input] array of handles to densities
      integer nactive           ! [input] number of active atoms
      logical oactive(nactive)  ! [input] key to active atoms
      double precision hess(*)  ! [output] XC hessian matrix
      integer g_rhs(ipol*3,*)   ! [output] XC CPKS RHS matrices
      integer calc_type         ! [input] calculation type
                                !   1 = XC gradient (*** not ready yet ***)
                                !   2 = XC hessian
                                !   3 = XC CPKS RHS
c
c     Local declarations
c
      logical do_hess, do_cpks_r,grid_on_file
      integer natoms
      integer l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &        l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom, i, j, i1
c !!! BGJ test
c
      integer l_ti_hess, k_ti_hess, g_ti(6)
c
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &     grid_on_file))then
         grid_on_file = .false.
      endif
      if (bgj_print() .gt. 0)
     &     write(*,*)'--------Entered xc_deriv_gen-------------'
      do_hess = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      if (.not. (do_hess .or. do_cpks_r))
     &     call errquit('xc_deriv_gen: illegal calculation type',0)
      if (bgj_print() .gt. 0)
     &     write(*,*)'logicals ', do_hess, do_cpks_r
      if (bgj_print() .gt. 0)
     &     write(*,*)'rtdb geom ipol',rtdb,geom,ipol
      if (.not. geom_ncent(geom, natoms))
     &   call errquit('xc_deriv_gen: geom_ncent failed',1)
c      write(*,*)'natoms',natoms
#if 0
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
      if(grid_on_file) then
         call dftgh_gridv0(rtdb,geom,basis, ipol,nbf,
     .        g_dens, 0d0, hess, g_rhs, calc_type,
     &        natoms, oactive, nactive, dbl_mb(k_rdens_atom),
     &        int_mb(k_cetobfr))
      else
         call dftgh_quadv0(rtdb,geom,basis, ipol,
     .        g_dens, 0d0, hess, g_rhs, calc_type,
     &        natoms, oactive, nactive, dbl_mb(k_rdens_atom),
     &        int_mb(k_cetobfr))
      endif
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
            if (oactive(i)) then
               if (i1 .lt. i) then
                  do j = 1, ipol*3
                     g_rhs(j,i) = g_rhs(j,i1)
                  enddo
               endif
               i1 = i1 - 1
            endif
         enddo
      endif
c
      if (bgj_print() .gt. 0) then
         if (do_hess) then
c
            write(6,*) 'xc_deriv_gen: xc hessian'
            call output(hess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
c
c     TI check
c
            if (.not. MA_push_get(MT_DBL, natoms*natoms*9, 'ti_hess',
     &           l_ti_hess, k_ti_hess))
     &           call errquit('j_nucder_gen: could not alloc ti_hess',1)
            call dcopy(natoms*natoms*9, hess, 1, dbl_mb(k_ti_hess), 1)
            call hessti(dbl_mb(k_ti_hess),natoms)
            if (.not. MA_pop_stack(l_ti_hess))
     &           call errquit('j_nucder_gen: could not pop ti_hess',1)
c
         else if (do_cpks_r) then
c
            do i = 1, natoms
               if (oactive(i)) then
                  do j = 1, ipol*3
                     write(*,*)'*** xc cpks rhs matrix',j,i
                     call ga_print(g_rhs(j,i))
                  enddo
               endif
            enddo
c
c     TI check
c
            do i = 1, ipol*3
              if (.not. ga_duplicate(g_dens(1), g_ti(i), 'RHS TI'))
     &          call errquit('xc_deriv_gen: could not create temp GA',0)
              call ga_zero(g_ti(i))
            enddo
            do i = 1, natoms
               if (oactive(i)) then
                  do j = 1, ipol*3
                      call ga_add(1d0, g_rhs(j,i), 1d0, g_ti(j),
     &                            g_ti(j))
                  enddo
               endif
            enddo
            do j = 1, ipol*3
               write(*,*)'*** TI check',j
               call ga_print(g_ti(j))
               if (.not. ga_destroy(g_ti(j)))
     &         call errquit('xc_deriv_gen: could not destroy temp GA',0)
            enddo
         endif
         write(*,*)'--------Leaving xc_deriv_gen-------------'
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
