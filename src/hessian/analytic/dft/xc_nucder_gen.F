      subroutine xc_hessian(g_dens, oactive, nactive, hess)
c
c     $Id: xc_nucder_gen.F,v 1.1 1998-09-09 13:21:15 bjohnson Exp $
c
c     Wrapper routine for XC hessian without fitting
c
c     BGJ - 8/98
c
      implicit none
c
      integer g_dens(*)         ! [input] array of handles to densities
      integer nactive           ! [input] number of active atoms
      integer oactive(nactive)  ! [input] key to active atoms
      double precision hess(*)  ! [output] XC hessian matrix
c
      call xc_deriv_gen(g_dens, oactive, nactive, hess, 0, 2)
c
      return
      end
c
      subroutine xc_cpks_rhs(g_dens, oactive, nactive, g_rhs)
c
c     $Id: xc_nucder_gen.F,v 1.1 1998-09-09 13:21:15 bjohnson Exp $
c
c     Wrapper routine for XC CPKS RHS without fitting
c
c     BGJ - 8/98
c
      implicit none
c
      integer g_dens(*)         ! [input] array of handles to densities
      integer nactive           ! [input] number of active atoms
      integer oactive(nactive)  ! [input] key to active atoms
      integer g_rhs(*)          ! [output] XC CPKS RHS matrices
c
      call xc_deriv_gen(g_dens, oactive, nactive, 0d0, g_rhs, 3)
c
      return
      end
c
      subroutine xc_deriv_gen(g_dens, oactive, nactive, hess, g_rhs,
     &                        calc_type)
c
c     $Id: xc_nucder_gen.F,v 1.1 1998-09-09 13:21:15 bjohnson Exp $
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
#include "../nwdft/include/cdft.fh"
c!!! BGJ
#include "bgj.fh"
c!!! BGJ
c
      integer g_dens(*)         ! [input] array of handles to densities
      integer nactive           ! [input] number of active atoms
      integer oactive(nactive)  ! [input] key to active atoms
      double precision hess(*)  ! [output] XC hessian matrix
      integer g_rhs(*)          ! [output] XC CPKS RHS matrices
      integer calc_type         ! [input] calculation type
                                !   1 = XC gradient (*** not ready yet ***)
                                !   2 = XC hessian
                                !   3 = XC CPKS RHS
c
c     Local declarations
c
      logical do_hess, do_cpks_r
      integer rtdb
      integer natoms
      integer l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &        l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom,
     &        nscr, l_scr, k_scr, i
      character*3 name
c
      if (bgj_print() .gt. 0)
     &     write(*,*)'--------Entered xc_deriv_gen-------------'
      do_hess = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      if (.not. (do_hess .or. do_cpks_r))
     &     call errquit('xc_deriv_gen: illegal calculation type',0)
      if (bgj_print() .gt. 0)
     &     write(*,*)'logicals ', do_hess, do_cpks_r
      rtdb = bgj_get_rtdb_handle()
      if (bgj_print() .gt. 0)
     &     write(*,*)'rtdb geom ipol',rtdb,geom,ipol
      if (.not. geom_ncent(geom, natoms))
     &   call errquit('xc_deriv_gen: geom_ncent failed',1)
c      write(*,*)'natoms',natoms
      if (.not. bgj_get_scf_method(name))
     $     call errquit('get_scf_method failed',0)
c
c     Initialize output quantities
c
      if (do_hess) then
         call dfill(9*natoms*natoms, 0.0d0, hess, 1)
      else if (do_cpks_r) then
         do i = 1, ipol*3*nactive
            call ga_zero(g_rhs(i))
         enddo
      endif
c
c     Set up indexing arrays needed by dftgh_quadv0
c
      call xc_init_index(g_dens, l_cntoce, k_cntoce, l_cntobfr,
     &                   k_cntobfr, l_cetobfr, k_cetobfr,
     &                   l_rdens_atom, k_rdens_atom)
c
      call dftgh_quadv0(rtdb, g_dens, 0d0, hess, g_rhs, calc_type,
     &                  natoms, oactive, nactive, dbl_mb(k_rdens_atom),
     &                  int_mb(k_cntoce), int_mb(k_cntobfr),
     &                  int_mb(k_cetobfr))
c
c     Clean up indexing arrays
c
      call xc_exit_index(l_cntoce, l_cntobfr, l_cetobfr, l_rdens_atom)
c
      if (bgj_print() .gt. 0) then
         if (do_hess) then
            write(6,*) 'xc_deriv_gen: xc hessian'
            call output(hess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
         else if (do_cpks_r) then
            do i = 1, ipol*3*nactive
               write(*,*)'*** xc cpks rhs matrix',i
               call ga_print(g_rhs(i))
            enddo
         endif
         write(*,*)'--------Leaving xc_deriv_gen-------------'
      endif
      return
      end
