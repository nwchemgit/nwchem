c
c     Driver routines which call generalized dftg_grad2 (dftg_grad2_gen)
c     to create derivative J matrices and J hessian (without fitting)
c
c     BGJ - 9/98
c
c     $Id: j_nucder_gen.F,v 1.2 1998-09-25 11:46:00 bjohnson Exp $
c
      subroutine j_hessian(g_dens, oactive, nactive, hess)
c
c     Wrapper routine for J hessian without fitting
c
      implicit none
c
      integer g_dens            ! [input] GA handle for total density
      integer nactive           ! [input] number of active atoms
      logical oactive(nactive)  ! [input] key to active atoms
      double precision hess(*)  ! [output] J hessian matrix
c
      call j_deriv_gen(g_dens, oactive, nactive, hess, 0, 2)
c
      return
      end
c
c
      subroutine j_cpks_rhs(g_dens, oactive, nactive, g_rhs)
c
c     Wrapper routine for J CPKS RHS without fitting
c
      implicit none
c
      integer g_dens            ! [input] handle to total density
      integer nactive           ! [input] number of active atoms
      logical oactive(nactive)  ! [input] key to active atoms
      integer g_rhs(*)          ! [output] J CPKS RHS matrices
c
      call j_deriv_gen(g_dens, oactive, nactive, 0d0, g_rhs, 3)
c
      return
      end
c
c
      subroutine j_deriv_gen(g_dens, oactive, nactive, hess, g_rhs,
     &                       calc_type)
c
c     $Id: j_nucder_gen.F,v 1.2 1998-09-25 11:46:00 bjohnson Exp $
c
c     General routine for J hessian and CPKS RHS without fitting
c     Should not be called directly by client, but through wrapper routines
c
      implicit none
c
#include "geom.fh"
#include "mafdecls.fh"
#include "cdft.fh"
c!!! BGJ
#include "bgj.fh"
#include "global.fh"
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c!!! BGJ
c
      integer g_dens            ! [input] handle to total density
      integer nactive           ! [input] number of active atoms
      logical oactive(nactive)  ! [input] key to active atoms
      double precision hess(*)  ! [output] J hessian matrix
      integer g_rhs(3,*)        ! [output] J CPKS RHS matrices
      integer calc_type         ! [input] calculation type
                                !   1 = J gradient (*** not ready yet ***)
                                !   2 = J hessian
                                !   3 = J CPKS RHS
c
c     Local declarations
c
      logical do_hess, do_cpks_r
      integer rtdb
      integer natoms
      integer l_d_ija, k_d_ija, l_d_kla, k_d_kla, ldim, l_scr, k_scr,
     &        lscr, l_buf, k_buf, lbuf, l_s34, k_s34
      integer i, nproc
      double precision tol2e
      integer g_t
c !!! BGJ test
      integer l_ti_hess, k_ti_hess, g_ti(3)
      integer j
c !!! BGJ test
c
      if (bgj_print() .gt. 0)
     &     write(*,*)'--------Entered j_deriv_gen-------------'
      do_hess = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      if (.not. (do_hess .or. do_cpks_r))
     &     call errquit('j_deriv_gen: illegal calculation type',0)
      if (bgj_print() .gt. 0)
     &     write(*,*)'logicals ', do_hess, do_cpks_r
      rtdb = bgj_get_rtdb_handle()
      if (bgj_print() .gt. 0)
     &     write(*,*)'rtdb geom',rtdb,geom
      if (.not. geom_ncent(geom, natoms))
     &   call errquit('j_deriv_gen: geom_ncent failed',1)
c      write(*,*)'natoms',natoms
c
c     Initialize output quantities
c
      if (do_hess) then
         call dfill(9*natoms*natoms, 0.0d0, hess, 1)
      else if (do_cpks_r) then
         do i = 1, natoms
            if (oactive(i)) then
               call ga_zero(g_rhs(1,i))
               call ga_zero(g_rhs(2,i))
               call ga_zero(g_rhs(3,i))
            endif
         enddo
      endif
c
c     Set up temporary arrays needed by dftg_grad2_gen
c
      if (calc_type.eq.2) then
         write(*,*)'*** j_nucder_alloc: is hess buf big enough?'
      endif
      call j_nucder_alloc(l_d_ija, k_d_ija, l_d_kla, k_d_kla,
     &                    ldim, l_scr, k_scr, lscr, l_buf, k_buf,
     &                    lbuf, l_s34, k_s34, natoms)
c
      nproc = ga_nnodes()
      tol2e = 10.d0**(-itol2e)
      call dftg_grad2_gen(dbl_mb(k_d_ija), dbl_mb(k_d_kla), 0d0, hess,
     &                    g_rhs, g_dens, ldim, ao_bas_han, nproc,
     &                    natoms, dbl_mb(k_scr), lscr, dbl_mb(k_buf),
     &                    lbuf, tol2e, oactive, nactive, geom,
     &                    dbl_mb(k_s34), oskel, calc_type)
c
c     Clean up temporary arrays
c
      call j_nucder_dealloc(l_d_ija, l_d_kla, l_scr, l_buf, l_s34)
c
c     Symmetrize matrices across the diagonal by adding the transpose;
c     the lower level routines have set things up so that this will
c     produce the correct result
c
      if (do_cpks_r) then
c
         if (.not. ga_duplicate(g_dens, g_t, 'RHS transpose'))
     &       call errquit('j_deriv_gen: could not create temp GA',0)
c
         do i = 1, natoms
            if (oactive(i)) then
               do j = 1, 3
                  call ga_transpose(g_rhs(j,i), g_t)
                  call ga_add(1.0d0, g_rhs(j,i), 1.0d0, g_t, g_rhs(j,i))
               enddo
            endif
         enddo
c
         if (.not. ga_destroy(g_t))
     &       call errquit('j_deriv_gen: could not destroy temp GA',0)
c
      endif
c
c !!! BGJ test !!!
      if (bgj_print() .gt. 0) then
         if (do_hess) then
c
            write(6,*) 'j_deriv_gen: j hessian'
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
                  do j = 1, 3
                     write(*,*)'*** j cpks rhs matrix',j,i
                     call ga_print(g_rhs(j,i))
                  enddo
               endif
            enddo
c
c     TI check
c
            g_ti(1) = ga_create_atom_blocked (geom, ao_bas_han,
     &           'TI check x')
            g_ti(2) = ga_create_atom_blocked (geom, ao_bas_han,
     &           'TI check y')
            g_ti(3) = ga_create_atom_blocked (geom, ao_bas_han,
     &           'TI check z')
            do j = 1, 3
               call ga_zero(g_ti(j))
            enddo
            do i = 1, natoms
               if (oactive(i)) then
                  do j = 1, 3
                     call ga_add(1d0, g_rhs(j,i), 1d0, g_ti(j), g_ti(j))
                  enddo
               endif
            enddo
            do j = 1, 3
               write(*,*)'*** TI check',j
               call ga_print(g_ti(j))
               if (.not.ga_destroy(g_ti(j))) then
                  call errquit('j_deriv_gen: could not destroy TI',j)
               endif
            enddo
c
         endif
         write(*,*)'--------Leaving j_deriv_gen-------------'
      endif
c !!! BGJ test !!!
c
      return
      end
c
c
      subroutine j_nucder_alloc(l_d_ija, k_d_ija, l_d_kla, k_d_kla,
     &                          ldim, l_scr, k_scr, lscr, l_buf, k_buf,
     &                          lbuf, l_s34, k_s34, nat)
c
c     Allocates temporary storage needed by J nuclear derivative routines
c
      implicit none
c
#include "mafdecls.fh"
#include "bas.fh"
#include "cdft.fh"
c
      integer l_d_ija, k_d_ija, l_d_kla, k_d_kla, ldim, l_scr, k_scr,
     &        lscr, l_buf, k_buf, lbuf, l_s34, k_s34
      integer nat
c
c     Local variables
c
      integer max_at_bf, i, ilo, ihi
      integer max1e, max2e, mscratch_1e, mscratch_2e
c
c     Some allocations are the same for all calculation types
c
      max_at_bf = 0
      do i = 1, nat
         if (.not. bas_ce2bfr(ao_bas_han, i, ilo, ihi))
     &        call errquit('j_nucder_alloc: bas_ce2bfr failed', i)
         max_at_bf = max(max_at_bf, ihi-ilo+1)
      enddo
      ldim = max_at_bf
      if (.not. MA_push_get(MT_DBL, ldim*ldim, 'd_ija', l_d_ija,
     &     k_d_ija))
     &     call errquit('j_nucder_alloc: could not alloc d_ija',0)
      if (.not. MA_push_get(MT_DBL, ldim*ldim, 'd_kla', l_d_kla,
     &     k_d_kla))
     &     call errquit('j_nucder_alloc: could not alloc d_kla',0)
c
      if (.not. MA_push_get(MT_DBL, nat, 's34', l_s34, k_s34))
     &     call errquit('j_nucder_alloc: could not alloc s34',0)
c
c     determine memory requirements for integral derivatives
c
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      lbuf = max(max1e, max2e)
      lscr = max(mscratch_1e, mscratch_2e)
      if (.not. MA_push_get(MT_DBL, lbuf, 'deriv buffer', l_buf, k_buf))
     &     call errquit('j_nucder_alloc: could not allocate buffer',1)
      if (.not. MA_push_get(MT_DBL, lscr, 'deriv scratch', l_scr,
     &     k_scr))
     &     call errquit('j_nucder_alloc: could not allocate scratch',1)
c
      return
      end
c
c
      subroutine j_nucder_dealloc(l_d_ija, l_d_kla, l_scr, l_buf, l_s34)
c
c     Deallocates temporary storage used by J nuclear derivative routines
c
      implicit none
c
#include "mafdecls.fh"
c
      integer l_d_ija, l_d_kla, l_scr, l_buf, l_s34
c
      if (.not. MA_pop_stack(l_scr))
     &     call errquit('j_nucder_dealloc: could not pop scratch',1)
      if (.not. MA_pop_stack(l_buf))
     &     call errquit('j_nucder_dealloc: could not pop buffer',1)
      if (.not. MA_pop_stack(l_s34))
     &     call errquit('j_nucder_dealloc: could not pop s34',0)
      if (.not. MA_pop_stack(l_d_kla))
     &     call errquit('j_nucder_dealloc: could not pop d_kla',0)
      if (.not. MA_pop_stack(l_d_ija))
     &     call errquit('j_nucder_dealloc: could not pop d_ija',0)
c
      return
      end
