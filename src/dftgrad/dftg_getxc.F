      Subroutine dftg_getxc(rtdb, natom, iga_dens, force, oactive, 
     &                      nactive, s34, rdens_atom, cetobfr)
*
* $Id: dftg_getxc.F,v 1.22 2000-06-09 19:11:55 edo Exp $
*
      implicit none
c      
c     compute XC component of energy gradient.  
c
      integer rtdb
      integer natom
      integer iga_dens(2)
      double precision force(3*natom), s34(natom)
      integer nactive
      logical oactive(natom)
c!!! BGJ test !!!
      integer l_hess, k_hess, g_rhs(100), j
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c!!! BGJ test !!!
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "sym.fh"
#include "stdio.fh"
#include "cdft.fh"
c!!! BGJ test !!!
#include "bgj.fh"
c!!! BGJ test !!!
#include "cgridfile.fh"
c 
      integer cetobfr(2,natom)
      double precision rdens_atom(ipol*natom*natom)    
      double precision zero, one, onem
      parameter(zero = 0.0d0, one = 1.0d0, onem = -1.0d0)
      logical status,grid_on_file
c
      integer nproc,i
      double precision  tot
c
      nproc = ga_nnodes()
c
c     Exact exchange is now computed by grad_force
c

      call ga_sync
      tot = -xfac(1)
      do i = 1, numfunc
         tot = tot + xfac(i) + cfac(i)
      enddo
c
      if (abs(tot).gt.1e-8)then
c     
c        Compute via direct numerical quadrature.
c     
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &     grid_on_file))then
         grid_on_file = .false.
      endif
      if(grid_on_file.and.grid_written) then
         call dftg_gridv0(rtdb, iga_dens, force, natom, oactive, 
     &                    nactive, rdens_atom,  cetobfr)
      else
         call dftg_quadv0(rtdb, iga_dens, force, natom, oactive, 
     &                    nactive, rdens_atom, cetobfr)
      endif

c!!! BGJ test !!!
c
c     XC hesssian test calculation done by setting bgj:xc_hessian
c     to true
c
         if (.not. rtdb_get(rtdb, 'bgj:xc_hessian', mt_log,
     &        1, status)) status = .false.
         if (status) then
            status = MA_push_get(MT_DBL, 3*natom*3*natom,
     &           'xc hessian', l_hess, k_hess)
            if (.not.status) 
     &         call errquit('dftg_getxc: could not allocate xc hessian',
     &                      1)
            call dfill(9*natom*natom, 0.0d0, dbl_mb(k_hess), 1)
            if (bgj_print() .gt. 0)
     &           write(LuOut,*)'*** In dftg_getxc: calling xc_hessian'
            call xc_hessian(iga_dens, oactive, nactive, dbl_mb(k_hess))
            status = MA_pop_stack(l_hess)
            if (.not.status) call
     &         errquit('dftg_get_xc: could not pop stack for hessian',
     &                   1)
         endif

c
c     XC CPKS RHS test calculation done by setting bgj:xc_cpks_rhs
c     to true
c
         if (.not. rtdb_get(rtdb, 'bgj:xc_cpks_rhs', mt_log,
     &        1, status)) status = .false.
         if (status) then
c
c           Allocate and initialize temp GA's for RHS
c
            if (bgj_print() .gt. 0)
     &           write(*,*)'*** xc cpks rhs test: nactive =',nactive
            if (ipol*3*natom.gt.100)
     &           call errquit('dftg_getxc: dimension error in test',0)
            do i = 1, natom
             if (oactive(i)) then
               j = ipol*3*(i-1)
               j = j + 1
               g_rhs(j) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test a x')
               call ga_zero(g_rhs(j))
               if (ipol.eq.2) then
                  j = j + 1
                  g_rhs(j) = ga_create_atom_blocked
     &                 (geom, ao_bas_han, 'CPKS RHS test b x')
                  call ga_zero(g_rhs(j))
               endif
               j = j + 1
               g_rhs(j) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test a y')
               call ga_zero(g_rhs(j))
               if (ipol.eq.2) then
                  j = j + 1
                  g_rhs(j) = ga_create_atom_blocked
     &                 (geom, ao_bas_han, 'CPKS RHS test b y')
                  call ga_zero(g_rhs(j))
               endif
               j = j + 1
               g_rhs(j) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test a z')
               call ga_zero(g_rhs(j))
               if (ipol.eq.2) then
                  j = j + 1
                  g_rhs(j) = ga_create_atom_blocked
     &                 (geom, ao_bas_han, 'CPKS RHS test b z')
                  call ga_zero(g_rhs(j))
               endif
             endif
            enddo

            if (bgj_print() .gt. 0)
     &           write(LuOut,*)'*** In dftg_getxc: calling xc_cpks_rhs'
            call xc_cpks_rhs(iga_dens, oactive, nactive, g_rhs)

            do i = 1, natom
              if (oactive(i)) then
                do j = 1, ipol*3
                  if (.not.ga_destroy(g_rhs(j+ipol*3*(i-1)))) then
                    call errquit('xc_cpks_rhs: problem destroying ga',1)
                  endif
                enddo
              endif
            enddo
         endif
c!!! BGJ test !!!
      endif
c
      return
      end
