      Subroutine dftg_getxc(rtdb, natom, iga_dens, force, oactive, 
     &                      nactive, s34, rdens_atom, cntoce,
     &                      cntobfr, cetobfr)
*
* $Id: dftg_getxc.F,v 1.17 1999-06-17 18:25:10 edo Exp $
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
#include "cdft.fh"
c!!! BGJ test !!!
#include "bgj.fh"
c!!! BGJ test !!!
c 
      integer cntoce(nshells_ao), cntobfr(2,nshells_ao), 
     &        cetobfr(2,natom)
      double precision rdens_atom(ipol*natom*natom)    
      double precision zero, one, onem
      parameter(zero = 0.0d0, one = 1.0d0, onem = -1.0d0)
      logical status
c
      integer nproc, max_at_bf, lbuf, max1e, max2e,
     &        lscratch, mscratch_1e, mscratch_2e, i, lsqa
      integer l_buf, k_buf, l_scr, k_scr, l_ik, i_ik, l_il, i_il,
     &        l_jk, i_jk, l_jl, i_jl, ilo, ihi
      integer ispin
      double precision tol2e, tot
c
      nproc = ga_nnodes()
c

c
c     Exact exchange is now computed by grad_force
c

c$$$      if (abs(xfac(1)).gt.1e-8)then
c$$$c       
c$$$c        Compute the exact exchange potential (as in Hartree-Fock 
c$$$c        calculations).
c$$$c       
c$$$         tol2e=10.d0**(-itol2e)
c$$$c
c$$$         max_at_bf = 0
c$$$         do i = 1, natom
c$$$            if (.not. bas_ce2bfr(ao_bas_han, i, ilo, ihi))
c$$$     &         call errquit('dftg_getxc: bas_ce2bfr failed', i)
c$$$            max_at_bf = max(max_at_bf, ihi-ilo+1)
c$$$         enddo
c$$$c
c$$$         call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
c$$$         lbuf = max(max1e, max2e)
c$$$         lscratch = max(mscratch_1e, mscratch_2e)
c$$$         status = MA_Push_get(MT_DBL, lbuf, 'deriv buffer', l_buf,
c$$$     &            k_buf)
c$$$         if (.not.status)call errquit('could not allocate buffer',1)
c$$$         status = MA_Push_get(MT_DBL, lscratch, 'deriv scratch', 
c$$$     &            l_scr, k_scr)
c$$$         if (.not.status)call errquit('could not allocate deriv buf',1)
c$$$c
c$$$c        allocate local density matrix block
c$$$c
c$$$         lsqa = max_at_bf * max_at_bf
c$$$c
c$$$         status = status.and.MA_Push_get(MT_DBL,lsqa,
c$$$     &                       'local_w_density', l_ik, i_ik)
c$$$         status = status.and.MA_Push_get(MT_DBL,lsqa,
c$$$     &                       'local_w_density', l_il, i_il)
c$$$         status = status.and.MA_Push_get(MT_DBL,lsqa,
c$$$     &                       'local_w_density', l_jk, i_jk)
c$$$         status = status.and.MA_Push_get(MT_DBL,lsqa,
c$$$     &                       'local_w_density', l_jl, i_jl)
c$$$         if (.not.status) 
c$$$     &      call errquit('force_dft: could not allocate scratch',1)
c$$$c
c$$$         do ispin=1,ipol
c$$$         call dftg_hf_x(xfac(1), dbl_mb(i_ik), dbl_mb(i_jl),
c$$$     &                  dbl_mb(i_il), dbl_mb(i_jk), force, 
c$$$     &        iga_dens(ispin),
c$$$     &                  max_at_bf, AO_bas_han, nproc, natom, 
c$$$     &                  dbl_mb(k_scr), lscratch, dbl_mb(k_buf), lbuf, 
c$$$     &                  tol2e, oactive, nactive, geom, s34, oskel,
c$$$     &          ipol)
c$$$         enddo
c$$$c
c$$$         status = MA_Pop_Stack (l_jl)
c$$$         status = MA_Pop_Stack (l_jk)
c$$$         status = MA_Pop_Stack (l_il)
c$$$         status = MA_Pop_Stack (l_ik)
c$$$         status = MA_Pop_Stack (l_scr)
c$$$         status = MA_Pop_Stack (l_buf)
c$$$         if (.not.status) 
c$$$     &      call errquit('force_dft: could not pop stack',1)
c$$$c 
c$$$         if (.not.status) 
c$$$     &      call errquit(' dftg_getcx: failed to deallocate',11)
c$$$c
c$$$      endif
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
         call dftg_quadv0(rtdb, iga_dens, force, natom, oactive, 
     &                    nactive, rdens_atom, cntoce,
     &                    cntobfr, cetobfr)

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
            if (bgj_print() .gt. 0)
     &           write(6,*)'*** In dftg_getxc: calling xc_hessian'
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
c           Allocate temp GA's for RHS
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
               if (ipol.eq.2) then
                  j = j + 1
                  g_rhs(j) = ga_create_atom_blocked
     &                 (geom, ao_bas_han, 'CPKS RHS test b x')
               endif
               j = j + 1
               g_rhs(j) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test a y')
               if (ipol.eq.2) then
                  j = j + 1
                  g_rhs(j) = ga_create_atom_blocked
     &                 (geom, ao_bas_han, 'CPKS RHS test b y')
               endif
               j = j + 1
               g_rhs(j) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test a z')
               if (ipol.eq.2) then
                  j = j + 1
                  g_rhs(j) = ga_create_atom_blocked
     &                 (geom, ao_bas_han, 'CPKS RHS test b z')
               endif
             endif
            enddo

            if (bgj_print() .gt. 0)
     &           write(6,*)'*** In dftg_getxc: calling xc_cpks_rhs'
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
