      Subroutine dftg_getxc(rtdb, natom, iga_dens, force, oactive, 
     &                      nactive,  rdens_atom, cetobfr, cd_coef)
*
* $Id$
*
      implicit none
c      
c     compute XC component of energy gradient.  
c
      integer rtdb
      integer natom
      integer iga_dens(2)
      double precision force(3*natom)
      integer nactive
      logical oactive(natom)
c!!! BGJ test !!!
      integer l_hess, k_hess, g_rhs(100), j
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c!!! BGJ test !!!

      integer icetobfr,lcetobfr,icntoce,lcntoce,icntobfr,lcntobfr
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "sym.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "errquit.fh"
c!!! BGJ test !!!
#include "bgj.fh"
c!!! BGJ test !!!
#include "cgridfile.fh"
#include "geom.fh"
c 
      integer cetobfr(2,natom)
      double precision rdens_atom(ipol*natom*natom)    
      double precision cd_coef(*)
      double precision zero, one, onem
      parameter(zero = 0.0d0, one = 1.0d0, onem = -1.0d0)
      logical status,xc_gotxc
      external xc_gotxc
c
      character*4 scftype
c
      integer nproc,i
      double precision tol_rho_scf ! store scf tol_rho before hack
      logical grid_reopen
      external grid_reopen

      integer nbf_bas,basis
c
c     hack for hessian
c
      scftype="RHF "
      if(ipol.eq.2) scftype="UHF "
      nproc = ga_nnodes()
c
c     Exact exchange is now computed by grad_force
c

      call ga_sync
c
      if (xc_gotxc())then
c     
c     hack tol_rho to improve gradient precision
c
      tol_rho_scf=tol_rho
      tol_rho=min(tol_rho*0.1d0,1d-11)


      if (adft) then
        nbf_bas = nbf_cd
        basis = cd_bas_han
        if (.not.ma_push_get
     $       (mt_int,natom*2,'cntoce map',lcetobfr,icetobfr))
     $       call errquit('dft_scf:push_get failed', 13, MA_ERR)
        if (.not.ma_push_get
     $       (mt_int,nshells_cd,'cntoce map',lcntoce,icntoce))
     $       call errquit('dft_scf:push_get failed', 13, MA_ERR)
        if (.not.ma_push_get
     $       (mt_int,nshells_cd*2,'cntoce map',lcntobfr,icntobfr))
     $       call errquit('dft_scf:push_get failed', 13, MA_ERR)
        call icopy(natom*2,cetobfr,1,int_mb(icetobfr),1)
        call build_maps(cd_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     $                  cetobfr, natom, nshells_cd)
        if (.not.ma_chop_stack(lcntoce))
     $     call errquit('dft_gradient: cannot pop stack',0, MA_ERR)
      else
        nbf_bas = nbf_ao
        basis = ao_bas_han
      endif
c
c        Compute via direct numerical quadrature.
         if(.not.grid_written) then
c
c        check if we can get the grid from a file
c         
            grid_written=grid_reopen(geom)
         endif
         if(grid_written) then
            call grid_setquad(rtdb)
            call grid_packing_info
         endif
         call dftg_gridv0(rtdb, geom, basis, ipol, nbf_bas,
     .        iga_dens, force, natom, oactive, 
     &        nactive, rdens_atom,  cetobfr, cd_coef, adft)
         oreopen=.false.

         if (adft) then
            call icopy(2*natom,int_mb(icetobfr),1,cetobfr,1)
            if (.not.ma_chop_stack(lcetobfr))
     $       call errquit('dft_gradient: cannot pop stack',0, MA_ERR)
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
     &                      1, MA_ERR)
            call dfill(9*natom*natom, 0.0d0, dbl_mb(k_hess), 1)
            if (bgj_print() .gt. 0)
     &           write(LuOut,*)'*** In dftg_getxc: calling xc_hessian'
            call xc_hessian(geom,rtdb,nbf_ao,ao_bas_han,
     .           iga_dens, oactive, nactive, dbl_mb(k_hess),
     .           scftype)
            status = MA_pop_stack(l_hess)
            if (.not.status) call
     &         errquit('dftg_get_xc: could not pop stack for hessian',
     &                   1, MA_ERR)
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
     &           call errquit('dftg_getxc: dimension error in test',0,
     &       INPUT_ERR)
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
            call xc_cpks_rhs(geom,rtdb,nbf_ao,ao_bas_han,
     .           iga_dens, oactive, nactive, g_rhs,
     .           scftype)

            do i = 1, natom
              if (oactive(i)) then
                do j = 1, ipol*3
                  if (.not.ga_destroy(g_rhs(j+ipol*3*(i-1)))) then
                    call errquit('xc_cpks_rhs: problem destroying ga',1,
     &       GA_ERR)
                  endif
                enddo
              endif
            enddo
         endif
c!!! BGJ test !!!
      endif
c
c    restore tol_rho value
c
      tol_rho = tol_rho_scf
c
      return
      end
