      subroutine fock_xc(geom, nbf_ao, ao_bas_han,
     ,  nfock, g_dens, g_xc,scftype,l3d)
c
c     $Id: fock_xc.F,v 1.17 2002-01-16 01:36:56 edo Exp $
c
c     Wrapper routine for AO-basis XC matrices without fitting
c
c     BGJ - 8/98
c
      implicit none
c
#include "global.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "../nwdft/grid/cgridfile.fh"
c!!! BGJ
#include "bgj.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c!!! BGJ
c
      integer geom              ! [in]
      integer nbf_ao,ao_bas_han ! [in]
      character*4 scftype       ! [in]
      integer nfock             ! [input] number of XC matrices
      integer ndens             ! [in] number of DM
      integer g_dens(*)        ! [input] array of handles to densities
      integer g_xc              ! [input] array of handles to XC matrices
      integer calc_type         ! [input] type of XC matrix calculation
      logical l3d
c                               !  *** from rtdb
c                               !  ******* !!! document here !!!
c
c     Local declarations
c
      integer rtdb
      integer natoms
      integer nExc_dum, g_wght_dum, g_xyz_dum, g_nq_dum,
     &        l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &        l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom
      integer g_dens_scf(2)
      logical wght_GA_dum,grid_reopen,grid_on_file,xc_gotxc
      external grid_reopen,xc_gotxc
      double precision Exc_dum(2), rho_n_dum,tdum
      logical cphf_tddft
cold      integer ioff, joff
cold      integer imat
      integer dims(3),chunk(3)
      integer dft_npatch,g_dens2(2)
      logical oprint
      integer alo(3), ahi(3), blo(3), bhi(3)
      external dft_npatch
c
      if(.not.xc_gotxc()) return
      oprint= util_print('fock_xc',print_debug)
      ndens=1
      if(scftype.eq.'UHF') ndens=2
      rtdb = bgj_get_rtdb_handle()
      if (.not. geom_ncent(geom, natoms))
     &   call errquit('fock_xc: geom_ncent failed',1)
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &     grid_on_file))then
         grid_on_file = .false.
      endif
c
c     Get fock_xc variables
c
      if (.not. rtdb_get(rtdb, 'fock_xc:calc_type', mt_int, 1,
     &     calc_type)) then
         if (oprint.and.ga_nodeid().eq.0)
     &        write(luout,*)' fockxc: calc_type not set: setting to 1'
         calc_type = 1          !!!
      endif
      if (calc_type .eq. 0) then
         if (oprint.and.ga_nodeid().eq.0)
     &        write(luout,*)'*** fock_xc: calc_type = 0: setting to 1'
         calc_type = 1          !!!
      endif
      if (oprint.and.ga_nodeid().eq.0) then
         write(luout,*) '--------Entered fock_xc-------------'
         write(luout,*) ' nfock ',nfock,' ndens ',ndens
         write(luout,*) ' calc_type =',calc_type
         write(luout,*) ' l3d ',l3d
         call util_flush(luout)
       endif
c
c     assumed l3d true AND calc_type=2
c
      if(l3d.and.calc_type.ne.2.or.
     .      (.not.l3d).and.calc_type.eq.2) call errquit('
     '      fxc: calc_type-l3d logic wrong ',0)
c
c     Set up local copies of density matrix handles. The number of 
c     density matrices needed is ndens*nvec+2. The last two are scf
c     (unperturbed) density matrices. nmat is actually 
c
      if(l3d) then
         if(ndens.eq.1) then
            dims(1) = nfock+1
         else
            dims(1) = nfock+2
         endif
         dims(2) = nbf_ao
         dims(3) = nbf_ao
         chunk(1) = dims(1)
         chunk(2) = -1
         chunk(3) = -1
         if (.not. nga_create (MT_DBL, 3, dims, 'Density',chunk,
     &     g_dens2)) call errquit
     &     ('fock_xc: could not allocate g_dens2',555)
         alo(1) = 1
         ahi(1) = nfock
         alo(2) = 1
         ahi(2) = nbf_ao
         alo(3) = 1
         ahi(3) = nbf_ao
         blo(1) = 1
         bhi(1) = nfock
         blo(2) = 1
         bhi(2) = nbf_ao
         blo(3) = 1
         bhi(3) = nbf_ao
         call ga_zero(g_dens2)
         call nga_add_patch(1d0,g_dens,alo,ahi,
     .           0d0,g_dens2,blo,bhi,
     .        g_dens2,blo,bhi)      
      else
         g_dens2(1)=g_dens(1)
         if(ndens.eq.2) g_dens2(2)=g_dens(2)
      endif
cold     do imat = 1, nfock/ndens
cold     ioff = k_g_dens_local + 2*(imat-1) - 1
cold     joff = 4*(imat-1)
cold     int_mb(ioff + 1) = g_dens(joff+1) ! alpha density
cold     int_mb(ioff + 2) = g_dens(joff+3) ! beta  density
cold     enddo
      if (calc_type .eq. 2) then
cold        if(ndens.eq.1)then
cold            k_g_dens_scf = k_g_dens_local+nfock
cold        else if(ndens.eq.2)then
cold            k_g_dens_scf = k_g_dens_local+nfock/2 
cold        endif
         if (.not. bgj_get_scf_dens(g_dens_scf))
     &        call errquit('fock_xc: cant get scf density handles',0)
         write(6,*) 'denscf ',g_dens_scf
         call util_flush(6)
         alo(1) = 1
         ahi(1) = nbf_ao
         blo(1) = nfock+1
         bhi(1) = nfock+1
         call nga_copy_patch('N',g_dens_scf(1),alo,ahi,
     .        g_dens2,blo,bhi)
         if(ndens.eq.2) then
            blo(1) = nfock+2
            bhi(1) = nfock+2
            call nga_copy_patch('N',g_dens_scf(2),alo,ahi,
     .           g_dens2,blo,bhi)
         endif
      else
         g_dens_scf(1) = g_dens(1)
         if(ndens.eq.2) g_dens_scf(2) = g_dens(2)
      endif
c
c     Prepare to call xc_quadv0
c
      nExc_dum = 1
      g_wght_dum = -1
      g_xyz_dum = -1
      g_nq_dum = -1
      wght_GA_dum = .false.
      rho_n_dum = 0
      exc_dum(1) = 0
      exc_dum(2) = 0
      tdum = 0
      call xc_init_index(geom,ndens,nbf_ao,ao_bas_han,
     ,     g_dens_scf, l_cntoce, k_cntoce,
     &                   l_cntobfr, k_cntobfr, l_cetobfr, k_cetobfr,
     &                   l_rdens_atom, k_rdens_atom)
c
c     for uhf calculations, nmat is 4*nvec. The nmat passed to _quadv0_gen
c     should be the number of vectors 
c
      if (.not. rtdb_get(bgj_get_rtdb_handle(), 
     &     'cphf_solve:cphf_tddft', mt_log, 1, cphf_tddft)) then
         cphf_tddft = .false.
      endif
      if(cphf_tddft .and. ndens.eq.2)nfock=nfock/4
c     
c     check if we can get the grid from a file
c         
      grid_written=grid_reopen(rtdb,geom)
      if(grid_on_file.and.grid_written) then
      call grid_quadv0_gen(rtdb, 
     ,         g_dens2,g_xc,nexc_dum,
     &     rho_n_dum, Exc_dum,
     &     nfock, calc_type ,tdum,.true.)
      else
         call errquit( 'fxc: old grid broken ',0)
!      call xc_quadv0_gen(rtdb, nExc_dum, g_xc, int_mb(g_dens_local),
!     &     g_wght_dum, g_xyz_dum, g_nq_dum, wght_GA_dum, Exc_dum,
!     &     rho_n_dum, dbl_mb(k_rdens_atom), 
!     &     int_mb(k_cetobfr), natoms, nfock,
!     &     calc_type)
      endif
c
c     Clean up
c
      if(l3d) then
         if (.not. ga_destroy(g_dens2))
     &        call errquit('fock_xc: could not destroy g_dens2',0)
      endif
      call xc_exit_index(l_cntoce, l_cntobfr, l_cetobfr, l_rdens_atom)
c
      if (calc_type .eq. 2) then
         if (.not. ga_destroy(g_dens_scf(1)))
     &      call errquit('fock_xc: could not destroy DM',1)
         if (ndens. eq. 2) then
            if (.not. ga_destroy(g_dens_scf(2)))
     &           call errquit('fock_xc: could not destroy DM',2)
         endif
      endif
c

      call grid_cleanup(rtdb, .true.)
c
      if (oprint) then
         if(ga_nodeid().eq.0)write(luout,*)'XC matrix at end fock_xc:'
         call ga_print(g_xc)
         call util_flush(luout)
      endif
      return
      end
      integer function dft_npatch(g_dens,imat,nbf,geom,basis)
      implicit none
      integer g_dens ! N-D GA
      integer imat
      integer g_dens1 ! 2-D GA
      integer nbf
      integer geom,basis
c
      
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer alo(3), ahi(3), blo(2), bhi(2)
c
c     2-d GA dens1 
c
      g_dens1 = ga_create_atom_blocked(geom,basis,'gdens1')
      alo(1) = imat
      ahi(1) = imat
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf

      call nga_copy_patch('N',g_dens,alo,ahi,g_dens1,blo,bhi)
      dft_npatch=g_dens1
      return
      end
