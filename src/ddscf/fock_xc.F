      subroutine fock_xc(geom, nbf_ao, ao_bas_han,
     ,  nfock, g_dens, g_xc,l3d)
c
c     $Id: fock_xc.F,v 1.20 2002-02-19 02:10:37 edo Exp $
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
      integer nExc_dum, g_wght_dum, g_xyz_dum, g_nq_dum,
     &        l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &        l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom
      integer g_dens_scf(2)
      logical wght_GA_dum,grid_reopen,xc_gotxc
      double precision Exc_dum(2), rho_n_dum,tdum
      integer dims(3),chunk(3)
      integer g_dens_xc(2),xc_getipol
      logical oprint
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer nfock_xc,i
      integer ityp, ndim
      external grid_reopen,xc_gotxc,xc_getipol
c
      if(.not.xc_gotxc()) return
      oprint= util_print('fock_xc',print_debug)
      ndens=1
      if(xc_getipol().eq.2) ndens=2
      rtdb = bgj_get_rtdb_handle()
      alo(2) = 1
      ahi(2) = nbf_ao
      alo(3) = 1
      ahi(3) = nbf_ao
      blo(2) = 1
      bhi(2) = nbf_ao
      blo(3) = 1
      bhi(3) = nbf_ao
c
c     Get fock_xc variables
c
      if (.not. rtdb_get(rtdb, 'fock_xc:calc_type', mt_int, 1,
     &     calc_type).or.calc_type .eq. 0) then
         if (oprint.and.ga_nodeid().eq.0)
     &        write(luout,*)' fockxc: calc_type not set: setting to 1'
         calc_type = 1          !!!
      endif
c
c     for uhf calculations, nmat is 4*nvec. The nmat passed to _quadv0_gen
c     should be the number of vectors 
c
      nfock_xc=nfock
      if(ndens.eq.2)  then
         nfock_xc=nfock/4
      endif
      if (oprint.and.ga_nodeid().eq.0) then
         write(luout,*) '--------Entered fock_xc-------------'
         write(luout,*) ' calc_type =',calc_type
         write(luout,*) ' nfock ',nfock,' ndens ',ndens
         write(luout,*) ' nfock_xc ',nfock_xc
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
c     g_dens_xc is the DM used in the XC routines
c
      if(l3d) then
         if(ndens.eq.1) then
            dims(1) = nfock_xc+1
         else
            dims(1) = 2*nfock_xc+2
         endif
         dims(2) = nbf_ao
         dims(3) = nbf_ao
         chunk(1) = dims(1)
         chunk(2) = -1
         chunk(3) = -1
         if (.not. nga_create (MT_DBL, 3, dims, 'XCDensity',chunk,
     &     g_dens_xc)) call errquit
     &     ('fock_xc: could not allocate g_dens_xc',555)
c
c     alpha bit
c
         alo(1) = 1
         ahi(1) = nfock_xc
         blo(1) = 1
         bhi(1) = nfock_xc
         call ga_zero(g_dens_xc)
         call nga_add_patch(1d0,g_dens,alo,ahi,
     .           0d0,g_dens_xc,blo,bhi,
     .        g_dens_xc,blo,bhi)      
c
c     beta bit
c
         if(ndens.eq.2)  then
            alo(1) = nfock_xc*2+1
            ahi(1) = nfock_xc*3
            blo(1) = nfock_xc+1+1
            bhi(1) = nfock_xc*2+1
            call nga_add_patch(1d0,g_dens,alo,ahi,
     .           0d0,g_dens_xc,blo,bhi,
     .           g_dens_xc,blo,bhi)      
         endif
      else
         g_dens_xc(1)=g_dens(1)
c
c     pick 3rd
c
         if(ndens.eq.2) g_dens_xc(2)=g_dens(3)
         if(oprint.and.ga_nodeid().eq.0) then
            write(luout,*) ' picked 3rd dens ',g_dens_xc
         endif
      endif
      if (calc_type .eq. 2) then
         if (.not. bgj_get_scf_dens(g_dens_scf))
     &        call errquit('fock_xc: cant get scf density handles',0)
         alo(1) = 1
         ahi(1) = nbf_ao
         blo(1) = nfock_xc+1
         bhi(1) = nfock_xc+1
         call nga_copy_patch('N',g_dens_scf(1),alo,ahi,
     .        g_dens_xc,blo,bhi)
         if(ndens.eq.2) then
            blo(1) = 2*nfock_xc+1+1
            bhi(1) = 2*nfock_xc+1+1
            call nga_copy_patch('N',g_dens_scf(2),alo,ahi,
     .           g_dens_xc,blo,bhi)
         endif
      else
         g_dens_scf(1) = g_dens_xc(1)
         if(ndens.eq.2) g_dens_scf(2) = g_dens_xc(2)
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
      if(oprint) then
         if(ga_nodeid().eq.0) then
            write(luout,*) 'gdenscf1',g_dens_scf(1)
            call util_flush(luout)
         endif
         call ga_print(g_dens_scf(1))
         if(ndens.eq.2) then
            if(ga_nodeid().eq.0) then
               write(luout,*) 'gdenscf2',g_dens_scf(2)
               call util_flush(luout)
            endif
            call ga_print(g_dens_scf(2))
         endif
         if(ga_nodeid().eq.0) then
            write(luout,*) 'g_dens_xc1',g_dens_xc(1)
            call util_flush(luout)
         endif
         if(l3d) then
            do i=1,ndens*(nfock_xc+1)
               alo(1)=i
               ahi(1)=i
               call nga_print_patch(g_dens_xc,alo,ahi,0)
            enddo
         else
            call ga_print(g_dens_xc(1))
         endif
         if(ndens.eq.2.and.(.not.l3d)) then
            if(ga_nodeid().eq.0) then
               write(luout,*) 'g_dens_xc2',g_dens_xc(2)
               call util_flush(luout)
            endif
            call ga_print(g_dens_xc(2))
         endif
      endif
      call xc_init_index(geom,ndens,nbf_ao,ao_bas_han,
     ,     g_dens_scf, l_cntoce, k_cntoce,
     &                   l_cntobfr, k_cntobfr, l_cetobfr, k_cetobfr,
     &                   l_rdens_atom, k_rdens_atom)
c     
c     check if we can get the grid from a file
c         
      grid_written=grid_reopen(rtdb,geom)
      if(grid_written) then
      call grid_quadv0_gen(rtdb, 
     ,         g_dens_xc,g_xc,nexc_dum,
     &     rho_n_dum, Exc_dum,
     &     nfock_xc, calc_type ,tdum,.true.)
      else
         call errquit( 'fxc: old grid broken ',0)
!      call xc_quadv0_gen(rtdb, nExc_dum, g_xc, int_mb(g_dens_local),
      endif
c
c     Clean up
c
      if(l3d) then
         if (.not. ga_destroy(g_dens_xc))
     &        call errquit('fock_xc: could not destroy g_dens_xc',0)
      endif
      if(ndens.eq.2)  then
c
c     need to go from 2nfock ga to 4nfock (and sign flip)
c
c     move beta from block 2 to block 4
            alo(1) = nfock_xc+1
            ahi(1) = nfock_xc*2
            blo(1) = nfock_xc*3+1
            bhi(1) = nfock_xc*4
            call nga_copy_patch('N',g_xc,alo,ahi,
     .           g_xc,blo,bhi)
c     move alpha from block 1 to block 2
            alo(1) = 1
            ahi(1) = nfock_xc
            blo(1) = nfock_xc+1
            bhi(1) = nfock_xc*2
            call nga_copy_patch('N',g_xc,alo,ahi,
     .           g_xc,blo,bhi)
c     zero block 1 that was occupied by alpha
            call nga_zero_patch(g_xc, alo, ahi)
c     flip sign to make it consistent with HF K (Exchange)
            call ga_scale(g_xc,-1d0)
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
         tdum=ga_ddot(g_xc,g_xc)
         if(ga_nodeid().eq.0)write(luout,*)'g_xc dotproduct',tdum
         do i=1,nfock
            alo(1)=i
            ahi(1)=i
            tdum=nga_ddot_patch(g_xc, 'N', alo, ahi, 
     .           g_xc, 'N', alo, ahi) 


            if(ga_nodeid().eq.0) then
               write(luout,*) ' matrix no. ',i,tdum
               call util_flush(luout)
            endif
            if(abs(tdum).gt.1d-6) then
            if(ga_nodeid().eq.0) then
               call nga_inquire(g_xc, ityp, ndim, dims)
               write(luout,*) ' ndim ',ndim,' dims ',dims
               write(luout,*) ' alo ',alo
               write(luout,*) ' ahi ',ahi
               call util_flush(luout)
            endif
            call nga_print_patch(g_xc,alo,ahi,0)
            endif
         enddo
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
