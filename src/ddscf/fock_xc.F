      subroutine fock_xc(geom, nbf_ao, ao_bas_han,
     ,     nmat, g_dens, g_xc,scftype)
c
c     $Id: fock_xc.F,v 1.14 2001-11-16 02:08:33 edo Exp $
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
c!!! BGJ
c
      integer geom              ! [in]
      integer nbf_ao,ao_bas_han ! [in]
      character*4 scftype       ! [in]
      integer nmat              ! [input] number of XC matrices
      integer g_dens(nmat)      ! [input] array of handles to densities
      integer g_xc              ! [input] array of handles to XC matrices
      integer calc_type         ! [input] type of XC matrix calculation
c                               !  *** from rtdb
c                               !  ******* !!! document here !!!
c
c     Local declarations
c
      integer imat
      integer rtdb
      integer natoms
      integer nExc_dum, g_wght_dum, g_xyz_dum, g_nq_dum,
     &        l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &        l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom,
     &        l_g_dens_local, k_g_dens_local, k_g_dens_scf
      integer ipol
      logical wght_GA_dum,grid_reopen,grid_on_file,xc_gotxc
      external grid_reopen,xc_gotxc
      double precision Exc_dum(2), rho_n_dum,tdum
      logical cphf_tddft
      integer ioff, joff
c
      if(.not.xc_gotxc()) return
      ipol=1
      if(scftype.eq.'UHF') ipol=2

      if (bgj_print() .gt. 0)
     &     write(*,*)'--------Entered fock_xc-------------'
      rtdb = bgj_get_rtdb_handle()
c      write(*,*)'nmat rtdb geom ipol',nmat,rtdb,geom,ipol
      if (.not. geom_ncent(geom, natoms))
     &   call errquit('fock_xc: geom_ncent failed',1)
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &     grid_on_file))then
         grid_on_file = .false.
      endif
c      write(*,*)'natoms',natoms
c
c     Get fock_xc variables
c
      if (.not. rtdb_get(rtdb, 'fock_xc:calc_type', mt_int, 1,
     &     calc_type)) then
         if (bgj_print() .gt. 0)
     &        write(*,*)'*** fock_xc: calc_type not set: setting to 1'
         calc_type = 1          !!!
      endif
      if (calc_type .eq. 0) then
         if (bgj_print() .gt. 0)
     &        write(*,*)'*** fock_xc: calc_type = 0: setting to 1'
         calc_type = 1          !!!
      endif
      if (bgj_print() .gt. 0)
     &     write(*,*)'fock_xc: calc_type =',calc_type
c
c     Set up local copies of density matrix handles. The number of 
c     density matrices needed is ipol*nvec+2. The last two are scf
c     (unperturbed) density matrices. nmat is actually 
c
      if (.not.MA_Push_Get(mt_int,nmat+2,'g_dens_local',
     &     l_g_dens_local,k_g_dens_local))
     &     call errquit('fock_xc: cannot allocate g_dens_local',0)
      if(ipol.eq.1)then
         do imat = 1, nmat
            int_mb(k_g_dens_local+imat-1) = g_dens(imat)
         enddo
      else if(ipol.eq.2)then 
         do imat = 1, nmat/ipol
            ioff = k_g_dens_local + 2*(imat-1) - 1
            joff = 4*(imat-1)
            int_mb(ioff + 1) = g_dens(joff+1) ! alpha density
            int_mb(ioff + 2) = g_dens(joff+3) ! beta  density
         enddo
      else
         call errquit("unknown ipol, fock_xc", 0)
      endif
      if (calc_type .eq. 2) then
         if(ipol.eq.1)then
            k_g_dens_scf = k_g_dens_local+nmat
         else if(ipol.eq.2)then
            k_g_dens_scf = k_g_dens_local+nmat/2 
         endif
         if (.not. bgj_get_scf_dens(int_mb(k_g_dens_scf)))
     &      call errquit('fock_xc: could not get scf density handles',0)
      else
         k_g_dens_scf = k_g_dens_local
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
      call xc_init_index(geom,ipol,nbf_ao,ao_bas_han,
     ,     int_mb(k_g_dens_scf), l_cntoce, k_cntoce,
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
      if(cphf_tddft .and. scftype.eq.'UHF')nmat=nmat/4
c     
c     check if we can get the grid from a file
c         
      grid_written=grid_reopen(rtdb,geom)
      if(grid_on_file.and.grid_written) then
      call grid_quadv0_gen(rtdb, 
     ,         int_mb(k_g_dens_local),g_xc,nexc_dum,
     &     rho_n_dum, Exc_dum,
     &     nmat, calc_type ,tdum,.true.)
      else
      call xc_quadv0_gen(rtdb, nExc_dum, g_xc, int_mb(k_g_dens_local),
     &     g_wght_dum, g_xyz_dum, g_nq_dum, wght_GA_dum, Exc_dum,
     &     rho_n_dum, dbl_mb(k_rdens_atom), 
     &     int_mb(k_cetobfr), natoms, nmat,
     &     calc_type)
      endif
c
c     Clean up
c
      call xc_exit_index(l_cntoce, l_cntobfr, l_cetobfr, l_rdens_atom)
c
      if (calc_type .eq. 2) then
         if (.not. ga_destroy(int_mb(k_g_dens_scf)))
     &      call errquit('fock_xc: could not destroy scf density',0)
         if (ipol. eq. 2) then
            if (.not. ga_destroy(int_mb(k_g_dens_scf+1)))
     &         call errquit('fock_xc: could not destroy scf density',0)
         endif
      endif
c
      if (.not.ma_pop_stack(l_g_dens_local))
     &   call errquit('fock_xc: cannot pop stack',0)

      call grid_cleanup(rtdb, .true.)
c
      if (bgj_print() .gt. 0) then
         write(*,*)'XC matrix',imat,' at end of fock_xc:'
         call ga_print(g_xc)
         write(*,*)'--------Leaving fock_xc-------------'
      endif
      return
      end
