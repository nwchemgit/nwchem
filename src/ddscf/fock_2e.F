      subroutine fock_2e( geom1, ao_basis, nfock, jfac, kfac,
     $                    tol2e1, oskel1, vg_dens, vg_fock )
c$Id: fock_2e.F,v 1.12 1999-06-22 19:22:15 edo Exp $
c
c     wrapper routine for ao_fock_2e and riscf_fock_2e
c
c
      implicit none
c     !!! BGJ test !!!
#include "global.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c     !!! BGJ test !!!
#include "cscf.fh"
#include "cfock.fh"

      integer geom1, ao_basis        ! [input] parameter handles
      integer nfock                  ! [input] number of Fock matrices
      double precision jfac(nfock)   ! [input] Coulomb prefactor
      double precision kfac(nfock)   ! [input] exchange prefactor
      double precision tol2e1        ! [input] integral selection threshold
      logical oskel1                 ! [input] toggle skeleton Fock matrix
      integer vg_dens(nfock)         ! [input] array of handles to densities
      integer vg_fock(nfock)         ! [input] array of handles to Fock matrices
c     !!! BGJ test !!!
      integer ifock
      integer g_xc(100) ! Assume no more than 100 fock matrices !!!
      integer ga_create_atom_blocked
      integer vg_fockc
      external ga_create_atom_blocked
      logical xc_active, jfit
c     !!! BGJ test !!!

c     !!! BGJ test !!!
c
c     See if XC is active for this call
c
      if (.not.rtdb_get(bgj_get_rtdb_handle(),'bgj:xc_active',MT_LOG,1,
     &     xc_active)) xc_active = .false.
      if (bgj_print().gt.0)
     &   write(*,*)'*** fock_2e: xc_active ',xc_active

c     Set K contribution to 0 for DFT - not good enough in general,
c     e.g. for hybrids, but will do for testing pure DFT.
c
      if (xc_active .and. bgj_have_xc()) then
         if (bgj_print().gt.0)
     &      write(*,*)'*** fock_2e: zeroing kfac ',kfac
         call dfill(nfock, 0.0d0, kfac, 1)
      endif
c
c     Determine whether J fitting is involved
c
c      jfit = bgj_have_j_fit().and.(jfac(1).ne.0d0)
      jfit=.false.

c     !!! BGJ test !!!
      if (jfit) then
c
c     get the HF exchange part (must be replaced by DFT X)
c
         vg_fockc = ga_create_atom_blocked
     $        (geom1, ao_basis, 'fock_2e: vg_fockc')
         jfac(1)=0d0
         kfac(1)=-0.5d0
         call ao_fock_2e( geom1, ao_basis, 1, jfac, kfac,
     $                    tol2e1, oskel1, vg_dens, vg_fock )

         call fock_j_fit(nfock, vg_dens, vg_fockc)
         call ga_add(1.0d0, vg_fockc, 1.0d0, vg_fock,
     &        vg_fock)
         if (.not.ga_destroy(vg_fockc))
     .      call errquit('fock_2e: problem destroying test array',1)
         return
      else if (nriscf.gt.1) then
        call riscf_fock_2e( geom1, ao_basis, riscf_basis, nfock, 
     $       jfac, kfac, tol2e1, oskel1, vg_dens, vg_fock )
      else if (nriscf.eq.1) then
        if (rifock) then
          call riscf_fock_2e( geom1, ao_basis, riscf_basis, nfock, 
     $         jfac, kfac, tol2e1, oskel1, vg_dens, vg_fock )
        else
          call ao_fock_2e( geom1, ao_basis, nfock, jfac, kfac,
     $         tol2e1, oskel1, vg_dens, vg_fock )
        endif
      else
          call ao_fock_2e( geom1, ao_basis, nfock, jfac, kfac,
     $                    tol2e1, oskel1, vg_dens, vg_fock )
      endif

c     !!! BGJ test !!!

c !!! Exit if no XC part to do - this can be changed to just do XC
c !!! part and print it out for comparison by removing xc_active
c !!! from the test below
      if (.not. (xc_active .and. bgj_have_xc())) then
c         write(*,*)'*** fock_2e: no xc, returning'
         return
      endif

c      write(*,*)'*** fock_2e: nfock =',nfock
      if (nfock.gt.100)
     &   call errquit('fock_2e: dimension error in test',0)

      do ifock = 1, nfock
         g_xc(ifock) = ga_create_atom_blocked
     $        (geom1, ao_basis, 'fock_2e:BGJ XC test')
      enddo

      call fock_xc(nfock, vg_dens, g_xc)

      do ifock = 1, nfock
         if (xc_active) then
            call ga_add(1.0d0, g_xc(ifock), 1.0d0, vg_fock(ifock),
     &           vg_fock(ifock))
         endif
         if (bgj_print() .gt. 1) then
            write(*,*)'*** fock_2e: fock matrix',ifock
            call ga_print(vg_fock(ifock))
         endif
         if (.not.ga_destroy(g_xc(ifock))) then
            call errquit('fock_2e: problem destroying test array',1)
         endif
      enddo

c     !!! BGJ test !!!

      return
      end
