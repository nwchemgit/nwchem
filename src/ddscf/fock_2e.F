      subroutine fock_2e( geom1, ao_basis, nfock, facj, fack,
     $                    tol2e1, oskel1, vg_dens, vg_fock, asym )
c$Id: fock_2e.F,v 1.16 2000-08-25 20:07:12 edo Exp $
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
      double precision facj(nfock)   ! [input] Coulomb prefactor
      double precision fack(nfock)   ! [input] exchange prefactor
      double precision tol2e1        ! [input] integral selection threshold
      logical oskel1                 ! [input] toggle skeleton Fock matrix
      integer vg_dens(nfock)         ! [input] array of handles to densities
      integer vg_fock(nfock)         ! [input] array of handles to Fock matrices
      logical asym                   ! [input] flag to antisymmetrize Fock matrices
c     !!! BGJ test !!!
      integer ifock
      integer maxfocks
      parameter(maxfocks=100)
      integer g_xc(maxfocks) ! Assume no more than 100 fock matrices !!!
      double precision kfac(maxfocks), jfac(maxfocks)
      integer ga_create_atom_blocked
      integer vg_fockc
      external ga_create_atom_blocked
      logical xc_active, jfit
c     !!! BGJ test !!!

c
c     Make local copy of J and K prefactors in case we need to manipulate them
c
      if (nfock .gt. maxfocks) then
         call errquit('Dimensioning assumption incorrect in fock_2e',0)
      endif
      call dcopy(nfock, fack, 1, kfac, 1)
      call dcopy(nfock, facj, 1, jfac, 1)
c
c     !!! BGJ test !!!
c
c     See if XC is active for this call
c
      if (.not.rtdb_get(bgj_get_rtdb_handle(),'bgj:xc_active',MT_LOG,1,
     &     xc_active)) xc_active = .false.
      if (bgj_print().gt.0)
     &   write(*,*)'*** fock_2e: xc_active ',xc_active
c
c     Set K contribution to be correct even for pure DFT and hybrids
c
      if (xc_active .and. bgj_have_xc()) then
         call dscal(nfock, bgj_kfac(), kfac, 1)
         if (bgj_print().gt.0)
     &      write(*,*)'*** fock_2e: set kfac to ',
     &        (kfac(ifock),ifock=1,nfock)
      endif
c
c     Determine whether J fitting is involved
c
      jfit = bgj_have_j_fit().and.(jfac(1).ne.0d0)
c !!! Edo commented this out before since it caused bugs - ask him about
c !!! it when it is necessary to turn back on for J fit CPKS
      jfit = .false.
c     !!! BGJ test !!!
      if (jfit) then
c
c     get the HF exchange part (must be replaced by DFT X)
c
         vg_fockc = ga_create_atom_blocked
     $        (geom1, ao_basis, 'fock_2e: vg_fockc')
         if (bgj_print().gt.0)
     &      write(*,*)'*** fock_2e: jfit case - original jfac = ',
     &        jfac(1)
         jfac(1)=0d0
         kfac(1)=-0.5d0*bgj_kfac()
         call ao_fock_2e( geom1, ao_basis, 1, jfac, kfac,
     $        tol2e1, oskel1, vg_dens, vg_fock, asym )

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
     $          tol2e1, oskel1, vg_dens, vg_fock,asym )
        endif
      else
         call ao_fock_2e( geom1, ao_basis, nfock, jfac, kfac,
     $                    tol2e1, oskel1, vg_dens, vg_fock, asym )
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
      if (nfock.gt.maxfocks)
     &   call errquit('fock_2e: dimension error in test',0)

      do ifock = 1, nfock
         g_xc(ifock) = ga_create_atom_blocked
     $        (geom1, ao_basis, 'fock_2e:BGJ XC test')
         call ga_zero(g_xc(ifock))
      enddo

      call fock_xc(geom1, nbf,ao_basis,
     ,     nfock, vg_dens, g_xc,scftype)

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
