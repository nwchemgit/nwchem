      subroutine fock_2e( geom1, ao_basis, mfock, facj, fack,
     $                    tol2e1, oskel1, vg_dens, vg_fock, asym )
c$Id: fock_2e.F,v 1.18 2001-06-21 18:49:37 zzhang Exp $
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
      integer nfock, mfock, sfock                  ! [input] number of Fock matrices
      double precision facj(mfock)   ! [input] Coulomb prefactor
      double precision fack(mfock)   ! [input] exchange prefactor
      double precision tol2e1        ! [input] integral selection threshold
      logical oskel1                 ! [input] toggle skeleton Fock matrix
      integer vg_dens(mfock)         ! [input] array of handles to densities
      integer vg_fock(mfock)         ! [input] array of handles to Fock matrices
      logical asym                   ! [input] flag to antisymmetrize Fock matrices
c     !!! BGJ test !!!
      integer ifock, ibase, ivec, nvec
      integer maxfocks
      parameter(maxfocks=100)
      integer g_xc(maxfocks) ! Assume no more than 100 fock matrices !!!
      integer g_xctmp(maxfocks)
      double precision kfac(maxfocks), jfac(maxfocks)
      integer ga_create_atom_blocked
      integer vg_fockc, g_tmp
      external ga_create_atom_blocked
      logical xc_active, jfit
      logical cphf_uhf, cphf_tddft
c     !!! BGJ test !!!
      cphf_uhf = .false.
      cphf_tddft = .false.
      if (.not. rtdb_get(bgj_get_rtdb_handle(), 
     &     'cphf_solve:cphf_uhf', mt_log, 1, cphf_uhf)) then
         cphf_uhf = .false.
      endif
      if (.not. rtdb_get(bgj_get_rtdb_handle(), 
     &     'cphf_solve:cphf_tddft', mt_log, 1, cphf_tddft)) then
         cphf_tddft = .false.
      endif
      nfock = mfock
      if(cphf_uhf)nfock = 4
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
      if (cphf_uhf)then 
c
c     cphf_uhf .eq. true: called from uhf_energy/cphf_solve.
c     vg_dens(2) should have the beta density and g_xc(2)
c     returns the beta xc matrix 
c
         g_tmp = vg_dens(2)
         vg_dens(2) = vg_dens(3)
      endif
      if(cphf_tddft)sfock=nfock 
      call fock_xc(geom1, nbf,ao_basis,
     ,     nfock, vg_dens, g_xc,scftype)
      if(cphf_tddft)then 
         nfock = sfock 
         write(*,*)"please verify this: output vg_focks"
         do ifock = 1,nfock
            call ga_print(vg_fock(ifock))
         enddo
c     
c     take the sum of vg_fock(1) and vg_fock(3). 
c     
         if(scftype.eq.'UHF')then
            do ifock = 1, nfock/4
               ibase = 4*(ifock-1)
               call ga_add(1.0d0, vg_fock(ibase+1), 1.0d0, 
     &              vg_fock(ibase+3), vg_fock(ibase+1))
               call ga_copy(vg_fock(ibase+1), vg_fock(ibase+3))
            enddo
            write(*,*)"xc matrices in the order of g_xc(1:nvec, 1:ipol)" 
         endif
         do ifock = 1, nfock
            call ga_print(g_xc(ifock))
         enddo
c     
c     The xc matrices are in the order of g_xc(1:nvec, 1:ipol). 
c     Map these to the same order as vg_focks(1:4, nvec)
c     
         if(scftype.eq.'UHF')then
            nvec = nfock/4
            do ivec = 1,nvec
               g_xctmp(4*(ivec-1)+1)=g_xc(ivec) ! alpha xc
               g_xctmp(4*(ivec-1)+2)=g_xc(2*nvec+ivec) ! zero xc 
               g_xctmp(4*(ivec-1)+3)=g_xc(nvec+ivec) !beta xc
               g_xctmp(4*(ivec-1)+4)=g_xc(3*nvec+ivec) ! zero xc 
            enddo
            do ifock = 1, nfock
               g_xc(ifock)=g_xctmp(ifock)
            enddo
            do ifock = 1, nfock
               call ga_print(vg_fock(ifock))
               call ga_print(g_xc(ifock))
            enddo
         else if(scftype.eq.'RHF')then
            do ifock = 1,nfock
               call ga_dscal(vg_fock(ifock),4.0d0)
               call ga_dscal(g_xc(ifock),4.0d0)
            enddo
         else
            call errquit("fock_2e: tddft RHF/UHF", 0)
         endif
      endif
      if (cphf_uhf) then 
         vg_dens(2) = g_tmp     ! restore vg_dens
         call ga_copy(g_xc(1), vg_fock(5))
         call ga_copy(g_xc(2), vg_fock(6))
      endif
      do ifock = 1, nfock
         if (xc_active) then
c     Scale the XC matrices by the same scaling factor as J in this case
c     since there is no explicit way to pass in an XC scaling factor, yet
c     the 2nd derivative code uses a scale factor of 8 for the OO part of
c     the CPSCF RHS that needs to be applied to XC also.
            if (facj(ifock) .ne. 1d0 .and. facj(ifock) .ne. 8d0 .and. 
     &           (.not.cphf_uhf) .and. (.not.cphf_tddft)) then
               write(*,*)'*** Scaling XC matrix',ifock,' by',facj(ifock)
               call errquit
     &              ('Unexpected XC scale factor in fock_2e',0)
            endif
            if (.not.cphf_uhf)then 
               call ga_add(facj(ifock), g_xc(ifock), 1.0d0, 
     &              vg_fock(ifock), vg_fock(ifock))
               call ga_print(vg_fock(ifock))
            endif
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
