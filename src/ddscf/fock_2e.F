      subroutine fock_2e( geom1, ao_basis, nfock, facj, fack,
     $                    tol2e1, oskel1, vg_dens, vg_fock, asym )
c$Id: fock_2e.F,v 1.22 2001-11-29 22:33:52 bert Exp $
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
#include "bas.fh"

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
      parameter(maxfocks=1000)
      integer g_xc ! Assume no more than 1000 fock matrices !!!
      double precision kfac(maxfocks), jfac(maxfocks),ffac
      double precision lfac(maxfocks)
      integer ga_create_atom_blocked
      integer vg_fockc
      integer dims(3), chunk(3)
      integer k_buf,l_buf
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
      else if (nriscf.eq.10) then
        do ifock=1,nfock
            lfac(ifock)=0d0
        enddo
c
c       get HF exchange part in normal way
c
        call ao_fock_2e( geom1, ao_basis, nfock, lfac, kfac,
     $        tol2e1, oskel1, vg_dens, vg_fock, asym )
c
c       get HF coulomb part from RISCF integrals
c
        call riscf_fock_2e( geom1, ao_basis, riscf_basis, nfock,
     $       jfac, lfac, tol2e1, oskel1, vg_dens, vg_fock )
c
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
      if (.not. bas_numbf(ao_basis, nbf))
     &   call errquit('fock_2e: basis set error:', 86)
      write(0,*) ' nfock ',nfock
      dims(1) = nfock
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'DFT CPKS RHS',chunk, 
     &     g_xc))
     &   call errquit('fock_2e: could not allocate g_xc',555)
      call ga_zero(g_xc)
cedo      do ifock = 1, nfock
cedo         g_xc(ifock) = ga_create_atom_blocked
cedo     $        (geom1, ao_basis, 'fock_2e:BGJ XC test')
cedo         call ga_zero(g_xc(ifock))
cedo      enddo

      call fock_xc(geom1, nbf,ao_basis,
     ,     nfock, vg_dens, g_xc,scftype)
chack
      if(ga_nodeid().eq.0) then
         if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'charge',l_buf,k_buf))
     &        call errquit('fock_2e: hack allocate failed',0)
      endif
chack
      do ifock = 1, nfock
!         if (xc_active.and.facj(ifock).ne.0d0) then
         if (xc_active) then
c     Scale the XC matrices by the same scaling factor as J in this case
c     since there is no explicit way to pass in an XC scaling factor, yet
c     the 2nd derivative code uses a scale factor of 8 for the OO part of
c     the CPSCF RHS that needs to be applied to XC also.
            if (facj(ifock) .ne. 1d0 .and. facj(ifock) .ne. 8d0) then
               write(*,*)'cant Scale XC matrix',ifock,' by',facj(ifock)
               ffac=fack(ifock)
               write(*,*)'*** Scaling XC matrix',ifock,' by',ffac
!               call errquit('Unexpected XC scale factor in fock_2e',0)
            else
               ffac=facj(ifock)
            endif
!            call ga_add(facj(ifock), g_xc(ifock), 1.0d0, vg_fock(ifock),
!     &           vg_fock(ifock))
            if(ga_nodeid().eq.0) then
               call dft_3dget(g_xc,dbl_mb(k_buf),ifock,nbf)
               call dscal(nbf, ffac, dbl_mb(k_buf), 1)
              call ga_acc(vg_fock(ifock), 1, nbf, 1, nbf, 
     ,              dbl_mb(k_buf), nbf,1d0)
            endif
         endif
         if (bgj_print() .gt. 1) then
            write(*,*)'*** fock_2e: fock matrix',ifock
            call ga_print(vg_fock(ifock))
         endif
      enddo
chack
      if(ga_nodeid().eq.0) then
         if (.not.MA_pop_stack(l_buf))
     &        call errquit('fock_2e: hack popstack failed',0)
      endif
chack
      if (.not.ga_destroy(g_xc)) then
         call errquit('fock_2e: problem destroying test array',1)
      endif

c     !!! BGJ test !!!

      return
      end
      subroutine dft_3dget(g_fock,buf,imat,nbf)
C$Id: fock_2e.F,v 1.22 2001-11-29 22:33:52 bert Exp $
      implicit none
      integer g_fock
      integer imat
      integer nbf
      double precision buf(*)
c      
      integer lo(3)  !array of starting indices for ga   [input]
      integer hi(3)  !  array of ending indices for ga     [input]
      integer ld(2) !array specifying leading dimensions for buffer 
      ld(1)=1  
      ld(2)=nbf
      lo(1)=imat
      hi(1)=imat
      lo(2)=1
      hi(2)=nbf
      lo(3)=1
      hi(3)=nbf
      call nga_get(g_fock, lo, hi, buf, ld)
      return
      end
