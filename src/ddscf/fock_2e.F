      subroutine fock_2e( geom1, ao_basis, mfock, jfac, kfac,
     $                    tol2e1, oskel1, vg_dens, vg_fock, asym )
c$Id$
c
c     wrapper routine for ao_fock_2e and riscf_fock_2e
c
c
      implicit none
#include "errquit.fh"
c     !!! BGJ test !!!
#include "global.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c     !!! BGJ test !!!
#include "cscf.fh"
#include "cfock.fh"
#include "bas.fh"
#include "util.fh"
#include "stdio.fh"

      integer geom1, ao_basis        ! [input] parameter handles
      integer nfock,mfock            ! [input] number of Fock matrices
      double precision tol2e1        ! [input] integral selection threshold
      logical oskel1                 ! [input] toggle skeleton Fock matrix
      integer vg_dens(mfock)         ! [input] array of handles to densities
      integer vg_fock(mfock)         ! [input] array of handles to Fock matrices
      logical asym                   ! [input] flag to antisymmetrize Fock matrices
      integer ifock      
      integer maxfocks
      parameter(maxfocks=1000)
      integer g_xc 
      double precision kfac(maxfocks), jfac(maxfocks),ffac
      double precision lfac(maxfocks)
      double precision Exc(2)
      integer nExc
      integer ga_create_atom_blocked
      integer vg_fockc
      integer dims(3), chunk(3),alo(3),ahi(3),blo(3),bhi(3)
      external ga_create_atom_blocked
      logical xc_active, jfit,oprint,odft,cphf_uhf
      oprint= util_print('fock_2e',print_debug)
      cphf_uhf = .false.
      if (.not. rtdb_get(bgj_get_rtdb_handle(), 
     &     'cphf_solve:cphf_uhf', mt_log, 1, cphf_uhf)) then
         cphf_uhf = .false.
      endif
      nfock = mfock
      if(cphf_uhf)nfock = 4
c
c     Make local copy of J and K prefactors in case we need to manipulate them
c
      if (nfock .gt. maxfocks) then
         call errquit('Dimens. assumption incorrect in fock_2e',nfock,
     &       UNKNOWN_ERR)
      endif
c
c     !!! BGJ test !!!
c
c     See if XC is active for this call
c
      if (.not.rtdb_get(bgj_get_rtdb_handle(),'bgj:xc_active',MT_LOG,1,
     &     xc_active)) xc_active = .false.
      if (oprint)
     &   write(luout,*)'*** fock_2e: xc_active ',xc_active
c
c     Set K contribution to be correct even for pure DFT and hybrids
c
      if (xc_active .and. bgj_have_xc()) then
         call dscal(nfock, bgj_kfac(), kfac, 1)
         if (oprint)
     &      write(luout,*)'*** fock_2e: set kfac to ',
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
         if (oprint)
     &      write(luout,*)'*** fock_2e: jfit case - original jfac = ',
     &        jfac(1)
cold         jfac(1)=0d0
cold         kfac(1)=-0.5d0*bgj_kfac()
         call ao_fock_2e( geom1, ao_basis, 1, 0d0, -0.5d0*bgj_kfac(),
     $        tol2e1, oskel1, vg_dens, vg_fock, asym )

         call fock_j_fit(nfock, vg_dens, vg_fockc)
         call ga_add(1.0d0, vg_fockc, 1.0d0, vg_fock,
     &        vg_fock)
         if (.not.ga_destroy(vg_fockc))
     .      call errquit('fock_2e: problem destroying test array',1,
     &       GA_ERR)
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
      if(nfock.ne.1.and.nfock.ne.4) then
         if(ga_nodeid().eq.0) then
            write(0,*) ' WARNING: likely restart: nfock = ',nfock
         endif
         if (.not.rtdb_put(bgj_get_rtdb_handle(),'bgj:xc_active',
     M        MT_LOG,1,.false.)) 
     C        call errquit(' fock2e: bgjrtdbput failed',0,0)
         return
      endif


c      write(*,*)'*** fock_2e: nfock =',nfock
      if (nfock.gt.maxfocks)
     &   call errquit('fock_2e: dimension error in test',0,
     &       UNKNOWN_ERR)
      if (.not. bas_numbf(ao_basis, nbf))
     &   call errquit('fock_2e: basis set error:', 86, BASIS_ERR)
      dims(1) = nfock
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'DFT CPKS RHS',chunk, 
     &     g_xc))
     &   call errquit('fock_2e: could not allocate g_xc',555, GA_ERR)
      call ga_zero(g_xc)
      if(oprint) then
         write(luout,*) 'fock2e: nfock ',nfock
         write(luout,*) 'fock2e: gdens ',vg_dens
         write(luout,*) 'fock2e: g_xc  ',g_xc
         write(luout,*) 'fock2e: g_xc  dims ',dims
      endif
c
c     works with nfock=1 for CS and nfock=4 for OS
c
      if(nfock.ne.1.and.nfock.ne.4) call errquit(
     &     '    fock_2e: unexpected nfock value = ',nfock, UNKNOWN_ERR)
c       
      Exc(1) = 0.0d0
      Exc(2) = 0.0d0
      nExc = 1
      call fock_xc(geom1, nbf,ao_basis,
     &             nfock, vg_dens, g_xc,Exc,nExc,.false.)
      ffac=1d0
      if (xc_active) then
         odft=nfock.eq.4
         ifock=1
         if(odft) ifock=2
         if (oprint)
     .        write(luout,*) 'f2e: Scaling XC matrix',
     ,        ifock,' by',ffac
         alo(1) = ifock
         ahi(1) = ifock
         alo(2) = 1
         ahi(2) = nbf
         alo(3) = 1
         ahi(3) = nbf
         blo(1) = 1
         bhi(1) = nbf
         blo(2) = 1
         bhi(2) = nbf
         blo(3)=-33
         bhi(3)=-33
         call nga_add_patch(ffac,g_xc,alo,ahi,
     .        1d0,vg_fock(ifock),blo,bhi,
     .        vg_fock(ifock),blo,bhi)
         if (oprint) then
            write(luout,*)'*** fock_2e: fock matrix',ifock
            call ga_print(vg_fock(ifock))
         endif
         if(odft) then
            ifock=4
            if (oprint)
     .           write(luout,*) 'f2e: Scaling XC matrix',
     ,           ifock,' by',ffac
            alo(1) = ifock
            ahi(1) = ifock
            call nga_add_patch(ffac,g_xc,alo,ahi,
     .           1d0,vg_fock(ifock),blo,bhi,
     .           vg_fock(ifock),blo,bhi)
            if (oprint) then
               write(luout,*)'*** fock_2e: fock matrix',ifock
               call ga_print(vg_fock(ifock))
            endif
         endif
         
         if (.not.ga_destroy(g_xc)) then
            call errquit('fock_2e: problem destroying test array',1,
     &       GA_ERR)
         endif
      endif
      return
      end
