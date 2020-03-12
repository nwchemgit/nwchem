c
c     Evaluates the TDDFT energy of each root using input data
c     HvD 9/2006, NG 11/2012, DS 1/2013
c
      subroutine tddft_grad_compute_energy(rtdb,ihdl_geom,ihdl_bfao,
     +  tol2e,tda,ipol,nao,nfc,naoc,nocc,nav,nfv,
     +  nroot,o_eig,g_mo,g_xpy,g_xmy,kfac,lhashf,otriplet)
c
c     This subroutine simply evaluates the TDDFT energy using the
c     information read in. The energy expression evaluated is
c
c         0.5 * (<X+Y|A+B|X+Y> + <X-Y|A-B|X-Y>)
c
c     or in more detail
c
c         0.5 * (sum_iasjbt (X+Y)ias (A+B)ias,ibt (X+Y)jbt
c               +sum_iasjbt (X-Y)ias (A-B)ias,ibt (X-Y)jbt)
c
c     Note that (A-B) only has exchange terms and therefore does not
c     couple alpha and beta spins. (A+B) on the other hand contains
c     Coulomb terms and therefore does couple alpha and beta spins.
c
c     This is the energy expression the gradient is based on.
c     The purpose of this routine is simply to test that the data has
c     been read in correctly and the routines invoked are properly 
c     understood.
c
      implicit none
c
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "tddft_grad_util.fh"
#include "rtdb.fh"
c
c     Input:
c
      integer rtdb      ! the runtime database handle
      integer ihdl_geom ! the geometry handle
      integer ihdl_bfao ! the AO basis set handle
      logical tda       ! True if Tamm-Dancoff approximation
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer nao       ! the number of AOs
      integer nfc(2)    ! the number of frozen cores
      integer naoc(2)   ! the number of active occupied orbitals
      integer nocc(2)   ! the number of occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nfv(2)    ! the number of frozen virtuals
      integer nroot     ! the number of roots to be treated
c
c     g_something(1) = alpha, g_something(2) = beta (even in RHF case)
      integer g_mo(2)   ! global arrays to the MOs
      integer g_xpy(2)  ! global arrays to hold X+Y
      integer g_xmy(2)  ! global arrays to hold X-Y
c
      logical lhashf        ! =.true.  hybrid functionals
                            ! =.false. otherwise
      logical otriplet      ! =.true.  triplet excited states
                            ! =.false. singlet excited states
                            ! value does not matter for TDUDFT
c
      double precision kfac       ! the weight of the Hartree-Fock
                                  ! exchange contributions
      double precision o_eig(nao,2) ! orbital eigenvalues
      double precision tol2e      ! 2-electron integral tolerance
c
c     Local for now...
c
      logical oskel
      parameter (oskel=.false.)
c
c     Local:
c
      double precision eroot(nroot) ! the excitation energies
      integer g_hib(2)           ! workspace global arrays
      integer ip                 ! counter for spin
      integer ir                 ! counter for root
      integer g_x                ! work space array for X (AO basis)
      integer g_y                ! work space array for Y (AO basis)
      integer g_apbx             ! work space array for (A+B)X
      integer g_ambx             ! work space array for (A-B)X
      integer g_apby             ! work space array for (A+B)Y
      integer g_amby             ! work space array for (A-B)Y
      integer idim(3)            ! dimensions
      integer ichnk(3)           ! chunking
c     integer ld(3)              ! leading dimensions
      integer alo(3),ahi(3)      ! lower and upper limits for GA patches
      integer iwhich
c     integer i,a                ! orbital counters
c     integer imo,amo            ! orbital labels
c     integer mini,maxi,numi     ! limits
c     integer mina,maxa,numa     ! limits
c     integer ihdl_xpy,iptr_xpy  ! handle and address of X+Y buffer
c     integer ihdl_xmy,iptr_xmy  ! handle and address of X-Y buffer
c
      character*32 pname
c
      pname="tddft_grad_compute_energy: "
      iwhich = 0  ! call to tddft_nga_cont()
c
c     Do the (X+Y) and (X-Y) contributions
c
      do ir = 1,nroot
        eroot(ir) = 0.0d0
      enddo
c
c     Create global arrays for X and Y in AO basis
c
      idim(1) = nroot*ipol
      idim(2) = nao
      idim(3) = nao
      ichnk(1) = nroot*ipol
      ichnk(2) = -1
      ichnk(3) = -1
c
      if (.not.nga_create(mt_dbl,3,idim,'vectors Xuv',ichnk,g_x))
     +    call errquit(pname//'failed to create g_x',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors Yuv',ichnk,g_y))
     +    call errquit(pname//'failed to create g_y',0,GA_ERR)
c
c     Transform (X+Y) to AO basis in g_x and (X-Y) to AO basis in g_y
c
      call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ib",g_mo,g_xpy,"ib",g_x)
      if (.not.tda) then
        call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ib",g_mo,g_xmy,"ib",g_y)
c
c    In order to ensure Y = 0 for later parts of this
c    code, we put g_xpy into g_y such that the lines in 2.3 work the same
c    way they do for RPA.  This is needed so that CIS and RPA (B = 0) are
c    identical.
      else
        call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ib",g_mo,g_xpy,"ib",g_y)
      end if ! tda
c
c     Compute X and Y from (X+Y) and (X-Y) in place
c
      alo(1) = 1
      alo(2) = 1
      alo(3) = 1
      ahi(1) = nroot*ipol
      ahi(2) = nao
      ahi(3) = nao
c
c     With the modification above, the lines here behave
c     like you'd expect (i.e. g_x = X and g_y = Y = 0 for CIS, rather than
c     what was here before which made g_x = 0.50*X and g_y = 0.50*X).
c
      call nga_add_patch(0.5d0,g_x,alo,ahi,0.5d0,g_y,alo,ahi,g_x,
     +                   alo,ahi)
      call nga_add_patch(1.0d0,g_x,alo,ahi,-1.0d0,g_y,alo,ahi,g_y,
     +                   alo,ahi)
c
c     Allocate various workspace arrays
c
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A+B)X',ichnk,g_apbx))
     +    call errquit(pname//'failed to create g_apbx',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A-B)X',ichnk,g_ambx))
     +    call errquit(pname//'failed to create g_ambx',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A+B)Y',ichnk,g_apby))
     +    call errquit(pname//'failed to create g_apby',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'vectors (A-B)Y',ichnk,g_amby))
     +    call errquit(pname//'failed to create g_amby',0,GA_ERR)
c
c     Compute (A+B)X, (A-B)X, (A+B)Y and (A-B)Y
c
      if (otriplet) then
        if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.true.))
     1    call errquit(pname//'failed to set triplet',0,
     2      RTDB_ERR)        
      endif
c
      call tddft_nga_cont(rtdb,ihdl_geom,ihdl_bfao,g_x,g_apbx,g_ambx,
     +     nao,ipol,tol2e,tda,oskel,kfac,lhashf,otriplet,nroot,iwhich)
      call tddft_nga_cont(rtdb,ihdl_geom,ihdl_bfao,g_y,g_apby,g_amby,
     +     nao,ipol,tol2e,tda,oskel,kfac,lhashf,otriplet,nroot,iwhich)
c
c     Dispose of X and Y
c
      if (.not.ga_destroy(g_x)) 
     +  call errquit(pname//'failed to destroy g_x',0,GA_ERR)
      if (.not.ga_destroy(g_y)) 
     +  call errquit(pname//'failed to destroy g_y',0,GA_ERR)
c
c     Compute (A+B)(X+Y) and (A-B)(X-Y)
c
      call nga_add_patch(1.0d0,g_apbx,alo,ahi, 1.0d0,g_apby,alo,ahi,
     +                   g_apbx,alo,ahi)
      call nga_add_patch(1.0d0,g_ambx,alo,ahi,-1.0d0,g_amby,alo,ahi,
     +                   g_ambx,alo,ahi)
c
c     Dispose of (A+B)Y and (A-B)Y
c
      if (.not.ga_destroy(g_apby)) 
     +  call errquit(pname//'failed to destroy g_apby',0,GA_ERR)
      if (.not.ga_destroy(g_amby)) 
     +  call errquit(pname//'failed to destroy g_amby',0,GA_ERR)
c
c     Allocate (A+-B)(X+-Y)ib
c
      do ip = 1, ipol
        idim(1) = nroot
        idim(2) = naoc(ip)
        idim(3) = nav(ip)
        ichnk(1) = nroot
        ichnk(2) = -1
        ichnk(3) = -1
        if (.not.nga_create(mt_dbl,3,idim,'vec hib',ichnk,g_hib(ip)))
     +      call errquit(pname//'failed to create g_hib',0,GA_ERR)
      enddo
c
c     Transform (A+B)(X+Y) to MO basis occupied-virtual block only
c
      call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ib",g_mo,g_apbx,g_hib,"ib")
c
c     Add the orbital energy term (e_a-e_i)(X+Y)ia
c     Drive the parallelization off the distribution of g_xpy.
c
      call tddft_grad_epsilon(g_xpy,g_hib,ipol,nroot,nao,nocc,nfc,o_eig)
c
c     Add sum_ib (X+Y)ib [(A+B)(X+Y)ib] to E
c
      do ip = 1, ipol
        do ir=1,nroot
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = naoc(ip)
          alo(3) = 1
          ahi(3) = nav(ip)
          eroot(ir) = eroot(ir) + nga_ddot_patch(g_xpy(ip),'n',alo,ahi,
     +                                     g_hib(ip),'n',alo,ahi)
        enddo
      enddo
c
      if (.not.tda) then
c
c     Transform (A-B)(X-Y) to MO basis occupied-virtual block only
c
      call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,nfv,
     +     nroot,1.0d0,0.0d0,"ib",g_mo,g_ambx,g_hib,"ib")
c
c     Add the orbital energy term (e_a-e_i)(X-Y)ia
c     Drive the parallelization off the distribution of g_xmy.
c
      call tddft_grad_epsilon(g_xmy,g_hib,ipol,nroot,nao,nocc,nfc,o_eig)
c
c     Add sum_ib (X-Y)ib [(A-B)(X-Y)ib] to E
c
      do ip = 1, ipol
        do ir=1,nroot
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = naoc(ip)
          alo(3) = 1
          ahi(3) = nav(ip)
          eroot(ir) = eroot(ir) + nga_ddot_patch(g_xmy(ip),'n',alo,ahi,
     +                                     g_hib(ip),'n',alo,ahi)
        enddo
      enddo
c
      endif  !tda
c
c     Dispose of g_hib
c
      do ip = 1, ipol
        if (.not.ga_destroy(g_hib(ip))) 
     +   call errquit(pname//'failed to destroy g_hib',0, GA_ERR)
      enddo
c
c     Dispose of (A+B)(X+Y) and (A+B)(X-Y)
c
      if (.not.ga_destroy(g_apbx)) 
     +  call errquit(pname//'failed to destroy g_apbx',0,GA_ERR)
       if (.not.ga_destroy(g_ambx)) 
     +  call errquit(pname//'failed to destroy g_ambx',0,GA_ERR)
c
c     Scale the energy by half for full TDDFT
c
      do ir = 1, nroot
       if (.not.tda) eroot(ir) = 0.5d0*eroot(ir)
       if (ga_nodeid().eq.0) then
         write(luout,"('TDDFT Energy Check(',i3,') = ',f22.14)")
     &     ir,eroot(ir)
       endif
      enddo
c
      end
c $Id$
