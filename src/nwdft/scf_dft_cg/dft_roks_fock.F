      subroutine dft_roks_fock(rtdb, geom, basis, nclosed, nopen, nmo,
     $     tol2e, g_vecs,
     $     eone, etwo, edisp, Exc, g_fcv, g_fpv,
     $     g_fcp, nExc, oskel, iter)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "util.fh"
#include "cscfps.fh"
#include "cosmo.fh"
#include "geom.fh"
#include "case.fh"
#include "frozemb.fh"
c     
c     Build fock matrices and energies for ROKS wavefunctions
c     -------------------------------------------------------
c
c     This routine builds blocks of the Fock matrix for ROHF and ROKS
c     problems. Originally this routine dealt with ROHF only, the ROKS
c     was added later and the corresponding expressions derived by 
c     analogy from the ROHF expressions.
c
c     In the derivation of the ROHF expressions the fact that the Fock
c     matrix is linear in the density matrices is exploited. As a result
c     a number of equivalent ways to express the Fock matrix can be
c     used. In DFT this is not possible as the Kohn-Sham matrix is 
c     not a linear function of the density. So for the Kohn-Sham matrix
c     the expressions have to be given in terms of the alpha- and beta-
c     density Da, Db and the alpha- and beta- Fock(Kohn-Sham) matrices
c     Fa and Fb. 
c
c     The original Fock matrix expressions are given as:
c
c       F_CO = F(C) + Fc(O)                                     (g_fcp)
c       F_CV = F(C) + Fc(O) + Fe(O)                             (g_fcv)
c       F_OV = F(C) + Fc(O) + 2Fe(O)                            (g_fpv)
c
c     In these expressions F is the total Fock matrix and Fc is the
c     coulomb and Fe the exchange terms. On the left-hand-side the 
c     various matrix blocks are given:
c
c       CO -- the closed shell orbital - open shell orbital block
c       CV -- the closed shell orbital - virtual orbital block
c       OV -- the open shell orbital   - virtual orbital block
c
c     On the right-hand-side the terms contributing to these blocks are
c     given. Each term is a function of the density matrix given in 
c     round brackets. 
c
c     To arrive at the ROKS expressions we express ROHF ones in terms
c     of Fa, Fb, Da, and Db. Assuming that in the high spin open shell
c     case there will be more alpha than beta electrons we can express
c     the density matrixes as follows:
c
c       C = Dc = 2Db   (the closed shell density matrix)
c       O = Do = Da-Db (the open shell density matrix)
c
c     Using the fact the ROHF Fock matrices are linear in the density
c     matrices the terms in the Fock matrix expressions can be recast
c     as:
c
c       ( F(2Db) + Fc(Da-Db)
c       ( F(2Db) + Fc(Da-Db) + Fe(Da-Db)
c       ( F(2Db) + Fc(Da-Db) + 2Fe(Da-Db)
c
c       ( Fc(2Db) + Fe(2Db) + Fc(Da-Db)
c       ( Fc(2Db) + Fe(2Db) + Fc(Da-Db) + Fe(Da-Db)
c       ( Fc(2Db) + Fe(2Db) + Fc(Da-Db) + 2Fe(Da-Db)
c
c       ( Fc(Da) + Fc(Db) + 2Fe(Db)          = F(CO)
c       ( Fc(Da) + Fc(Db) + Fe(Db) + Fe(Da)  = F(CV)
c       ( Fc(Da) + Fc(Db) + 2Fe(Da)          = F(OV)
c
c     In DFT the Fe maps onto the exchange-correlation contributions. 
c     Obviously as there are no exchange cross terms between Da and Db
c     this means we can equate Fe(Da) with the alpha spin Kohn-Sham
c     matrix and Fe(Db) with the beta spin Kohn-Sham. For the 
c     correlation contributions there still is some ambiguity about this
c     but we'll assume that we can do the same in that case.
c
c     Arguments
c     
      integer rtdb
      integer geom, basis       ! [input] Handles to environment
      integer nclosed, nopen    ! [input] no. of closed and open shells
      integer nmo               ! [input] no. of MO vectors
      integer g_vecs            ! [input] MO vectors
      double precision tol2e    ! [input] Integral tolerance
      double precision eone, etwo ! [output] Energy terms
      integer g_fcv             ! [output] Closed-virtual Fock matrix
      integer g_fpv             ! [output] Open-virtual Fock matrix
      integer g_fcp             ! [output] Closed-open Fock matrix
      logical oskel             ! [input] use skeleton symmetry
      integer iter              ! [input] iteration number
      integer nExc              ! [input] no. exch-corr components
      double precision Exc(2)   ! [output] exchange-correlation energy
      double precision Edisp    ! [input] dispersion correction
c     
c     Local variables
c     
      double precision eone_closed, eone_open, etwo_closed, etwo_open
      double precision errmax, dum, Ebq
      double precision time1_2e, time2_2e, time1_xc, time2_xc
      integer nset
      parameter(nset=7)
      integer g_cdens, g_odens, g_cfock, g_hcore, g_ocoul, g_oexch
      integer g_dens
      integer g_c_xc ! closed shell KS contributions (closed shell only)
      integer g_a_xc ! alpha spin KS contributions (open shell only)
      integer g_b_xc ! beta spin KS contributions (open shell only)
      integer g_adens ! alpha density matrix for DFT (open shell only)
      integer g_bdens ! beta density matrix for DFT (open shell only)
      integer nbf, nvir, nfock, nfockxc
      double precision two, one, zero, mone, half
      parameter (two=2.0d0, one=1.0d0, zero=0.0d0, mone=-1.0d0)
      parameter (half=0.5d0)
      double precision sym_tol  ! Threshold for detecting symmetry
                                ! breaking
      double precision jfactor(nset),kfactor(nset)
      integer iv_dens(nset),iv_fock(nset)
      integer g_tmp(nset)
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical  xc_chkdispauto
      external xc_chkdispauto
c
      integer g_cos1e
      dimension g_cos1e(3)
      double precision ecosmo
      logical osome
      logical disp
      logical xc_gotxc
      external xc_gotxc
      double precision xc_hfexch
      external xc_hfexch
c
      integer itype_wfn
c     AOR begin                                                                 
      double precision exdm
      integer ixdm_v, ixdm_a, ixdm_ml, natoms
      external xc_xdm_lxdm
      integer xc_xdm_lxdm
c     AOR end         
      logical odebug
      character*255 cosmo_file
c
c     Initialise
c     
      odebug = util_print('dft_rohf_fock',print_debug)
c
      if (.not. bas_numbf(basis, nbf))
     $   call errquit('dft_roks_fock: bas_numbf failed',0, BASIS_ERR)
      call rohf_jkfac(jfactor,kfactor)
      nvir = nmo - nclosed - nopen
c     
c     hcore = 1-e integrals
c     cdens = closed shell AO density
c     cfock = closed shell AO fock
c     odens = open shell density
c     ocoul = open shell coulomb
c     oexch = open shell exchange
c
      g_a_xc  = 0
      g_b_xc  = 0
      g_adens = 0
      g_bdens = 0
c
      g_cdens = ga_create_atom_blocked(geom, basis, 'rohf: cdens')
      g_cfock = ga_create_atom_blocked(geom, basis, 'rohf: cfock')
      g_hcore = ga_create_atom_blocked(geom, basis, 'rohf: hcore')
      if (xc_gotxc()) then
        if (nopen .eq. 0 ) then
          g_c_xc = ga_create_atom_blocked(geom, basis, 'roks: c xc')
        else
          g_a_xc = ga_create_atom_blocked(geom, basis, 'roks: a xc')
          g_b_xc = ga_create_atom_blocked(geom, basis, 'roks: b xc')
          g_adens = ga_create_atom_blocked(geom, basis, 'roks: adens')
          g_bdens = ga_create_atom_blocked(geom, basis, 'roks: bdens')
        endif
      endif
      if (nopen .gt. 0) then
         g_odens = ga_create_atom_blocked(geom, basis, 'rohf: odens')
         g_ocoul = ga_create_atom_blocked(geom, basis, 'rohf: oexch')
         g_oexch = ga_create_atom_blocked(geom, basis, 'rohf: ocoul')
      endif
c     
c     Compute closed and open shell density matrices
c     
      call rohf_dens(g_vecs, nbf, nclosed, nopen, g_cdens, g_odens)
      if (odebug) then
         call ga_print(g_vecs)
         call ga_print(g_cdens)
         if (nopen .gt. 0) call ga_print(g_odens)
      endif
c
c     Symmetrize the density matrix.
c     If the HOMO and LUMO orbitals are of the same irrep and degenerate
c     then they can rotate in arbitrary ways among eachother. As a
c     result it cannot be guaranteed that a symmetric density matrix
c     can be calculated even if the orbitals are symmetry adapted.
c     Hence density matrix symmetrization is essential. This phenomenon
c     was observed when running the c2h4 test case using D2d symmetry.
c
      if (oskel) then
         call ga_copy(g_cdens, g_hcore)
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .true., g_cdens)
         if (oscfps) call pstat_off(ps_sym_sym)
         call ga_dadd(one, g_hcore, mone, g_cdens, g_hcore)
         call ga_maxelt(g_hcore, errmax)
         sym_tol = tol2e*10.0d0
         if (errmax .gt. sym_tol) then
           if (ga_nodeid().eq.0) then
             write(luout,*)'Warning: symmetry breaking in RKS: ',
     &                     errmax
             call util_flush(luout)
           endif
         endif
         if (nopen .gt. 0) then
            call ga_copy(g_odens, g_hcore)
            if (oscfps) call pstat_on(ps_sym_sym)
            call sym_symmetrize(geom, basis, .true., g_odens)
            if (oscfps) call pstat_off(ps_sym_sym)
            call ga_dadd(one, g_hcore, mone, g_odens, g_hcore)
            call ga_maxelt(g_hcore, errmax)
            if (errmax .gt. sym_tol) then
              if (ga_nodeid().eq.0) then
                write(luout,*)'Warning: symmetry breaking in ROKS: ',
     &                        errmax
                call util_flush(luout)
              endif
            endif
         endif
      endif
      if (xc_gotxc().and.nopen.gt.0) then
        call ga_copy(g_cdens,g_bdens)
        call ga_dscal(g_bdens,0.5d0)
        call ga_add(half,g_cdens,one,g_odens,g_adens)
      endif
c
c     One-electron integrals and energies
c     
      call ga_zero(g_hcore)
      if (rtdb_get(rtdb,'dft:bq_energy',mt_dbl,1,Ebq)) then
        call int_1e_ga(basis, basis, g_hcore, 'bq_pot', oskel)
        Ebq = ga_ddot(g_cdens,g_hcore)
        if (nopen.gt.0) then
          Ebq = Ebq + ga_ddot(g_odens,g_hcore)
        endif
        if (.not. rtdb_put(rtdb,'dft:bq_energy',mt_dbl,1,Ebq))
     $     call errquit('dft_roks_fock: writing bq_energy failed',
     $     0, RTDB_ERR)
      endif
      call int_1e_ga(basis, basis, g_hcore, 'kinetic', oskel)
      call int_1e_ga(basis, basis, g_hcore, 'potential', oskel)
      if (odebug)  call ga_print(g_hcore)
c
c     cosmo charges, potential and energy contribution
c
      osome=.false.
      osome=osome.and.ga_nodeid().eq.0
      if(cosmo_on.and.cosmo_phase.eq.2) then
       cosmo_file = "cosmo.xyz"
       g_dens = ga_create_atom_blocked(geom,basis,'density matrix')
       call ga_zero(g_dens)
       call ga_copy(g_cdens,g_dens)
       if (nopen.gt.0) call ga_add(1.d0,g_cdens,one,g_odens,g_dens)
        call cosmo_charges_from_dmat(rtdb, basis, geom, ecosmo, odebug,
     &                                  1,           ! 1 = closed  shell
     &                                  g_dens,      ! input density
     &                                  cosmo_file)  ! cosmo charges file name
       call int_1e_ga(basis,basis,g_hcore,'cos_chg_pot',.false.)
       if (.not.ga_destroy(g_dens))
     &  call errquit('dft_roks_fock: ga_destroy failed g_dens',0,GA_ERR)
      endif  ! cosmo check
c
c     add in frozen embedding
c
      if (frozemb) then
        call ga_dadd(1.d0,g_hcore,1.d0,g_frozemb,g_hcore)
      end if
c
      call ga_sync()
c
      eone_closed = ga_ddot(g_hcore, g_cdens)
      eone_open = 0.0d0
      if (nopen .gt. 0) eone_open = ga_ddot(g_hcore, g_odens)
      eone = eone_closed + eone_open
c     
c     Generate two-electron closed- and open-shell Fock components
c     
c     Fock              1.0*F[ D ]      (g_cfock)
c     
c     Open coulomb      1.0*J[ D ]      (g_ocoul) 
c     
c     Open exchange    -0.5*K[ D ]      (g_oexch) 
c     
      iv_dens(1) =  g_cdens
      iv_fock(1) =  g_cfock
      call ga_zero(g_cfock)
      if (xc_gotxc()) then
        if (nopen.eq.0) then
          iv_dens(4) =  g_cdens
          iv_fock(4) =  g_c_xc
          call ga_zero(g_c_xc)
          nfockxc = 1
        else
          iv_fock(4) =  g_a_xc
          iv_fock(5) =  g_b_xc
          call ga_zero(g_a_xc)
          call ga_zero(g_b_xc)
          iv_dens(4) =  g_adens
          iv_dens(5) =  g_adens
          iv_dens(6) =  g_bdens
          iv_dens(7) =  g_bdens
          nfockxc = 8
        endif
      endif
      nfock = 1
      if (nopen .gt. 0) then
         iv_dens(2) =  g_odens
         iv_dens(3) =  g_odens
         iv_fock(2) =  g_ocoul
         iv_fock(3) =  g_oexch
         call ga_zero(g_oexch)
         call ga_zero(g_ocoul)
         nfock = 3
      endif
c
      if(odebug) time1_2e=util_cpusec()   ! start 2e build time
c
C     don't do it if RI HESSIAN is active
      call do_riscf (.false.)
C     Dyalls method for RHF
      if ( nopen .eq. 0 ) call riscf_dovecs(.true., g_vecs)
c
      if (.not.cam_exch) then
        jfactor(1)=1d0
        jfactor(2)=1d0
        jfactor(3)=0d0
        kfactor(1)=-0.5d0*xc_hfexch()
        kfactor(2)=0d0
        kfactor(3)=-0.5d0*xc_hfexch()
        call fock_2e(geom, basis, nfock, jfactor, kfactor,
     $     tol2e, oskel, iv_dens, iv_fock, .false.)
      else 
c
c       For attenuated calculations
c
c       get some work space
c
        g_tmp(1) = ga_create_atom_blocked(geom, basis,
     $                                    'dft_roks_fock: tmp1')
        call ga_zero(g_tmp(1))
        if (nopen.gt.0) then
          g_tmp(2) = ga_create_atom_blocked(geom, basis,
     $                                      'dft_roks_fock: tmp2')
          g_tmp(3) = ga_create_atom_blocked(geom, basis,
     $                                      'dft_roks_fock: tmp3')
          call ga_zero(g_tmp(2))
          call ga_zero(g_tmp(3))
        endif
c
c       calculate the CAM exchange
c
        call case_setflags(.true.)
        jfactor(1)=0d0
        jfactor(2)=0d0
        jfactor(3)=0d0
        kfactor(1)=-0.5d0*xc_hfexch()
        kfactor(2)=0d0
        kfactor(3)=-0.5d0*xc_hfexch()
        call fock_2e(geom, basis, nfock, jfactor, kfactor,
     $     tol2e, oskel, iv_dens, g_tmp, .false., .false.)
        call ga_dadd(1d0,iv_fock,1d0,g_tmp(1),iv_fock)
        if (nopen.gt.0) then
          call ga_dadd(1d0,iv_fock(2),1d0,g_tmp(2),iv_fock(2))
          call ga_dadd(1d0,iv_fock(3),1d0,g_tmp(3),iv_fock(3))
        endif
c
c       calculate the full Coulomb
c
        call ga_zero(g_tmp(1))
        if (nopen.gt.0) then
          call ga_zero(g_tmp(2))
          call ga_zero(g_tmp(3))
        endif
c
        call case_setflags(.false.)
        jfactor(1)=1d0
        jfactor(2)=1d0
        jfactor(3)=0d0
        kfactor(1)=0d0
        kfactor(2)=0d0
        kfactor(3)=0d0
        call fock_2e(geom, basis, nfock, jfactor, kfactor,
     $     tol2e, oskel, iv_dens, g_tmp, .false., .true.)
        call ga_dadd(1d0,iv_fock,1d0,g_tmp(1),iv_fock)
        if (nopen.gt.0) then
          call ga_dadd(1d0,iv_fock(2),1d0,g_tmp(2),iv_fock(2))
          call ga_dadd(1d0,iv_fock(3),1d0,g_tmp(3),iv_fock(3))
        endif
c
c       destroy work space
        if (.not. ga_destroy(g_tmp(1))) call errquit
     $             ('xc_getv: ga corrupt?',0, GA_ERR)
        if (nopen.gt.0) then
          if (.not. ga_destroy(g_tmp(2))) call errquit
     $               ('xc_getv: ga corrupt?',0, GA_ERR)
          if (.not. ga_destroy(g_tmp(3))) call errquit
     $               ('xc_getv: ga corrupt?',0, GA_ERR)
        endif
c
      end if
c
      call do_riscf (.true.)
      if ( nopen .eq. 0 ) call riscf_dovecs(.false., g_vecs)
      if(util_print('dft timings', print_high))
     &      time2_2e=util_cpusec()   ! end 2e build time
c
c     print fock_2e build time
c
      if(util_print('dft timings', print_high)) then
       if (ga_nodeid().eq.0) then
         write(*,"(4x,'Fock_2e Build Time:',F13.1,'s')")
     &              time2_2e-time1_2e
       endif
      end if
c
c     do DFT stuff
c
      if(util_print('dft timings', print_high))
     &      time1_xc=util_cpusec()   ! start xc build time
      Exc(1) = 0.0d0
      Exc(2) = 0.0d0
      if (xc_gotxc()) then
        call fock_xc(geom, nbf, basis, nfockxc, iv_dens(4), iv_fock(4),
     +               Exc, nExc,.false.)
c     AOR begin
        call xc_xdm_init(rtdb,ixdm_v,ixdm_ml)
        if (.not. rtdb_put(rtdb,'dft:xdmsave', mt_log, 1, .false.)) 
     &  call errquit('dftrocks: rtdbput failed ',0,RTDB_ERR)

        if (xc_xdm_lxdm() .eq. 1) then
           if(.not.geom_ncent(geom, natoms)) call 
     c          errquit('geomncent failed ',0,GEOM_ERR)
           call xc_xdm(rtdb,iv_dens(4),iv_fock(4),natoms,nexc,exdm,
     D          dum,dbl_mb(ixdm_v),dbl_mb(ixdm_ml),'energy')
         exc = exc + exdm
c     deallocate xdm arrays
         call xc_xdm_cleanup(rtdb)
      endif
c     AOR end
      endif
      if(util_print('dft timings', print_high))
     &      time2_xc=util_cpusec()   ! end xc build time
c
c     print fock_xc build time
c
      if(util_print('dft timings', print_high)) then
        if (ga_nodeid().eq.0) then
          write(*,"(4x,'Fock_xc Build Time:',F13.1,'s')")
     &              time2_xc-time1_xc
        endif
      end if
c
c     enddo DFT stuff
c     
c     Energy evaluation
c     
      etwo_closed = 0.5d0*ga_ddot(g_cdens,g_cfock)
      etwo_open = 0.0d0
      if (nopen .gt. 0) etwo_open = 
     $     ga_ddot(g_cdens,g_ocoul) + 
     $     ga_ddot(g_cdens,g_oexch) +
     $     0.5d0*ga_ddot(g_odens,g_ocoul) +
     $     ga_ddot(g_odens,g_oexch)
      etwo = etwo_closed + etwo_open + Edisp
c     
      call ga_dadd(One, g_hcore, One, g_cfock, g_cfock)
      if (xc_gotxc().and.nopen.eq.0) then
        call ga_dadd(One, g_c_xc, One, g_cfock, g_cfock)
      endif
c     
c     Assemble Fock components
c     F(closed-open)    = F(closed) + Fcoul(open)                  (fcp)
c     F(closed-virtual) = F(closed) + Fcoul(open) + Fexch(open)    (fcv)
c     F(open-virtual)   = F(closed) + Fcoul(open) + 2*Fexch(open)  (fpv)
c     
      if (nopen .gt. 0) then
c     
c        Strive for g_cfock=FCV(AO), g_ocoul=FCP(AO), g_oexch=FPV(AO)
c     
         call ga_dadd(one, g_ocoul, one, g_cfock, g_ocoul)
         call ga_dadd(one, g_ocoul, one, g_oexch, g_cfock)
         call ga_dadd(one, g_cfock, one, g_oexch, g_oexch)
         if (xc_gotxc()) then
c
c          The factors below were originally derived as twice the size.
c          The adjustment was introduced by analogy with the factors
c          1/2 that appear in the kfac arrays. For now this seems to
c          work but the proper justification for this is not established
c          at present.
c
           call ga_dadd(One, g_b_xc, One, g_ocoul, g_ocoul)
           call ga_dadd(Half, g_a_xc, One, g_cfock, g_cfock)
           call ga_dadd(Half, g_b_xc, One, g_cfock, g_cfock)
           call ga_dadd(One, g_a_xc, One, g_oexch, g_oexch)
         endif
      endif
c     
c     Symmetrize the AO components
c     
      if (oskel) then
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .false., g_cfock)
         if (oscfps) call pstat_off(ps_sym_sym)
         if (nopen .gt. 0) then
            if (oscfps) call pstat_on(ps_sym_sym)
            call sym_symmetrize(geom, basis, .false., g_ocoul)
            if (oscfps) call pstat_off(ps_sym_sym)
            if (oscfps) call pstat_on(ps_sym_sym)
            call sym_symmetrize(geom, basis, .false., g_oexch)
            if (oscfps) call pstat_off(ps_sym_sym)
         endif
      endif
c
      if (util_print('dft_roks_fock',print_debug)) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' Symmetrized combined AO fock matrices '
            call util_flush(6)
         endif
         call ga_print(g_cfock)
         if (nopen .gt. 0) then
            call ga_print(g_ocoul)
            call ga_print(g_oexch)
         endif
      endif
c     
c     Transform Fock matrices to MO representation using hcore as
c     scratch space
c     
      call two_index_transf(g_cfock, g_vecs, g_vecs, g_hcore, g_fcv)
      
      if (nopen .gt. 0) then
         call two_index_transf(g_ocoul, g_vecs, g_vecs, g_hcore, g_fcp)
         call two_index_transf(g_oexch, g_vecs, g_vecs, g_hcore, g_fpv)
      endif
c     
      if (util_print('dft_roks_fock',print_debug)) then
         call ga_print(g_fcv)
         if (nopen .gt. 0) then
            call ga_print(g_fcp)
            call ga_print(g_fpv)
         endif
         if (ga_nodeid().eq.0) then
            write(6,77) eone_open, eone_closed, etwo_open, etwo_closed
 77         format(1p,'  e1o=',d16.8,', e1c=',d16.8,', e2o=',d16.8,
     $           ', e2c=',d16.8)
            call util_flush(6)
         endif
      endif
c     
      if (.not. ga_destroy(g_hcore)) call errquit
     $     ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_cdens)) call errquit
     $     ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_cfock)) call errquit
     $     ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
      if (xc_gotxc()) then
        if (nopen.eq.0) then
          if (.not. ga_destroy(g_c_xc)) call errquit
     $       ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
        else
          if (.not. ga_destroy(g_a_xc)) call errquit
     $       ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
          if (.not. ga_destroy(g_b_xc)) call errquit
     $       ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
          if (.not. ga_destroy(g_adens)) call errquit
     $       ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
          if (.not. ga_destroy(g_bdens)) call errquit
     $       ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
        endif
      endif
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_ocoul)) call errquit
     $        ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
         if (.not. ga_destroy(g_oexch)) call errquit
     $        ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
         if (.not. ga_destroy(g_odens)) call errquit
     $        ('dft_roks_fock: ga_destroy failed',0, GA_ERR)
      endif
c     
      end
