c
c     Computes the excited-state gradient for every root  
c     HvD 11/2007, NG 10/2012, DWS 4/2013
c
      subroutine tddft_grad_compute_g(rtdb,ihdl_geom,ihdl_bfao,tda,nat,
     +   nao,nfc,naoc,nocc,nav,nfv,ipol,nroots,oskel,tol2e,g_mo,
     +   g_p,g_wp,g_xpy,g_xmy,g_g,kfac,lhashf,otriplet)
c
c Eq. 25 (Furche and Ahlrichs) plus ground state gradient contributions.
c   As of 1-14-2013, the algorithm is:
c     1.  Form the XC-kernel (fxc) for evaluation of the XC-potential
c         (Vxc) gradient term.  This is reused for all excited states.
c     2.  Form the ground state density matrix (D).
c     3.  Convert the relaxed one-particle difference density matrix 
c         (P) in the AO basis.
c *** G.1 Gradient of the XC-potential (dVxc*P) 
c     4.  Form the excited state density matrix (D + P). 
c     5.  Convert (X+Y) and (X-Y) to the AO basis.
c     6.  Transform the energy-weighted difference density matrix (W) 
c         to the AO basis.
c *** G.T Alternative gradient of the XC-potential (dVxc*P)
c *** G.2 Gradient of the XC-energy (ground state contribution).
c *** G.3 Gradient of the core Hamiltonian (dh*P).
c *** G.4 Gradient of the XC-kernel (dfxc*(X+Y)*(X+Y)).
c *** G.5 Gradient of the two-electron integrals.
c *** G.6 Sum of the gradient contributions (Eq. 25 + ground state
c         simultaneously).
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "tddft_grad_util.fh"
#include "case.fh"
#include "cosmo.fh"
#include "cgridfile.fh"
c
c     Input:
c
      integer rtdb      ! the runtime database handle
      integer ihdl_geom ! the geometry handle
      integer ihdl_bfao ! the basis set handleA
      logical tda       ! True if Tamm-Dancoff approximation
      integer g_mo(2)   ! the molecular orbitals
      integer g_p(2)    ! one-particle difference density matrix in MO basis
      integer g_wp(2)   ! the energy weighted difference density matrix in MO basis
      integer g_xpy(2)  ! global arrays to hold X+Y
      integer g_xmy(2)  ! global arrays to hold X-Y
      integer nat       ! the number of atoms in the molecule
      integer nao       ! the number of atomic basis functions
      integer nfc(2)    ! the number of frozen cores
      integer naoc(2)   ! the number of active occupied orbitals
      integer nocc(2)   ! the number of occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nfv(2)    ! the number of frozen virtuals
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer nroots    ! the number of roots
      logical oskel     ! Use of symmetry in Fock build
c
      logical lhashf    ! =.true.  hybrid functionals
                        ! =.false. otherwise
      logical otriplet  ! =.true.  triplet excited states
                        ! =.false. singlet excited states
                        ! set to .false. for TDUDFT
c
      double precision kfac  ! the weight of the Hartree-Fock exchange
                             ! contributions
      double precision tol2e ! 2-electron integral tolerance
c
c     Input/Workspace:
c
      integer g_d(8)    ! On input g_d(1:2) are the handles for the
                        ! the ground state density matrices in AO basis
                        ! The entries g_d(3:8) are used to hold global
                        ! array handles for the transition density 
                        ! matrices, (X+Y) and (X-Y) matrices
                        ! respectively, all in AO basis. These matrices
                        ! are recalculated for every root.
c
      integer g_dtmp(4) ! Temporary storage for D and (X+Y)
c
c     Output:
c
      integer g_g       ! the global array for the gradient of 
                        ! the energy
c
c     Local:
c
      integer g_ovlp   ! overlap
      integer g_p_ao   ! the total density in AO basis
      integer g_wp_ao  ! the total weighted density matrix AO basis
      integer g_rhs_xc(2) ! the XC component of the der. Fock matrix
c
      integer idim(3)  ! the dimensions of the global arrays
      integer ichnk(3) ! the chunk sizes
c
      integer nproc    ! the number of processors
      integer iproc    ! the rank of this processor
c
      integer max1e       ! the maximum buffer size for 1-el integrals
      integer max2e       ! the maximum buffer size for 2-el integrals
      integer mscratch_1e ! the maximum scratch size for 1-el integrals
      integer mscratch_2e ! the maximum scratch size for 2-el integrals
      integer lbuf        ! the actual buffer size for integrals
      integer lscratch    ! the actual scratch size for integrals
      integer lforce      ! the size of the forces array
      integer max_at_bf   ! maximum number of basis functions on atoms
      integer lsqatom     ! the size of a max_at_bf square
c
      integer l_force,  k_force        ! handle and index for total force
      integer l_frc_nuc,k_frc_nuc      ! handle and index for nuclear force
      integer l_frc_kin,k_frc_kin      ! handle and index for kinetic force
      integer l_frc_wgh,k_frc_wgh      ! handle and index for weighted
                                       ! density force
      integer l_frc_2el,k_frc_2el      ! handle and index for 2-electron 
                                       ! force
      integer l_frc_2el_j, k_frc_2el_j ! handle and index for Coulomb
                                       ! integral force (CAM)
      integer l_frc_2el_k, k_frc_2el_k ! handle and index for Exchange
                                       ! integral force (CAM)
      integer l_frc_df0,k_frc_df0      ! handle and index for ground state 
                                       ! DFT force
      integer l_frc_df1,k_frc_df1      ! handle and index for transition 
                                       ! density TDDFT force
      integer l_frc_dft,k_frc_dft      ! handle and index for 3rd derivative
                                       ! TDDFT force
      integer calc_type                ! what do the DFT routines have to 
                                       ! compute? See subroutine dftgh_gridv0
                                       ! for definitions of valid values.
c
      integer l_buf,k_buf ! handle and index for the integral buffer
      integer l_scr,k_scr ! handle and index for the integral scratch
c
      integer l_dens, k_dens  ! handle and index for density buffer
      integer l_wdens,k_wdens ! handle and index for weighted density 
                              ! buffer
c
      integer ip, ir ! counters for spin components and roots resp.
c
      integer ilo(3),ihi(3) ! block limit for GA operations
      integer ld(3)         ! leading dimensions for local buffers
c
      integer blen    ! basis function blocking size for shell ordering
      integer maxblen ! maximum number of basis functions in group
      integer maxsh   ! maximum no. of shells in a group
      integer maxq    ! max quartets in a request
      parameter (maxblen = 36)
      parameter (maxsh   = 10)
      parameter (maxq    = maxsh**4)
      integer nopen   ! the number of singly occupied ROHF orbitals
      parameter (nopen   = 0)
c
      double precision jfac
      parameter (jfac=1.0d0) ! Hartree-Fock for now
c
      double precision dipmom(4) ! excited state dipole moment
      double precision cntr(3)   ! center (arbitrary)
      logical odipole            ! dipole flag
      double precision ecosmo
      character*255 cosmo_file,dmat_file
      logical odbug
c
      integer lsqa    ! the size of one square block
c
      integer nsh     ! the number of shells
c
      integer l_labels, k_labels
      integer l_list,   k_list
      integer l_q4,     k_q4
      integer l_shmap,  k_shmap
      integer l_shglo,  k_shglo,  l_shghi,  k_shghi
      integer l_shbflo, k_shbflo, l_shbfhi, k_shbfhi
      integer l_bfglo,  k_bfglo,  l_bfghi,  k_bfghi
      integer l_bfmap,  k_bfmap
      integer l_rbfmap, k_rbfmap
      integer l_bftoat, k_bftoat
      integer nshblocks
c
      integer ld_ija, kd_ija, ld_ijb, kd_ijb ! ground state density
      integer ld_ila, kd_ila, ld_ilb, kd_ilb ! ground state density
      integer ld_ika, kd_ika, ld_ikb, kd_ikb ! ground state density
      integer ld_kla, kd_kla, ld_klb, kd_klb ! ground state density
      integer ld_jla, kd_jla, ld_jlb, kd_jlb ! ground state density
      integer ld_jka, kd_jka, ld_jkb, kd_jkb ! ground state density
c
      integer lp_ija, kp_ija, lp_ijb, kp_ijb ! transition density
      integer lp_ila, kp_ila, lp_ilb, kp_ilb ! transition density
      integer lp_ika, kp_ika, lp_ikb, kp_ikb ! transition density
      integer lp_kla, kp_kla, lp_klb, kp_klb ! transition density
      integer lp_jla, kp_jla, lp_jlb, kp_jlb ! transition density
      integer lp_jka, kp_jka, lp_jkb, kp_jkb ! transition density
c
      integer lxpy_ija, kxpy_ija, lxpy_ijb, kxpy_ijb ! (X+Y)
      integer lxpy_kla, kxpy_kla, lxpy_klb, kxpy_klb ! (X+Y)
      integer lxpy_ila, kxpy_ila, lxpy_ilb, kxpy_ilb ! (X+Y)
      integer lxpy_jka, kxpy_jka, lxpy_jkb, kxpy_jkb ! (X+Y)
      integer lxpy_ika, kxpy_ika, lxpy_ikb, kxpy_ikb ! (X+Y)
      integer lxpy_jla, kxpy_jla, lxpy_jlb, kxpy_jlb ! (X+Y)
c
      integer lxmy_ila, kxmy_ila, lxmy_ilb, kxmy_ilb ! (X-Y)
      integer lxmy_jka, kxmy_jka, lxmy_jkb, kxmy_jkb ! (X-Y)
      integer lxmy_ika, kxmy_ika, lxmy_ikb, kxmy_ikb ! (X-Y)
      integer lxmy_jla, kxmy_jla, lxmy_jlb, kxmy_jlb ! (X-Y)
c
      integer l_pdm2d,  k_pdm2d
c
      integer l_act, k_act ! table for active atom administration
c     integer nactive      ! the number of active atoms
c
      integer l_cntoce    , k_cntoce     ! tables for quadrature
      integer l_cntobfr   , k_cntobfr    ! tables for quadrature
      integer l_cetobfr   , k_cetobfr    ! tables for quadrature
      integer l_rdens_atom, k_rdens_atom ! tables for quadrature
c
      integer i   ! arbitrary counter
      integer iat ! counter over atoms
c
      double precision pstrace
      double precision dum ! dummy variable to complete argument lists
      integer idum         ! dummy variable to complete argument lists
c
      character*8 scftype
      logical oroot
c
      integer msg_grad_nuc
      integer msg_grad_wgh
      integer msg_grad_kin
      integer msg_grad_2el
      integer msg_grad_df0
      integer msg_grad_df1
      integer msg_grad_dft
      integer msg_grad_2el_j ! cam j
      integer msg_grad_2el_k ! cam k
c
      parameter(msg_grad_nuc = 193)
      parameter(msg_grad_wgh = 194)
      parameter(msg_grad_kin = 195)
      parameter(msg_grad_2el = 196)
      parameter(msg_grad_df0 = 197)
      parameter(msg_grad_df1 = 198)
      parameter(msg_grad_dft = 199)
      parameter(msg_grad_2el_j = 200)
      parameter(msg_grad_2el_k = 201)
c
      character*32 pname
      logical status
c
c     Functions
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  grid_reopen
      external grid_reopen
      logical  xc_gotxc   ! do we have a density functional?
      external xc_gotxc
c
c     Preliminaries
      pname = "tddft_grad_compute_g: "
      odbug = util_print('cosmo', print_high)
      do i =1,3
        cntr(i) = 0.d0
      end do
      do i =1,4
        dipmom(i) = 0.d0
      end do
c
      scftype = ""
      if (ipol.eq.1) then
        scftype = "STDDFT"
        if (tda) scftype = "STDADFT"
        if (otriplet) then
          scftype = "TTDDFT"
          if (tda) scftype = "TTDADFT"
        endif
      else if (ipol.eq.2) then
        scftype = "UTDDFT"
        if (tda) scftype = "UTDADFT"
      else
        call errquit(pname//'illegal value of ipol',0,UERR)
      endif
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
      oroot = iproc.eq.0
c
      idim(1)  = nao
      idim(2)  = nao
      ichnk(1) = -1
      ichnk(2) = -1
      g_p_ao  = ga_create_atom_blocked(ihdl_geom,ihdl_bfao,"g_p_ao")
      g_wp_ao = ga_create_atom_blocked(ihdl_geom,ihdl_bfao,"g_wp_ao")
c
c     The (X-Y) vector has no contribution to the overall
c     gradients if a pure DFT calculation is done.
c
      do ip = 1, ipol
       g_d(0+ip)=ga_create_atom_blocked(ihdl_geom,ihdl_bfao,"g_d_ao")
       g_d(2+ip)=ga_create_atom_blocked(ihdl_geom,ihdl_bfao,"g_t_ao")
       g_d(4+ip)=ga_create_atom_blocked(ihdl_geom,ihdl_bfao,"g_xpy_ao")
       if (lhashf)
     1  g_d(6+ip)=ga_create_atom_blocked(ihdl_geom,ihdl_bfao,
     2                                   "g_xmy_ao")
      enddo
c
      if (xc_gotxc()) then
        do ip = 1, 2*ipol
          g_dtmp(ip) = ga_create_atom_blocked(ihdl_geom,ihdl_bfao,
     1                                        "g_dtmp_ao")
        enddo
      endif
      idim(1)  = 3*nat*ipol
      idim(2)  = nao
      idim(3)  = nao
      ichnk(1) = idim(1)
      ichnk(2) = -1
      ichnk(3) = -1
      if (.not. nga_create(MT_DBL, 3, idim, 'DFT CPKS RHS',ichnk,
     &   g_rhs_xc(1)))
     &  call errquit(pname//'could not create g_rhs_xc',555,GA_ERR)
      call ga_zero(g_rhs_xc(1))
c
c     Initialise the integral gradient code
c
      call int_init(rtdb,1,ihdl_bfao)
      call schwarz_init(ihdl_geom,ihdl_bfao)
      call int_terminate()
      call intd_init(rtdb,1,ihdl_bfao)
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      call intb_mem_2e4c(max2e, mscratch_2e) ! blocking algorithm
      max2e = max(max2e,1296*100)            ! 100 D quartets
      lbuf = max(max1e, max2e)
      lscratch = max(mscratch_1e, mscratch_2e)
c
c     Pick up some dimensions
c
      if (.not. bas_numcont(ihdl_bfao,nsh))
     +  call errquit(pname//'could not get nsh',0, BASIS_ERR)
c
c     Set up scratch arrays
c
      blen = min(nao,maxblen)
      lforce = nat * 3
      if (.not. ma_push_get(mt_dbl,lforce,'forces',l_force,k_force))
     +  call errquit(pname//'could not allocate l_force',lforce,MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'nuc forces', l_frc_nuc,
     +     k_frc_nuc)) call errquit
     +  (pname//'could not allocate l_frc_nuc',l_frc_nuc, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'kin forces', l_frc_kin,
     +     k_frc_kin)) call errquit
     +  (pname//'could not allocate l_frc_kin',l_frc_kin, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'wgh forces', l_frc_wgh,
     +     k_frc_wgh)) call errquit
     +  (pname//'could not allocate l_frc_wgh',l_frc_wgh, MA_ERR)
      if (.not.cam_exch) then
        if (.not. ma_push_get(mt_dbl, lforce, '2el forces', l_frc_2el,
     +       k_frc_2el)) call errquit
     +    (pname//'could not allocate l_frc_2el',l_frc_2el, MA_ERR)
      else
        if (.not. ma_push_get(mt_dbl, lforce, '2el forces', l_frc_2el_j,
     1      k_frc_2el_j)) call errquit
     2     (pname//'could not allocate l_frc_2el_j',lforce, MA_ERR)
        if (.not. ma_push_get(mt_dbl, lforce, '2el forces', 
     1      l_frc_2el_k, k_frc_2el_k)) call errquit
     2     (pname//'could not allocate l_frc_2el_k',lforce, MA_ERR)
      endif
      if (.not. ma_push_get(mt_dbl, lforce, 'TDDFT forces', l_frc_df0,
     +     k_frc_df0)) call errquit
     +  (pname//'could not allocate l_frc_df0',l_frc_df0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'TDDFT forces', l_frc_df1,
     +     k_frc_df1)) call errquit
     +  (pname//'could not allocate l_frc_df1',l_frc_df1, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'TDDFT forces', l_frc_dft,
     +     k_frc_dft)) call errquit
     +  (pname//'could not allocate l_frc_dft',l_frc_dft, MA_ERR)
c
      if (.not.ma_push_get(mt_dbl,lbuf,'deriv buffer',l_buf,k_buf))
     +  call errquit(pname//'could not allocate buffer',lbuf,MA_ERR)
      if (.not.ma_push_get(mt_dbl,lscratch,'deriv scratch',l_scr,k_scr))
     +  call errquit(pname//'scratch alloc failed',lscratch, MA_ERR)
c
c     Set up local density matrix blocks
c
      if (.not. bas_nbf_ce_max(ihdl_bfao,max_at_bf)) 
     + call errquit(pname//'could not get max_at_bf',0,BASIS_ERR)
      lsqatom = max_at_bf * max_at_bf
      if (.not.ma_push_get(mt_dbl,lsqatom,'local_density',
     +                       l_dens,k_dens)) 
     +  call errquit(pname//'could not allocate l_dens',lsqatom,MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqatom,'local_w_density',
     +                     l_wdens,k_wdens)) 
     +  call errquit(pname//'could not allocate l_wdens',lsqatom,MA_ERR)
c
c     Now set up local density matrix blocks for the 2-electron
c     contributions
c
      lsqa = blen*blen
      if (.not.ma_push_get(mt_dbl,lsqa,'d_ija',ld_ija,kd_ija)) call
     +  errquit(pname//'could not allocate d_ija',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'d_ila',ld_ila,kd_ila)) call
     +  errquit(pname//'could not allocate d_ila',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'d_ika',ld_ika,kd_ika)) call
     +  errquit(pname//'could not allocate d_ika',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'d_kla',ld_kla,kd_kla)) call
     +  errquit(pname//'could not allocate d_kla',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'d_jla',ld_jla,kd_jla)) call
     +  errquit(pname//'could not allocate d_jla',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'d_jka',ld_jka,kd_jka)) call
     +  errquit(pname//'could not allocate d_jka',0, MA_ERR)
c
      if (.not.ma_push_get(mt_dbl,lsqa,'p_ija',lp_ija,kp_ija)) call
     +  errquit(pname//'could not allocate p_ija',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'p_ila',lp_ila,kp_ila)) call
     +  errquit(pname//'could not allocate p_ila',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'p_ika',lp_ika,kp_ika)) call
     +  errquit(pname//'could not allocate p_ika',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'p_kla',lp_kla,kp_kla)) call
     +  errquit(pname//'could not allocate p_kla',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'p_jla',lp_jla,kp_jla)) call
     +  errquit(pname//'could not allocate p_jla',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'p_jka',lp_jka,kp_jka)) call
     +  errquit(pname//'could not allocate p_jka',0, MA_ERR)
c
      if (.not.ma_push_get(mt_dbl,lsqa,'xpy_ija',lxpy_ija,kxpy_ija))
     +  call errquit(pname//'could not allocate xpy_ija',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'xpy_kla',lxpy_kla,kxpy_kla))
     +  call errquit(pname//'could not allocate xpy_kla',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'xpy_ila',lxpy_ila,kxpy_ila))
     +  call errquit(pname//'could not allocate xpy_ila',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'xpy_jka',lxpy_jka,kxpy_jka))
     +  call errquit(pname//'could not allocate xpy_jka',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'xpy_ika',lxpy_ika,kxpy_ika))
     +  call errquit(pname//'could not allocate xpy_ika',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,lsqa,'xpy_jla',lxpy_jla,kxpy_jla))
     +  call errquit(pname//'could not allocate xpy_jla',0,MA_ERR)
c
      if (lhashf) then
        if (.not.ma_push_get(mt_dbl,lsqa,'xmy_ila',lxmy_ila,kxmy_ila))
     +    call errquit(pname//'could not allocate xmy_ila',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'xmy_jka',lxmy_jka,kxmy_jka))
     +    call errquit(pname//'could not allocate xmy_jka',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'xmy_ika',lxmy_ika,kxmy_ika))
     +    call errquit(pname//'could not allocate xmy_ika',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'xmy_jla',lxmy_jla,kxmy_jla))
     +    call errquit(pname//'could not allocate xmy_jla',0,MA_ERR)
      endif
c
      if (ipol.gt.1) then
        if (.not.ma_push_get(mt_dbl,lsqa,'d_ijb',ld_ijb,kd_ijb)) call
     +    errquit(pname//'could not allocate d_ijb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'d_ilb',ld_ilb,kd_ilb)) call
     +    errquit(pname//'could not allocate d_ilb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'d_ikb',ld_ikb,kd_ikb)) call
     +    errquit(pname//'could not allocate d_ikb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'d_klb',ld_klb,kd_klb)) call
     +    errquit(pname//'could not allocate d_klb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'d_jlb',ld_jlb,kd_jlb)) call
     +    errquit('tddft_grad_1e: could not allocate d_jlb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'d_jkb',ld_jkb,kd_jkb)) call
     +    errquit(pname//'could not allocate d_jkb',0, MA_ERR)
c
        if (.not.ma_push_get(mt_dbl,lsqa,'p_ijb',lp_ijb,kp_ijb)) call
     +    errquit(pname//'could not allocate p_ijb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'p_ilb',lp_ilb,kp_ilb)) call
     +    errquit(pname//'could not allocate p_ilb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'p_ikb',lp_ikb,kp_ikb)) call
     +    errquit(pname//'could not allocate p_ikb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'p_klb',lp_klb,kp_klb)) call
     +    errquit(pname//'could not allocate p_klb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'p_jlb',lp_jlb,kp_jlb)) call
     +    errquit(pname//'could not allocate p_jlb',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'p_jkb',lp_jkb,kp_jkb)) call
     +    errquit(pname//'could not allocate p_jkb',0, MA_ERR)
c
        if (.not.ma_push_get(mt_dbl,lsqa,'xpy_ijb',lxpy_ijb,kxpy_ijb))
     +    call errquit(pname//'could not allocate xpy_ijb',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'xpy_klb',lxpy_klb,kxpy_klb))
     +    call errquit(pname//'could not allocate xpy_klb',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'xpy_ilb',lxpy_ilb,kxpy_ilb))
     +    call errquit(pname//'could not allocate xpy_ilb',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'xpy_jkb',lxpy_jkb,kxpy_jkb))
     +    call errquit(pname//'could not allocate xpy_jkb',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'xpy_ikb',lxpy_ikb,kxpy_ikb))
     +    call errquit(pname//'could not allocate xpy_ikb',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,lsqa,'xpy_jlb',lxpy_jlb,kxpy_jlb))
     +    call errquit(pname//'could not allocate xpy_jlb',0,MA_ERR)
c
        if (lhashf) then
          if (.not.ma_push_get(mt_dbl,lsqa,'xmy_ilb',lxmy_ilb,kxmy_ilb))
     +     call errquit(pname//'could not allocate xmy_ilb',lsqa,MA_ERR)
          if (.not.ma_push_get(mt_dbl,lsqa,'xmy_jkb',lxmy_jkb,kxmy_jkb))
     +     call errquit(pname//'could not allocate xmy_jkb',lsqa,MA_ERR)
          if (.not.ma_push_get(mt_dbl,lsqa,'xmy_ikb',lxmy_ikb,kxmy_ikb))
     +     call errquit(pname//'could not allocate xmy_ikb',lsqa,MA_ERR)
          if (.not.ma_push_get(mt_dbl,lsqa,'xmy_jlb',lxmy_jlb,kxmy_jlb))
     +     call errquit(pname//'could not allocate xmy_jlb',lsqa,MA_ERR)
        endif
      endif
c
      if (.not.ma_push_get(mt_int,lbuf/3,'labels',l_labels,k_labels))
     +  call errquit(pname//'could not allocate labels',lbuf/3,MA_ERR)
      if (.not.ma_push_get(mt_int,4*maxq,'list',l_list,k_list))
     +  call errquit(pname//'could not allocate list',4*maxq,MA_ERR)
      if (.not.ma_push_get(mt_dbl,maxq,'q4',l_q4,k_q4))
     +  call errquit(pname//'could not allocate q4',maxq,MA_ERR)
      if (.not.ma_push_get(mt_log,nat,'active atoms',l_act,k_act))
     +  call errquit(pname//'could not allocate act',nat,MA_ERR)
c
      if (.not. ma_push_get(mt_int, nsh, 'shmap', l_shmap, k_shmap))
     +     call errquit(pname//'could not allocate shmap',nsh,MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shglo', l_shglo, k_shglo))
     +     call errquit(pname//'could not allocate blo',nsh,MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shghi', l_shghi, k_shghi))
     +     call errquit(pname//'could not allocate bhi',nsh,MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbflo', l_shbflo, k_shbflo))
     +     call errquit(pname//'could not allocate bflo',nsh,MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbfhi', l_shbfhi, k_shbfhi))
     +     call errquit(pname//'could not allocate bfhi',nsh,MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'bfglo', l_bfglo, k_bfglo))
     +     call errquit(pname//'could not allocate blo',nao,MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'bfghi', l_bfghi, k_bfghi))
     +     call errquit(pname//'could not allocate bhi',nao,MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'bfmap', l_bfmap, k_bfmap))
     +     call errquit(pname//'could not allocate bfmap',nao,MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'rbfmap', l_rbfmap, k_rbfmap))
     +     call errquit(pname//'could not allocate rbfmap',nao,MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'bftoat', l_bftoat, k_bftoat))
     +     call errquit(pname//'could not allocate bftoat',nao,MA_ERR)
c
      if (.not. ma_push_get(mt_dbl, blen**4, 'pdm2d',
     +     l_pdm2d, k_pdm2d)) 
     +  call errquit(pname//'failed allocating pdm2d',blen**4, MA_ERR)
c
      call grad_shorder(ihdl_bfao, nsh, nao, maxsh, blen,
     +     nshblocks, int_mb(k_shglo), int_mb(k_shghi),
     +     int_mb(k_shmap),
     +     int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bfglo),
     +     int_mb(k_bfghi), int_mb(k_shbflo), int_mb(k_shbfhi))

      do i = 1, nao
       if (bas_bf2ce(ihdl_bfao,i,iat)) then
          int_mb(k_bftoat+i-1) = iat
       else
          call errquit(pname//'bas_bf2ce failed',i,BASIS_ERR)
       end if
      end do
c
c     Make all atoms active for the DFT contributions and the 
c     2-electron contributions.
c
      do i = 0, nat-1
        log_mb(k_act+i) = .true.
      enddo
c
c 1.  Compute the XC components for the derivative Kohn-Sham matrix.
c     These are needed for every excited state.
c
      if (xc_gotxc()) then
c
c Daniel (12/1/12): This part calculates the XC-functional 2nd
c derivatives for the relaxed one-particle difference density
c contribution (known because of calc_type = 3).
c
        calc_type = 3
        call tddft_grad_compute_dao(ipol,nao,nocc,g_mo,g_d)
        if (ipol.eq.1) call ga_scale(g_d(1),2.0d0)
        call xc_init_index(ihdl_geom,ipol,nao,ihdl_bfao,g_d,
     &       l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &       l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom)
c
         if(.not.grid_written) then
            grid_written=grid_reopen(ihdl_geom)
         endif
c
c Daniel (1-14-13): This call gives the gradient of the XC-potential,
c which is the XC-kernel multiplied by the gradient of the ground state
c electron density.  After implementing this term in an alternative way,
c I have confirmed that the quadrature weight gradients are built in 
c the CPKS RHS routines.  The specific location is still puzzling. 
c
c Daniel (2-23-13): It turns out that the triplet flag is not needed
c for this gradient term in restricted triplet gradients.
        call dftgh_gridv0(rtdb,ihdl_geom,ihdl_bfao, ipol,nao,
     &       g_d, dum, dum, g_rhs_xc, calc_type,
     &       nat, log_mb(k_act), nat, dbl_mb(k_rdens_atom),
     &       int_mb(k_cetobfr), 0, .false., 0d0, .false.)
        call xc_exit_index(l_cntoce, l_cntobfr, l_cetobfr,
     &                     l_rdens_atom)
      endif
c
      do ir = 1, nroots
c
c       Clear the scratch arrays
c
        call dfill(lforce, 0.0D0, dbl_mb(k_force), 1)
        call dfill(lforce, 0.0D0, dbl_mb(k_frc_nuc), 1)
        call dfill(lforce, 0.0D0, dbl_mb(k_frc_kin), 1)
        call dfill(lforce, 0.0D0, dbl_mb(k_frc_wgh), 1)
        call dfill(lforce, 0.0D0, dbl_mb(k_frc_df0), 1)
        call dfill(lforce, 0.0D0, dbl_mb(k_frc_df1), 1)
        call dfill(lforce, 0.0D0, dbl_mb(k_frc_dft), 1)
        if (.not.cam_exch) then
          call dfill(lforce, 0.0D0, dbl_mb(k_frc_2el), 1)
        else
          call dfill(lforce, 0.0D0, dbl_mb(k_frc_2el_j), 1)
          call dfill(lforce, 0.0D0, dbl_mb(k_frc_2el_k), 1)
        endif
c
c 2.    Compute the ground state density matrix (this gets messed up
c       when the basis functions are reordered for the 2-electron 
c       contributions, so we have to recalculate it every time).
c
        call tddft_grad_compute_dao(ipol,nao,nocc,g_mo,g_d)
c DEBUG
c        call ga_print(g_d)
        if (util_print('tddft grad g',print_debug)) then
          if (oroot) write(LuOut,*)'DEBUG: '//pname//'D'
          call tddft_grad_print_array(ipol,1,g_d,dble(ipol))
        endif
c DEBUG
c
c
c 3.    Transform the transition density matrix adding alpha and beta 
c       spin components together
c
        call ga_zero(g_p_ao)
        do ip = 1, ipol
          call tddft_grad_trans_mo2ao(1,nao,nfc(ip),naoc(ip),nocc(ip),
     +          nav(ip),nfv(ip),ir,1.0d0,0.0d0,"pq",g_mo(ip),g_p(ip),
     +          "pq",g_d(2+ip))
c
c         Add the alpha and beta relaxed one-particle difference density 
c         matrices together to get the total relaxed one-particle 
c         difference density matrix.
c
          call ga_add(1.0d0,g_d(2+ip),1.0d0,g_p_ao,g_p_ao)
        enddo
c DEBUG
c        call ga_print(g_p_ao)
        if (util_print('tddft grad g',print_debug)) then
          if (oroot) write(LuOut,*)'DEBUG: '//pname//'P'
c          call tddft_grad_print_array(ipol,1,g_p_ao,dble(ipol))
          call tddft_grad_print_array(ipol,1,g_d(3),dble(ipol))
        endif
c DEBUG
c
c G.1   Form the gradient of Vxc multiplied by the relaxed one-particle
c       difference density matrix   
c
c Daniel (1/14/13): This is the XC-potential gradient term (i.e. dVxc*P).
        if (xc_gotxc()) then
          do ip = 1, ipol
            do i = 1, 3*nat
              ilo(1) = i+(ip-1)*3*nat
              ihi(1) = i+(ip-1)*3*nat
              ilo(2) = 1
              ilo(3) = 1
              ihi(2) = nao
              ihi(3) = nao
              dbl_mb(k_frc_df1+i-1) = dbl_mb(k_frc_df1+i-1)
     &        +          nga_ddot_patch(g_rhs_xc(1),'N',ilo,ihi,
     &                                  g_d(2+ip),'N',ilo(2),ihi(2))/
     &                                  nproc
            enddo
          enddo
        endif
c
c 4.    Add the ground state density matrix to the relaxed one-particle
c       difference density matrix.  
c
c Here, g_p_ao contains the total density matrix (alpha plus beta), 
c while g_d(3) (and g_d(4)) contain the alpha and beta density matrices 
c separately.  In both cases we get the excited state density matrix.
c
c       Calculate and write out the excited-state dipole moment
c
        call tddft_grad_compute_dao(ipol,nao,nocc,g_mo,g_d) ! calculate ground state density
        if (ipol.eq.1) call ga_scale(g_d(1),2.0d0) ! scale
        call tddft_grad_mpole(rtdb,ihdl_bfao,ipol,g_d(1),g_d(2),g_p_ao) ! multipole analysis
c
c       Write difference density matrix in the ao basis
c
        call util_file_name('diffdmat',.false.,.false.,dmat_file) ! get filename
        call ao_1prdm_write(nao,g_p_ao,dmat_file)
c
c       Add the ground state density to the difference density to get the 
c       total excited state relaxed density
c
        call tddft_grad_compute_dao(ipol,nao,nocc,g_mo,g_d) ! calculate ground state density
        if (ipol.eq.1) then
          call ga_add(2.0d0,g_d(1),1.0d0,g_p_ao,g_p_ao)
          call ga_add(1.0d0,g_d(1),1.0d0,g_d(3),g_d(3))
        else
          call ga_add(1.0d0,g_d(1),1.0d0,g_p_ao,g_p_ao)
          call ga_add(1.0d0,g_d(2),1.0d0,g_p_ao,g_p_ao)
c Skip scaling here so that we can output alpha and beta 
c excited state densities
c          call ga_scale(g_d(1),0.5d0)
c          call ga_scale(g_d(2),0.5d0)
          call ga_add(1.0d0,g_d(1),1.0d0,g_d(3),g_d(3))
          call ga_add(1.0d0,g_d(2),1.0d0,g_d(4),g_d(4))
        endif
c
c       Write total excited-state density in the ao basis 
c       (g_p_ao = total excited state density)
c
        call util_file_name('dmat',.false.,.false.,dmat_file) ! get filename
        call ao_1prdm_write(nao,g_p_ao,dmat_file)
c
c       Write excited-state density components in the ao basis
c
        if (ipol.gt.1) then
          call util_file_name('dmatA',.false.,.false.,dmat_file) ! get filename
          call ao_1prdm_write(nao,g_d(3),dmat_file)
          call util_file_name('dmatB',.false.,.false.,dmat_file) ! get filename
          call ao_1prdm_write(nao,g_d(4),dmat_file)
c Now scale g_d(1) and g_d(2) and subtract the result from
c g_d(3) and g_d(4) so that the gradients still work out correctly
          call ga_scale(g_d(1),0.5d0)
          call ga_scale(g_d(2),0.5d0)
          call ga_add(-1.d0,g_d(1),1.0d0,g_d(3),g_d(3))
          call ga_add(-1.d0,g_d(2),1.0d0,g_d(4),g_d(4))
        end if
c
c       Calculate the trace of P*S: Check normalization. Total number of electrons
c
        if (.not.ga_create(mt_dbl,nao,nao,'AO overlap',-1,-1,g_ovlp))
     +       call errquit(pname//'failed to create g_ovlp',0, GA_ERR)
        call ga_zero(g_ovlp)
        call int_1e_ga(ihdl_bfao,ihdl_bfao,g_ovlp,'overlap',oskel)
         if (oskel) call sym_symmetrize(ihdl_geom,ihdl_bfao,.false.,
     &          g_ovlp)
         pstrace=ga_ddot(g_p_ao,g_ovlp)
         if(ga_nodeid().eq.0) then
           write (luout,*)
           write (luout,'(5x,a,1x,e15.7)') 'No. of electrons (tr(P*S)): 
     & ',pstrace
           write (luout,*)
         end if
c
c        Calculate the cosmo charges using the excited-state density
c
c        if (oroot) write(luout,*) "cosmo_on:" ,cosmo_on
         if (cosmo_on) then
           cosmo_file = "cosmo-exci.xyz"
           ecosmo = 0.d0
c          if (oroot) write(luout,*) "calling cosmo_charges_from_dmat"
           call cosmo_charges_from_dmat(rtdb,
     &                                  ihdl_bfao,
     &                                  ihdl_geom,
     &                                  ecosmo, 
     &                                  odbug,       ! debug
     &                                  1,           ! up and down are combined
     &                                  g_p_ao,      ! input density
     &                                  cosmo_file)  ! cosmo charges file name
c
c          if (oroot) write(luout,*) "ecosmo:" ,ecosmo
         end if ! cosmo check
c
c        end if ! odipole
c
c 5.    Transform (X+Y) and (X-Y) from MO to AO basis
c
        do ip = 1, ipol
          call tddft_grad_trans_mo2ao(1,nao,nfc(ip),naoc(ip),nocc(ip),
     +         nav(ip),nfv(ip),ir,1.0d0,0.0d0,"ib",g_mo(ip),g_xpy(ip),
     +         "ib",g_d(4+ip))
c
c         tddft_grad_grad_2e assumes symmetric quantities
c
          call ga_symmetrize(g_d(4+ip))
c
c Daniel (1-7-13): For calculations with exact exchange, the X vector
c has symmetric and anti-symmetric parts when doing CIS.  We need to
c create the "(X-Y)" density matrix as a result, even for CIS.  This
c part can still be skipped for pure DFT calculations, where the X
c vector is symmetric. 
          if (lhashf) then
            if (.not.tda) then
              call tddft_grad_trans_mo2ao(1,nao,nfc(ip),naoc(ip),
     +             nocc(ip),nav(ip),nfv(ip),ir,1.0d0,0.0d0,"ib",
     +             g_mo(ip),g_xmy(ip),"ib",g_d(6+ip))
            else
              call tddft_grad_trans_mo2ao(1,nao,nfc(ip),naoc(ip),
     +             nocc(ip),nav(ip),nfv(ip),ir,1.0d0,0.0d0,"ib",
     +             g_mo(ip),g_xpy(ip),"ib",g_d(6+ip))
            endif
c
c           tddft_grad_grad_2e assumes symmetric quantities
c
            call tddft_grad_anti_symmetrize(g_d(6+ip))
c
          endif
        enddo
c Daniel (11-30-12): Print X+Y and X-Y in the AO basis
c DEBUG
c        call ga_print(g_d(5))
        if (util_print('tddft grad g',print_debug)) then
          if (oroot) write(LuOut,*)'DEBUG: '//pname//'(X+Y)'
          call tddft_grad_print_array(ipol,1,g_d(5),dble(ipol))
        endif
c        if (lhashf) call ga_print(g_d(7))
c DEBUG
c
c 6.    Transform the weighted density matrix adding alpha and beta spin
c       components together
c
        call ga_zero(g_wp_ao)
        do ip = 1, ipol
          call tddft_grad_trans_mo2ao(1,nao,nfc(ip),naoc(ip),nocc(ip),
     +         nav(ip),nfv(ip),ir,1.0d0,1.0d0,"pq",g_mo(ip),g_wp(ip),
     +          "pq",g_wp_ao)
        enddo
c Daniel (11-30-12): Print W (note, this contains the ground state
c density matrix contribution as well).
c DEBUG
c        call ga_print(g_wp_ao)
        if (util_print('tddft grad g',print_debug)) then
          if (oroot) write(LuOut,*)'DEBUG: '//pname//'W'
          call tddft_grad_print_array(1,1,g_wp_ao,dble(ipol))
        endif
c DEBUG
c
c G.2   Calculate the XC-energy functional gradient contributions
c
c Daniel (1/14/13): XC-energy gradient term from the ground state 
c energy.  This part is identical to the ground state gradient code,
c and therefore includes the quadrature weight gradients as well.
        if (xc_gotxc()) then
          calc_type = 1
          if (ipol.eq.1) then
            call ga_scale(g_d(1),2.0d0)
          else
            call ga_scale(g_d(1),2.0d0)
            call ga_scale(g_d(2),2.0d0)
          endif
          call xc_init_index(ihdl_geom,ipol,nao,ihdl_bfao,g_d,
     &         l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &         l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom)
c
          if(.not.grid_written) then
            grid_written=grid_reopen(ihdl_geom)
          endif
c
          call dftgh_gridv0(rtdb,ihdl_geom,ihdl_bfao, ipol,nao,
     &         g_d, 
     &         dbl_mb(k_frc_df0), dum, idum, calc_type,
     &         nat, log_mb(k_act), nat, dbl_mb(k_rdens_atom),
     &         int_mb(k_cetobfr), 0, otriplet, 0d0, .false.)
          call xc_exit_index(l_cntoce, l_cntobfr, l_cetobfr,
     &                       l_rdens_atom)
          if (ipol.eq.1) then
            call ga_scale(g_d(1),0.5d0)
c           call ga_scale(g_d(3),0.5d0)
          else
            call ga_scale(g_d(1),0.5d0)
            call ga_scale(g_d(2),0.5d0)
          endif
        endif ! xc_gotxc()
c
c G.3   Calculate the 1-electron contributions
c
c       - k_frc_nuc : will be set to the gradient of the nuclear
c                     repulsion energy.
        call tddft_grad_grad_1e(dbl_mb(k_buf),lbuf,dbl_mb(k_scr),
     +    lscratch,dbl_mb(k_dens), dbl_mb(k_wdens), dbl_mb(k_frc_nuc),
     +    dbl_mb(k_frc_kin), dbl_mb(k_frc_wgh), 0,
     +    g_p_ao, g_wp_ao, ihdl_bfao, ihdl_geom, nproc, nat, 
     +    max_at_bf, rtdb, oskel)
c
c G.4   Calculate the TDDFT excitation contributions
c
c Daniel (1-14-13): Calculate the XC-kernel gradient term.  This
c part needs work because calc_type 4 doesn't exist yet.  It may
c be a good idea to make this part behave like the XC-potential
c gradient in the future so that we don't need to recalculate
c the XC-third derivatives repeatedly for users that request 
c multiple roots.
c Daniel (1-16-13): We need to do quadrature weight gradients here
c for the XC-kernel contribution.  I'm going to take the "brute force"
c approach here, but there is probably a clever way of incorporating
c these...
        if (xc_gotxc()) then
          calc_type = 4
          if (ipol.eq.1) then
            call ga_scale(g_d(1),2.0d0)
            call ga_scale(g_d(5),2.0d0)
          else
            call ga_scale(g_d(1),2.0d0)
            call ga_scale(g_d(2),2.0d0)
c For unrestricted jobs, the total perturbed density is given as:
c
c prho_tot = -1/sqrt(2)*(prho_alpha + prho_beta)
c 
c Keep this in mind if you ever compare the values of the densities.
          endif
c Daniel (1-16-13): Density matrices:
c g_d(1,2) => GS density matrix (D)
c g_d(3,4) => Relaxed one-particle difference density matrix (P)
c g_d(5,6) => (X+Y) density matrix
c g_d(7,8) => (X-Y) density matrix
c Daniel (1-17-13): Here we reorder the excited state density matrices
c for xc_rho_gen.  There is probably a more memory efficient way to do 
c this. 
          if (ipol.eq.1) then
            call ga_copy(g_d(1), g_dtmp(2)) ! D in g_dtmp(2)
            call ga_copy(g_d(5), g_dtmp(1)) ! (X+Y) in g_dtmp(1)
          else
            call ga_copy(g_d(2), g_dtmp(4)) ! g_dtmp(4) stores D_beta
            call ga_copy(g_d(1), g_dtmp(2)) ! g_dtmp(2) stores D_alpha
            call ga_copy(g_d(5), g_dtmp(1)) ! g_dtmp(1) stores (X+Y)_alpha
            call ga_copy(g_d(6), g_dtmp(3)) ! g_dtmp(3) stores (X+Y)_beta
          endif
          call ga_sync()
          call xc_init_index(ihdl_geom,ipol,nao,ihdl_bfao,g_dtmp,
     &         l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &         l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom)
c
          if(.not.grid_written) then
            grid_written=grid_reopen(ihdl_geom)
          endif
c
c Daniel (1-16-13): Here I modified the gradient routine to include
c the number of perturbed density matrices.  This is necessary so that we
c can define a value for looping in xc_3rd_deriv and build the
c perturbed density in xc_rhogen.  Since we loop over the roots, 
c the number of Fock matrices is always equal to 1.
          call dftgh_gridv0(rtdb,ihdl_geom,ihdl_bfao, ipol,nao,
     &         g_dtmp, dbl_mb(k_frc_dft), dum, idum, calc_type,
     &         nat, log_mb(k_act), nat, dbl_mb(k_rdens_atom),
     &         int_mb(k_cetobfr), 1, otriplet, 0d0, .false.)
c          call dftgh_gridv0(rtdb,ihdl_geom,ihdl_bfao, ipol,nao,
c     &         g_dtmp, dbl_mb(k_frc_dft), dum, idum, calc_type,
c     &         nat, log_mb(k_act), nat, dbl_mb(k_rdens_atom),
c     &         int_mb(k_cetobfr), ipol, otriplet)
          call xc_exit_index(l_cntoce, l_cntobfr, l_cetobfr,
     &                       l_rdens_atom)
c Daniel (1-16-13): Here we put D back into its proper location for
c the 2-electron part.  For unrestricted jobs, we recover P from 
c g_p_ao and g_wp_ao in unrestricted jobs.
c          if (ipol.eq.1) then
c            g_d(1) = g_d(2)  ! g_d(1) stores D
c          else
c            g_d(1) = g_d(2)  ! g_d(1) stores D_alpha
c            g_d(2) = g_d(4)  ! g_d(2) stores D_beta
c            g_d(3) = g_p_ao  ! g_d(3) stores P_alpha
c            g_d(4) = g_wp_ao ! g_d(4) stroes P_beta
c          endif
          if (ipol.eq.1) then
            call ga_scale(g_d(1),0.5d0)
            call ga_scale(g_d(5),0.5d0)
          else
            call ga_scale(g_d(1),0.5d0)
            call ga_scale(g_d(2),0.5d0)
          endif
        endif  ! xc_gotxc()
c
c G.5   Calculate the 2-electron contributions
c
        do ip = 1, ipol
          call ga_reorder(g_d(0+ip), .true., int_mb(k_rbfmap),
     +                    .true., int_mb(k_rbfmap))
          call ga_reorder(g_d(2+ip), .true., int_mb(k_rbfmap),
     +                    .true., int_mb(k_rbfmap))
          call ga_reorder(g_d(4+ip), .true., int_mb(k_rbfmap),
     +                    .true., int_mb(k_rbfmap))
c Daniel (1-7-13): Need to have 2-electron contributions for (X-Y), even
c if CIS/TDA is used.  This can be avoided if pure DFT calculations are
c done.
          if (lhashf) then
            call ga_reorder(g_d(6+ip), .true., int_mb(k_rbfmap),
     +                      .true., int_mb(k_rbfmap))
          end if
        end do
c
c       Eq. 27 (Furche & Ahlrichs)
c       tddft_grad_get_dens.F & tddft_grad_make_2pdm.F are called within
c       we should have another version of this function
c Daniel (4-8-13): Need to have conditional statement here for CAM
c functionals.
        if (.not.cam_exch) then
          call tddft_grad_grad_2e( 
     +         dbl_mb(kd_ija), dbl_mb(kd_ijb), 
     +         dbl_mb(kd_ila), dbl_mb(kd_ilb), 
     +         dbl_mb(kd_ika), dbl_mb(kd_ikb),
     +         dbl_mb(kd_kla), dbl_mb(kd_klb),
     +         dbl_mb(kd_jla), dbl_mb(kd_jlb),
     +         dbl_mb(kd_jka), dbl_mb(kd_jkb),
     +         dbl_mb(kp_ija), dbl_mb(kp_ijb),      ! p matrix
     +         dbl_mb(kp_ila), dbl_mb(kp_ilb),
     +         dbl_mb(kp_ika), dbl_mb(kp_ikb),
     +         dbl_mb(kp_kla), dbl_mb(kp_klb),
     +         dbl_mb(kp_jla), dbl_mb(kp_jlb),
     +         dbl_mb(kp_jka), dbl_mb(kp_jkb),
     +         dbl_mb(kxpy_ija), dbl_mb(kxpy_ijb),  ! arrays with x+y
     +         dbl_mb(kxpy_kla), dbl_mb(kxpy_klb),
     +         dbl_mb(kxpy_ila), dbl_mb(kxpy_ilb),
     +         dbl_mb(kxpy_jka), dbl_mb(kxpy_jkb),
     +         dbl_mb(kxpy_ika), dbl_mb(kxpy_ikb),
     +         dbl_mb(kxpy_jla), dbl_mb(kxpy_jlb),
     +         dbl_mb(kxmy_ila), dbl_mb(kxmy_ilb),  ! arrays with x-y
     +         dbl_mb(kxmy_jka), dbl_mb(kxmy_jkb), 
     +         dbl_mb(kxmy_ika), dbl_mb(kxmy_ikb),
     +         dbl_mb(kxmy_jla), dbl_mb(kxmy_jlb), 
     +         dbl_mb(k_frc_2el), g_d, 0, blen, ihdl_geom, ihdl_bfao,
     +         nproc, nat, lscratch, dbl_mb(k_scr), 
     +         lbuf/12, dbl_mb(k_buf), int_mb(k_labels), maxq, 
     +         int_mb(k_list), dbl_mb(k_q4), tol2e, nsh, log_mb(k_act),
     +         oskel, scftype, .false., nopen, nao, 
     +         dum, dum, dum,  dum, dbl_mb(k_pdm2d), dum,
     +         nshblocks, int_mb(k_shmap), int_mb(k_shglo), 
     +         int_mb(k_shghi), int_mb(k_bfglo), int_mb(k_bfghi),
     +         int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bftoat),
     +         int_mb(k_shbflo), int_mb(k_shbfhi),jfac,kfac,.false.)
        else
c Coulomb part
          call case_setflags(.false.)  ! No screening for the J part
c
          call tddft_grad_grad_2e(
     +         dbl_mb(kd_ija), dbl_mb(kd_ijb),
     +         dbl_mb(kd_ila), dbl_mb(kd_ilb),
     +         dbl_mb(kd_ika), dbl_mb(kd_ikb),
     +         dbl_mb(kd_kla), dbl_mb(kd_klb),
     +         dbl_mb(kd_jla), dbl_mb(kd_jlb),
     +         dbl_mb(kd_jka), dbl_mb(kd_jkb),
     +         dbl_mb(kp_ija), dbl_mb(kp_ijb),      ! p matrix
     +         dbl_mb(kp_ila), dbl_mb(kp_ilb),
     +         dbl_mb(kp_ika), dbl_mb(kp_ikb),
     +         dbl_mb(kp_kla), dbl_mb(kp_klb),
     +         dbl_mb(kp_jla), dbl_mb(kp_jlb),
     +         dbl_mb(kp_jka), dbl_mb(kp_jkb),
     +         dbl_mb(kxpy_ija), dbl_mb(kxpy_ijb),  ! arrays with x+y
     +         dbl_mb(kxpy_kla), dbl_mb(kxpy_klb),
     +         dbl_mb(kxpy_ila), dbl_mb(kxpy_ilb),
     +         dbl_mb(kxpy_jka), dbl_mb(kxpy_jkb),
     +         dbl_mb(kxpy_ika), dbl_mb(kxpy_ikb),
     +         dbl_mb(kxpy_jla), dbl_mb(kxpy_jlb),
     +         dbl_mb(kxmy_ila), dbl_mb(kxmy_ilb),  ! arrays with x-y
     +         dbl_mb(kxmy_jka), dbl_mb(kxmy_jkb),
     +         dbl_mb(kxmy_ika), dbl_mb(kxmy_ikb),
     +         dbl_mb(kxmy_jla), dbl_mb(kxmy_jlb),
     +         dbl_mb(k_frc_2el_j), g_d, 0, blen, ihdl_geom, ihdl_bfao,
     +         nproc, nat, lscratch, dbl_mb(k_scr),
     +         lbuf/12, dbl_mb(k_buf), int_mb(k_labels), maxq,
     +         int_mb(k_list), dbl_mb(k_q4), tol2e, nsh, log_mb(k_act),
     +         oskel, scftype, .false., nopen, nao,
     +         dum, dum, dum,  dum, dbl_mb(k_pdm2d), dum,
     +         nshblocks, int_mb(k_shmap), int_mb(k_shglo),
     +         int_mb(k_shghi), int_mb(k_bfglo), int_mb(k_bfghi),
     +         int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bftoat),
     +         int_mb(k_shbflo), int_mb(k_shbfhi),jfac,0.0d0,.false.)
c Exchange part
          call case_setflags(.true.)  ! Screen the K part 
c
          call tddft_grad_grad_2e(
     +         dbl_mb(kd_ija), dbl_mb(kd_ijb),
     +         dbl_mb(kd_ila), dbl_mb(kd_ilb),
     +         dbl_mb(kd_ika), dbl_mb(kd_ikb),
     +         dbl_mb(kd_kla), dbl_mb(kd_klb),
     +         dbl_mb(kd_jla), dbl_mb(kd_jlb),
     +         dbl_mb(kd_jka), dbl_mb(kd_jkb),
     +         dbl_mb(kp_ija), dbl_mb(kp_ijb),      ! p matrix
     +         dbl_mb(kp_ila), dbl_mb(kp_ilb),
     +         dbl_mb(kp_ika), dbl_mb(kp_ikb),
     +         dbl_mb(kp_kla), dbl_mb(kp_klb),
     +         dbl_mb(kp_jla), dbl_mb(kp_jlb),
     +         dbl_mb(kp_jka), dbl_mb(kp_jkb),
     +         dbl_mb(kxpy_ija), dbl_mb(kxpy_ijb),  ! arrays with x+y
     +         dbl_mb(kxpy_kla), dbl_mb(kxpy_klb),
     +         dbl_mb(kxpy_ila), dbl_mb(kxpy_ilb),
     +         dbl_mb(kxpy_jka), dbl_mb(kxpy_jkb),
     +         dbl_mb(kxpy_ika), dbl_mb(kxpy_ikb),
     +         dbl_mb(kxpy_jla), dbl_mb(kxpy_jlb),
     +         dbl_mb(kxmy_ila), dbl_mb(kxmy_ilb),  ! arrays with x-y
     +         dbl_mb(kxmy_jka), dbl_mb(kxmy_jkb),
     +         dbl_mb(kxmy_ika), dbl_mb(kxmy_ikb),
     +         dbl_mb(kxmy_jla), dbl_mb(kxmy_jlb),
     +         dbl_mb(k_frc_2el_k), g_d, 0, blen, ihdl_geom, ihdl_bfao,
     +         nproc, nat, lscratch, dbl_mb(k_scr),
     +         lbuf/12, dbl_mb(k_buf), int_mb(k_labels), maxq,
     +         int_mb(k_list), dbl_mb(k_q4), tol2e, nsh, log_mb(k_act),
     +         oskel, scftype, .false., nopen, nao,
     +         dum, dum, dum,  dum, dbl_mb(k_pdm2d), dum,
     +         nshblocks, int_mb(k_shmap), int_mb(k_shglo),
     +         int_mb(k_shghi), int_mb(k_bfglo), int_mb(k_bfghi),
     +         int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bftoat),
     +         int_mb(k_shbflo), int_mb(k_shbfhi),0.0d0,kfac,.false.)
c
          call case_setflags(.false.)
        endif
c
c G.6   Accumulate the contributions across all processors
c
        call ga_sync()
        call ga_dgop(msg_grad_nuc, dbl_mb(k_frc_nuc), lforce, '+')
        call ga_dgop(msg_grad_wgh, dbl_mb(k_frc_wgh), lforce, '+')
        call ga_dgop(msg_grad_kin, dbl_mb(k_frc_kin), lforce, '+')
c Daniel (4-8-13): Need conditional statement here for CAM functionals
        if (.not.cam_exch) then
          call ga_dgop(msg_grad_2el, dbl_mb(k_frc_2el), lforce, '+')
        else
          call ga_dgop(msg_grad_2el_j, dbl_mb(k_frc_2el_j), lforce, '+')
          call ga_dgop(msg_grad_2el_k, dbl_mb(k_frc_2el_k), lforce, '+')
        endif
        call ga_dgop(msg_grad_df0, dbl_mb(k_frc_df0), lforce, '+')
        call ga_dgop(msg_grad_df1, dbl_mb(k_frc_df1), lforce, '+')
        call ga_dgop(msg_grad_dft, dbl_mb(k_frc_dft), lforce, '+')
        call ga_sync()
c
c       Add the gradient contributions to the global array
c
        call daxpy(lforce,1.0d0,dbl_mb(k_frc_nuc),1,dbl_mb(k_force),1)
        call daxpy(lforce,1.0d0,dbl_mb(k_frc_kin),1,dbl_mb(k_force),1)
        call daxpy(lforce,1.0d0,dbl_mb(k_frc_wgh),1,dbl_mb(k_force),1)
c Daniel (4-8-13): Need conditional statement here for CAM functionals
        if (.not.cam_exch) then
          call daxpy(lforce,1.0d0,dbl_mb(k_frc_2el),1,dbl_mb(k_force),1)
        else
          call daxpy(lforce,1.0d0,dbl_mb(k_frc_2el_j),1,
     1               dbl_mb(k_force), 1)
          call daxpy(lforce,1.0d0,dbl_mb(k_frc_2el_k),1,
     1               dbl_mb(k_force), 1)
        endif
        call daxpy(lforce,1.0d0,dbl_mb(k_frc_df0),1,dbl_mb(k_force),1)
        call daxpy(lforce,1.0d0,dbl_mb(k_frc_df1),1,dbl_mb(k_force),1)
        call daxpy(lforce,1.0d0,dbl_mb(k_frc_dft),1,dbl_mb(k_force),1)
        ilo(1) = ir
        ihi(1) = ir
        ilo(2) = 1
        ihi(2) = 3
        ilo(3) = 1
        ihi(3) = nat
        ld(1)  = 1
        ld(2)  = 3
        ld(3)  = nat
        call nga_acc(g_g,ilo,ihi,dbl_mb(k_force),ld,1.0d0/nproc)
c
c       Print results if so required
c
        if (util_print('tddft grad terms',print_high)) then
          if (oroot) then
            write(LuOut,5)ir
    5       format(' The Excited State Energy Gradient by Terms',
     +             ' for State ',i3)
            write(LuOut,*)'- Nuclear Energy Gradient:'
            call print_matrix(dbl_mb(k_frc_nuc),3,nat,3)
            write(LuOut,*)'- One-Electron Gradient:'
            call print_matrix(dbl_mb(k_frc_kin),3,nat,3)
            write(LuOut,*)'- Overlap Gradient:'
            call print_matrix(dbl_mb(k_frc_wgh),3,nat,3)
c Daniel (4-8-13): For CAM functionals
            if (.not.cam_exch) then
              write(LuOut,*)'- Two-Electron Gradient:'
              call print_matrix(dbl_mb(k_frc_2el),3,nat,3)
            else
              write(LuOut,*)'- Coulomb Gradient (CAM):'
              call print_matrix(dbl_mb(k_frc_2el_j),3,nat,3)
              write(LuOut,*)'- Exchange Gradient (CAM):'
              call print_matrix(dbl_mb(k_frc_2el_k),3,nat,3)
            endif
            write(LuOut,*)'- Ground State DFT Gradient:'
            call print_matrix(dbl_mb(k_frc_df0),3,nat,3)
c            write(LuOut,*)'- Transition Density DFT Gradient:'
            write(LuOut,*)'- XC-potential TDDFT Gradient:'
            call print_matrix(dbl_mb(k_frc_df1),3,nat,3)
c            write(LuOut,*)'- Third-Order Functional Derivative Term:'
            write(LuOut,*)'- XC-kernel TDDFT Gradient:'
            call print_matrix(dbl_mb(k_frc_dft),3,nat,3)
            write(LuOut,*)'* Total Excited State Gradient:'
            call print_matrix(dbl_mb(k_force),3,nat,3)
          endif
        endif
      enddo

 9200 format(5x,'Dipole Moment    X',f9.5,'   Y',f9.5,'   Z',f9.5)
c
c     Tidy up
c
      call int_terminate()
      call schwarz_tidy()
c
c     Release the memory
c
      if (.not. ma_pop_stack(l_pdm2d))
     +     call errquit(pname//'could not deallocate pdm2d',nao,MA_ERR)
      if (.not. ma_pop_stack(l_bftoat))
     +     call errquit(pname//'could not deallocate bftoat',nao,MA_ERR)
      if (.not. ma_pop_stack(l_rbfmap))
     +     call errquit(pname//'could not deallocate rbfmap',nao,MA_ERR)
      if (.not. ma_pop_stack(l_bfmap))
     +     call errquit(pname//'could not deallocate bfmap',nao,MA_ERR)
      if (.not. ma_pop_stack(l_bfghi))
     +     call errquit(pname//'could not deallocate bhi',nao,MA_ERR)
      if (.not. ma_pop_stack(l_bfglo))
     +     call errquit(pname//'could not deallocate blo',nao,MA_ERR)
      if (.not. ma_pop_stack(l_shbfhi))
     +     call errquit(pname//'could not deallocate bfhi',nsh,MA_ERR)
      if (.not. ma_pop_stack(l_shbflo))
     +     call errquit(pname//'could not deallocate bflo',nsh,MA_ERR)
      if (.not. ma_pop_stack(l_shghi))
     +     call errquit(pname//'could not deallocate bhi',nsh,MA_ERR)
      if (.not. ma_pop_stack(l_shglo))
     +     call errquit(pname//'could not deallocate blo',nsh,MA_ERR)
      if (.not. ma_pop_stack(l_shmap))
     +     call errquit(pname//'could not deallocate shmap',nsh,MA_ERR)
c
      if (.not.ma_pop_stack(l_act))
     +  call errquit(pname//'could not deallocate act',nat,MA_ERR)
      if (.not.ma_pop_stack(l_q4))
     +  call errquit(pname//'could not deallocate q4', maxq,MA_ERR)
      if (.not.ma_pop_stack(l_list))
     +  call errquit(pname//'could not deallocate list',4*maxq,MA_ERR)
      if (.not.ma_pop_stack(l_labels))
     +  call errquit(pname//'could not deallocate labels',lbuf/3,MA_ERR)
c
      if (ipol.gt.1) then
c Daniel (1-7-13): We can avoid this part if pure DFT is used.
c Otherwise, we need the (X-Y) contribution, even if CIS/TDA is
c being performed.
        if (lhashf) then
          if (.not.ma_pop_stack(lxmy_jlb))
     +      call errquit(pname//'could not deallocate xmy_jlb',0,MA_ERR)
          if (.not.ma_pop_stack(lxmy_ikb))
     +      call errquit(pname//'could not deallocate xmy_ikb',0,MA_ERR)
          if (.not.ma_pop_stack(lxmy_jkb))
     +      call errquit(pname//'could not deallocate xmy_jkb',0,MA_ERR)
          if (.not.ma_pop_stack(lxmy_ilb))
     +      call errquit(pname//'could not deallocate xmy_ilb',0,MA_ERR)
        endif
c
        if (.not.ma_pop_stack(lxpy_jlb))
     +    call errquit(pname//'could not deallocate xpy_jlb',0,MA_ERR)
        if (.not.ma_pop_stack(lxpy_ikb))
     +    call errquit(pname//'could not deallocate xpy_ikb',0,MA_ERR)
        if (.not.ma_pop_stack(lxpy_jkb))
     +    call errquit(pname//'could not deallocate xpy_jkb',0,MA_ERR)
        if (.not.ma_pop_stack(lxpy_ilb))
     +    call errquit(pname//'could not deallocate xpy_ilb',0,MA_ERR)
        if (.not.ma_pop_stack(lxpy_klb))
     +    call errquit(pname//'could not deallocate xpy_klb',0,MA_ERR)
        if (.not.ma_pop_stack(lxpy_ijb))
     +    call errquit(pname//'could not deallocate xpy_ijb',0,MA_ERR)
c
        if (.not.ma_pop_stack(lp_jkb)) 
     +    call errquit(pname//'could not deallocate p_jkb',0,MA_ERR)
        if (.not.ma_pop_stack(lp_jlb)) 
     +    call errquit(pname//'could not deallocate p_jlb',0,MA_ERR)
        if (.not.ma_pop_stack(lp_klb)) 
     +    call errquit(pname//'could not deallocate p_klb',0,MA_ERR)
        if (.not.ma_pop_stack(lp_ikb)) 
     +    call errquit(pname//'could not deallocate p_ikb',0,MA_ERR)
        if (.not.ma_pop_stack(lp_ilb)) 
     +    call errquit(pname//'could not deallocate p_ilb',0,MA_ERR)
        if (.not.ma_pop_stack(lp_ijb)) 
     +    call errquit(pname//'could not deallocate p_ijb',0,MA_ERR)
c
        if (.not.ma_pop_stack(ld_jkb)) 
     +    call errquit(pname//'could not deallocate d_jkb',0,MA_ERR)
        if (.not.ma_pop_stack(ld_jlb)) 
     +    call errquit(pname//'could not deallocate d_jlb',0,MA_ERR)
        if (.not.ma_pop_stack(ld_klb)) 
     +    call errquit(pname//'could not deallocate d_klb',0,MA_ERR)
        if (.not.ma_pop_stack(ld_ikb)) 
     +    call errquit(pname//'could not deallocate d_ikb',0,MA_ERR)
        if (.not.ma_pop_stack(ld_ilb)) 
     +    call errquit(pname//'could not deallocate d_ilb',0,MA_ERR)
        if (.not.ma_pop_stack(ld_ijb)) 
     +    call errquit(pname//'could not deallocate d_ijb',0,MA_ERR)
      endif
c
c Daniel (1-7-13): We can avoid this part if pure DFT is used.
c Otherwise, we need the (X-Y) contribution, even if CIS/TDA is
c being performed.
      if (lhashf) then
        if (.not.ma_pop_stack(lxmy_jla))
     +    call errquit(pname//'could not deallocate xmy_jla',0,MA_ERR)
        if (.not.ma_pop_stack(lxmy_ika))
     +    call errquit(pname//'could not deallocate xmy_ika',0,MA_ERR)
        if (.not.ma_pop_stack(lxmy_jka))
     +    call errquit(pname//'could not deallocate xmy_jka',0,MA_ERR)
        if (.not.ma_pop_stack(lxmy_ila))
     +    call errquit(pname//'could not deallocate xmy_ila',0,MA_ERR)
      endif
c
      if (.not.ma_pop_stack(lxpy_jla))
     +  call errquit(pname//'could not deallocate xpy_jla',0,MA_ERR)
      if (.not.ma_pop_stack(lxpy_ika))
     +  call errquit(pname//'could not deallocate xpy_ika',0,MA_ERR)
      if (.not.ma_pop_stack(lxpy_jka))
     +  call errquit(pname//'could not deallocate xpy_jka',0,MA_ERR)
      if (.not.ma_pop_stack(lxpy_ila))
     +  call errquit(pname//'could not deallocate xpy_ila',0,MA_ERR)
      if (.not.ma_pop_stack(lxpy_kla))
     +  call errquit(pname//'could not deallocate xpy_kla',0,MA_ERR)
      if (.not.ma_pop_stack(lxpy_ija))
     +  call errquit(pname//'could not deallocate xpy_ija',0,MA_ERR)
c
      if (.not.ma_pop_stack(lp_jka)) 
     +  call errquit(pname//'could not deallocate p_jka',0,MA_ERR)
      if (.not.ma_pop_stack(lp_jla)) 
     +  call errquit(pname//'could not deallocate p_jla',0,MA_ERR)
      if (.not.ma_pop_stack(lp_kla)) 
     +  call errquit(pname//'could not deallocate p_kla',0,MA_ERR)
      if (.not.ma_pop_stack(lp_ika)) 
     +  call errquit(pname//'could not deallocate p_ika',0,MA_ERR)
      if (.not.ma_pop_stack(lp_ila)) 
     +  call errquit(pname//'could not deallocate p_ila',0,MA_ERR)
      if (.not.ma_pop_stack(lp_ija)) 
     +  call errquit(pname//'could not deallocate p_ija',0,MA_ERR)
c
      if (.not.ma_pop_stack(ld_jka)) 
     +  call errquit(pname//'could not deallocate d_jka',0,MA_ERR)
      if (.not.ma_pop_stack(ld_jla)) 
     +  call errquit(pname//'could not deallocate d_jla',0,MA_ERR)
      if (.not.ma_pop_stack(ld_kla)) 
     +  call errquit(pname//'could not deallocate d_kla',0,MA_ERR)
      if (.not.ma_pop_stack(ld_ika)) 
     +  call errquit(pname//'could not deallocate d_ika',0,MA_ERR)
      if (.not.ma_pop_stack(ld_ila)) 
     +  call errquit(pname//'could not deallocate d_ila',0,MA_ERR)
      if (.not.ma_pop_stack(ld_ija)) 
     +  call errquit(pname//'could not deallocate d_ija',0,MA_ERR)
c
c     Other stuff
c
      if (.not.ma_pop_stack(l_wdens)) 
     +  call errquit(pname//'could not deallocate l_wdens',0,MA_ERR)
      if (.not.ma_pop_stack(l_dens)) 
     +  call errquit(pname//'could not deallocate l_dens',0,MA_ERR)
      if (.not.ma_pop_stack(l_scr)) 
     +  call errquit(pname//'could not deallocate l_scr',0,MA_ERR)
      if (.not.ma_pop_stack(l_buf)) 
     +  call errquit(pname//'could not deallocate l_buf',0,MA_ERR)
      if (.not.ma_pop_stack(l_frc_dft)) 
     +  call errquit(pname//"could not deallocate l_frc_dft",0,MA_ERR)
      if (.not.ma_pop_stack(l_frc_df1)) 
     +  call errquit(pname//"could not deallocate l_frc_df1",0,MA_ERR)
      if (.not.ma_pop_stack(l_frc_df0)) 
     +  call errquit(pname//"could not deallocate l_frc_df0",0,MA_ERR)
c Daniel (4-8-13): For CAM functionals
      if (cam_exch) then
        if (.not.ma_pop_stack(l_frc_2el_k)) 
     +    call errquit(pname//"could not deallocate l_frc_2el_k",
     +      0,MA_ERR)
        if (.not.ma_pop_stack(l_frc_2el_j)) 
     +    call errquit(pname//"could not deallocate l_frc_2el_j",
     +      0,MA_ERR)
      else
        if (.not.ma_pop_stack(l_frc_2el)) 
     +    call errquit(pname//"could not deallocate l_frc_2el",0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_frc_wgh)) 
     +  call errquit(pname//"could not deallocate l_frc_wgh",0,MA_ERR)
      if (.not.ma_pop_stack(l_frc_kin)) 
     +  call errquit(pname//"could not deallocate l_frc_kin",0,MA_ERR)
      if (.not.ma_pop_stack(l_frc_nuc)) 
     +  call errquit(pname//"could not deallocate l_frc_nuc",0,MA_ERR)
      if (.not.ma_pop_stack(l_force)) 
     +  call errquit(pname//"could not deallocate l_force",0,MA_ERR)
c
c     Destroy scratch global arrays
c
      if (.not.ga_destroy(g_rhs_xc(1))) 
     +  call errquit(pname//"could not destroy g_rhs_xc",0,GA_ERR)
c
      if (xc_gotxc()) then
        do ip = 1, 2*ipol
          if (.not.ga_destroy(g_dtmp(ip)))
     1      call errquit(pname//"could not destroy g_dtmp_ao", 0,
     2        GA_ERR) 
        enddo
      endif
c
      do ip = 1, ipol
        if (.not.ga_destroy(g_d(0+ip)))
     +   call errquit(pname//"could not destroy g_d_ao",0,GA_ERR)
        if (.not.ga_destroy(g_d(2+ip)))
     +   call errquit(pname//"could not destroy g_t_ao",0,GA_ERR)
        if (.not.ga_destroy(g_d(4+ip)))
     +   call errquit(pname//"could not destroy g_xpy_ao",0,GA_ERR)
c
c    Note: The (X-Y) vector has no contribution to the overall
c    gradients if a pure DFT calculation is done.  We don't waste space
c    allocating memory for it.
        if (lhashf) then
          if (.not.ga_destroy(g_d(6+ip)))
     +     call errquit(pname//"could not destroy g_xmy_ao",0,GA_ERR)
        endif
      enddo
c
      if (.not.ga_destroy(g_wp_ao)) 
     + call errquit(pname//"could not destroy g_wp_ao",0,GA_ERR)
      if (.not.ga_destroy(g_p_ao)) 
     + call errquit(pname//"could not destroy g_p_ao",0,GA_ERR)
      if (.not.ga_destroy(g_ovlp)) 
     + call errquit(pname//"could not destroy g_ovlp",0,GA_ERR)
c
      end
c $Id$
