      logical function tce_energy(rtdb,excitedfragment)
!
! $Id$
! 
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "sf.fh"
#include "inp.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_hetio.fh"
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
!#include "paw_params.fh"
!
!     CI, CC, & MBPT
!
      integer rtdb             ! Run-time database
      logical nodezero         ! True if node 0
      logical recompf          ! True if recompute Fock
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      integer g_ao1e(2)        ! GA handle for AO Fock matrices
      integer d_ao2e           ! SF handle for AO 2e integrals
      integer d_f1,d_f1temp    ! SF handle for MO 1e integrals
      integer k_f1temp,l_f1temp
      integer l_f1_offset      ! Offset for 1e integral file
      integer k_f1_offset      ! Offset for 1e integral file
      integer size_1e          ! File size in doubles
      integer d_v2             ! SF handle for MO 2e integrals
      integer l_v2_offset      ! Offset for 2e integral file
      integer k_v2_offset      ! Offset for 2e integral file
      integer size_2e          ! File size in doubles
      integer size_2e_ga,size_2e_sf
! --- relevant if intorb = .true. ----
      integer lax_v2_alpha_offset     ! Offset for 2e integrals if intorb=.true.
      integer kax_v2_alpha_offset     ! Offset for 2e integrals if intorb=.true.
      integer lax_v2_alpha_offset_ga  ! Offset for 2e integrals if intorb=.true.
      integer kax_v2_alpha_offset_ga  ! Offset for 2e integrals if intorb=.true.
      integer lax_v2_alpha_offset_sf  ! Offset for 2e integrals if intorb=.true.
      integer kax_v2_alpha_offset_sf  ! Offset for 2e integrals if intorb=.true.
!      integer l_o2o_offset       ! Offset for offset to offset table
!      integer k_o2o_offset       ! Offset for offset to offset table
!      integer l_o2ot_offset      ! orbitals for offset to offset tables
!      integer k_o2ot_offset      ! orbitals for offset to offset tables
      integer lx_o2o_offset,kx_o2o_offset
      integer lx_o2ot_offset,kx_o2ot_offset
      integer lx_v2spin_offset,kx_v2spin_offset
      integer length
!      integer d_v2orb
! ------------------------------------
      integer d_t1             ! SF handle for t1 amplitudes
      integer l_t1_offset      ! Offset for t1 file
      integer k_t1_offset      ! Offset for t1 file
      integer size_t1          ! File size in doubles
      integer d_t2             ! SF handle for t2 amplitudes
      integer l_t2_offset      ! Offset for t2 file
      integer k_t2_offset      ! Offset for t2 file
      integer size_t2          ! File size in doubles
      integer d_t3             ! SF handle for t3 amplitudes
      integer l_t3_offset      ! Offset for t3 file
      integer k_t3_offset      ! Offset for t3 file
      integer size_t3          ! File size in doubles
      integer d_t4             ! SF handle for t4 amplitudes
      integer l_t4_offset      ! Offset for t4 file
      integer k_t4_offset      ! Offset for t4 file
      integer size_t4          ! File size in doubles
      integer d_r1             ! SF handle for r1 amplitudes
      integer l_r1_offset      ! Offset for r1 file
      integer k_r1_offset      ! Offset for r1 file
      integer size_r1          ! File size in doubles
      integer d_r2             ! SF handle for r2 amplitudes
      integer l_r2_offset      ! Offset for r2 file
      integer k_r2_offset      ! Offset for r2 file
      integer size_r2          ! File size in doubles
      integer d_r3             ! SF handle for r3 amplitudes
      integer l_r3_offset      ! Offset for r3 file
      integer k_r3_offset      ! Offset for r3 file
      integer size_r3          ! File size in doubles
      integer d_r4             ! SF handle for r4 amplitudes
      integer l_r4_offset      ! Offset for r4 file
      integer k_r4_offset      ! Offset for r4 file
      integer size_r4          ! File size in doubles
      integer d_e              ! SF handle for e file
      integer l_e_offset       ! Offset for e file
      integer k_e_offset       ! Offset for e file
      integer size_e           ! File size in doubles
      integer d_u1             ! SF handle for u1 file
      integer l_u1_offset      ! Offset for u1 file
      integer k_u1_offset      ! Offset for u1 file
      integer size_u1          ! File size in doubles
      integer d_u2             ! SF handle for u2 file
      integer l_u2_offset      ! Offset for u2 file
      integer k_u2_offset      ! Offset for u2 file
      integer size_u2          ! File size in doubles
      integer d_w1             ! SF handle for w1 file
      integer l_w1_offset      ! Offset for w1 file
      integer k_w1_offset      ! Offset for w1 file
      integer size_w1          ! File size in doubles
      integer d_w2             ! SF handle for w2 file
      integer l_w2_offset      ! Offset for w2 file
      integer k_w2_offset      ! Offset for w2 file
      integer size_w2          ! File size in doubles
      integer d_w3             ! SF handle for w3 file
      integer l_w3_offset      ! Offset for w3 file
      integer k_w3_offset      ! Offset for w3 file
      integer size_w3          ! File size in doubles
      integer d_w4             ! SF handle for w4 file
      integer l_w4_offset      ! Offset for w4 file
      integer k_w4_offset      ! Offset for w4 file
      integer size_w4          ! File size in doubles
      integer d_e2             ! SF handle for e2 file
      integer l_e2_offset      ! Offset for e2 file
      integer k_e2_offset      ! Offset for e2 file
      integer size_e2          ! File size in doubles
      double precision ref     ! Ground state energy
      double precision corr    ! Correlation energy
      double precision mbpt2   ! Correlation energy
      double precision mbpt3   ! Correlation energy
      double precision mbpt4   ! Correlation energy
      double precision r1      ! Residual t1
      double precision r2      ! Residual t2
      double precision r3      ! Residual t3
      double precision r4      ! Residual t4
      double precision residual! Largest residual
      integer dummy            ! Dummy argument for DIIS
!!!     Sigma2 
!!!      logical sig2
!!!      integer l_s2_re,k_s2_re
!!!      integer l_s2_im,k_s2_im
!!!      integer l_fs2_a,k_fs2_a
!!!      integer size_s2
!!!      double precision eta
! FOR LRH-SD
      integer d_j1,d_j2
      integer l_j1_offset,l_j2_offset
      integer k_j1_offset,k_j2_offset
      integer size_j1,size_j2
!QuantumComputing
      integer g3b,g4b,g1b,g2b,g5b,g6b
      integer g5b_orb,g6b_orb
      integer igi3,igi4,igi1,igi2,igi5,igi6
      integer tg1,tg2,tg3,tg4
      integer tirow,ticol
      integer offset_test,size_test
      integer irow,icol,index_pair
      integer k_aux,l_aux
      integer range_ov,range_occ
      integer nos,nvs,noas,nobs,nvas,nvbs,nactv,mov
      integer nacto,nrotv,nrot
      integer q_ele_a,q_ele_b,q_orb
      double precision xxx,xnorm1,xnorm2,xnorm
      double precision l1norm,ehfa,enrep_tce
c      double precision v4m(50,50,50,50)
c      double precision f2m(50,50)
      integer k_f2m,l_f2m,size_f2m
      integer k_v4m,l_v4m,size_v4m
      integer l_xrot,k_xrot
      logical quantumc,spust
      logical qrotation   !quantum rotation in active orbital space
c
      logical pwtce
c
      logical qducc
c ucc4 with intorb option
      integer l_v2aux_offset,k_v2aux_offset,size_2eaux
      integer d_v2copy,k_v2copy_offset
      integer d_v2aux
c
      logical ucc4
      character*4 hftype
!QuantumComputing
! FOR CR-EOMCCSD(T) ---------------------
      integer d_ex1,d_ex2
      integer l_ex1_offset,l_ex2_offset
      integer k_ex1_offset,k_ex2_offset
      integer size_ex1,size_ex2
      integer d_c1,d_c2
      integer l_c1_offset,l_c2_offset
      integer k_c1_offset,k_c2_offset
      integer size_c1,size_c2
      double precision excit
      double precision energy_grcr
!  NEW SOLVER ---------------------------
      integer l_hbar,k_hbar,ihbar
      integer imicro,xdiis,iroot
      integer l_c_eff,k_c_eff
      integer l_h_eff,k_h_eff
      integer ipos,ipos1,ipos2
! - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer l_x1_local,k_x1_local
! ---------------------------------------
! - intermediates reuse --
      integer l_aux_off,k_aux_off,size_aux
      integer ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5
      integer rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,rec_d_6
      integer rec_d_7,rec_d_8,rec_d_9
! - new restarts and eom iter.control ---------
      integer unitn,recforio
      integer xmem
      logical res_log
      logical active_log
      logical eom_start
      logical slice
      logical d4d5
      double precision r0xx
      double precision eom_excit(4),eom_r0(4)
      double precision thresh_l,aux_x
      double precision thresh_eom
! ---------------------------------------------
      logical needt1
      logical needt2
      logical needt3
      logical needt3a
      logical needt4
      logical needu1
      logical needu2
      logical needw1
      logical needw2
      logical needw3
      logical needw4
      logical neede2
      logical do_mbpt
      character*255 filename
      character*4 irrepname
      logical dft_energy
      external dft_energy
      logical scf
      external scf
      integer i
      integer j
#if defined(TCE_CUDA) || defined(TCE_HIP)
      integer icuda
#endif
!#ifdef EACCSD
      logical eaccsd
!#endif
!#ifdef IPCCSD
      logical ipccsd
!#endif
      logical need_mutex

! --- debug ---
!      integer dup1,dup2,dup3,dup4,dup5,dup6
! -------------
!
!     CC Lambda equation
!
#if !defined(NEWDIPOLE)
      integer g_aod1(3)         ! GA handle for AO dipole matrices
      integer g_aoq1(6)         ! GA handle for AO quadrupole matrices
      integer g_aoo1(10)        ! GA handle for AO octupole matrices
#endif
      integer d_d1(3)           ! MO dipole integral file
      integer l_d1_offset(3)    ! Offset for d1 integral file
      integer k_d1_offset(3)    ! Offset for d1 integral file
      integer size_d1(3)        ! File size in doubles
      character*3 dfilename(3)  ! File name stub
      data dfilename/'d1x','d1y','d1z'/
      integer d_q1(6)           ! MO quadrupole integral file
      integer l_q1_offset(6)    ! Offset for q1 integral file
      integer k_q1_offset(6)    ! Offset for q1 integral file
      integer size_q1(6)        ! File size in doubles
      character*4 qfilename(6)  ! File name stub
      data qfilename/'q1xx','q1yy','q1zz','q1xy','q1xz','q1yz'/
      integer d_o1(10)          ! MO octupole integral file
      integer l_o1_offset(10)   ! Offset for o1 integral file
      integer k_o1_offset(10)   ! Offset for o1 integral file
      integer size_o1(10)       ! File size in doubles
      character*5 ofilename(10) ! File name stub
      data ofilename/'o1xxx','o1yyy','o1zzz','o1xxy','o1xxz',
     &               'o1yyx','o1yyz','o1zzx','o1zzy','o1xyz'/
      character*3 axisname(19)  ! Axis
      data axisname/'X','Y','Z','XX','YY','ZZ','XY','XZ','YZ',
     & 'XXX','YYY','ZZZ','XXY','XXZ','YYX','YYZ','ZZX','ZZY','XYZ'/
      double precision debye,buck,octo
      parameter (debye=2.54176568d0)
      parameter (buck=1.344911d0)
      parameter (octo=0.711688d0)
      integer d_lambda1        ! Lambda 1 file
      integer d_lambda2        ! Lambda 2 file
      integer d_lambda3        ! Lambda 3 file
      integer d_lambda4        ! Lambda 4 file
      integer d_lr1            ! Lambda 1 residual file
      integer d_lr2            ! Lambda 2 residual file
      integer d_lr3            ! Lambda 3 residual file
      integer d_lr4            ! Lambda 4 residual file
      integer l_l1_offset      ! Offset for Lambda 1 file
      integer k_l1_offset      ! Offset for Lambda 1 file
      integer size_l1          ! File size in doubles
      integer l_l2_offset      ! Offset for Lambda 2 file
      integer k_l2_offset      ! Offset for Lambda 2 file
      integer size_l2          ! File size in doubles
      integer l_l3_offset      ! Offset for Lambda 3 file
      integer k_l3_offset      ! Offset for Lambda 3 file
      integer size_l3          ! File size in doubles
      integer l_l4_offset      ! Offset for Lambda 4 file
      integer k_l4_offset      ! Offset for Lambda 4 file
      integer size_l4          ! File size in doubles
      logical idens            ! true if calculates density matrix
      integer d_d0             ! Dipole moment
      integer l_d0_offset      ! Offset for d0 file
      integer k_d0_offset      ! Offset for d0 file
      integer size_d0          ! File size in doubles
      integer d_q0             ! Quadrupole moment
      integer l_q0_offset      ! Offset for q0 file
      integer k_q0_offset      ! Offset for q0 file
      integer size_q0          ! File size in doubles
      integer d_o0             ! Octupole moment
      integer l_o0_offset      ! Offset for o0 file
      integer k_o0_offset      ! Offset for o0 file
      integer size_o0          ! File size in doubles
      integer axis             ! X=1, Y=2, & Z=3
      double precision dipole_scf_t      ! Dipole moment - SCF total
      double precision dipole_cor_t      ! Dipole moment - CC total
      double precision dipole_exc_t      ! Dipole moment - EOM-CC total
      double precision dipole_scf(3)     ! Dipole moment - SCF XYZ
      double precision dipole_cor(3)     ! Dipole moment - CC XYZ
      double precision dipole_exc(3)     ! Dipole moment - EOM-CC
      double precision cdipole           ! Dipole moment correlation
      double precision quadrupole_scf(6) ! Quadrupole moment
      double precision quadrupole_cor(6) ! Quadrupole moment
      double precision quadrupole_exc(6) ! Quadrupole moment
      double precision cquadrupole       ! Quadrupole moment correlation
      double precision octupole_scf(10)  ! Octupole moment
      double precision octupole_cor(10)  ! Octupole moment
      double precision octupole_exc(10)  ! Octupole moment
      double precision coctupole         ! Octupole moment correlation
      double precision ddotfile
      external ddotfile
      integer sym_abelian_axis ! Move this and the following to "sym.fh"
      external sym_abelian_axis
!
!     CC Response
!
      character*5 rr1filename(3) ! File name stub
      data rr1filename/'rr1x ','rr1y ','rr1z '/
      character*5 rr2filename(3) ! File name stub
      data rr2filename/'rr2x ','rr2y ','rr2z '/
      character*5 rr3filename(3) ! File name stub
      data rr3filename/'rr3x ','rr3y ','rr3z '/
      character*5 rr4filename(3) ! File name stub
      data rr4filename/'rr4x ','rr4y ','rr4z '/
      character*5 tr1filename(9) ! File name stub
      data tr1filename/'tr1x+','tr1y+','tr1z+',
     1                 'tr1x-','tr1y-','tr1z-',
     2                 'tr1xi','tr1yi','tr1zi'/
      character*5 tr2filename(9) ! File name stub
      data tr2filename/'tr2x+','tr2y+','tr2z+',
     1                 'tr2x-','tr2y-','tr2z-',
     2                 'tr2xi','tr2yi','tr2zi'/
      character*5 tr3filename(9) ! File name stub
      data tr3filename/'tr3x+','tr3y+','tr3z+',
     1                 'tr3x-','tr3y-','tr3z-',
     2                 'tr3xi','tr3yi','tr3zi'/
      character*5 tr4filename(9) ! File name stub
      data tr4filename/'tr4x ','tr4y ','tr4z ',
     1                 'tr4x-','tr4y-','tr4z-',
     2                 'tr4xi','tr4yi','tr4zi'/
!
      character*5 yrr1filename(3) ! File name stub
      data yrr1filename/'yrr1x','yrr1y','yrr1z'/
      character*5 yrr2filename(3) ! File name stub
      data yrr2filename/'yrr2x','yrr2y','yrr2z'/
      character*5 yrr3filename(3) ! File name stub
      data yrr3filename/'yrr3x','yrr3y','yrr3z'/
      character*5 yrr4filename(3) ! File name stub
      data yrr4filename/'yrr4x','yrr4y','yrr4z'/
      character*5 yr1filename(6) ! File name stub
      data yr1filename/'yr1x+','yr1y+','yr1z+','yr1x-','yr1y-','yr1z-'/
      character*5 yr2filename(6) ! File name stub
      data yr2filename/'yr2x+','yr2y+','yr2z+','yr2x-','yr2y-','yr2z-'/
      character*5 yr3filename(6) ! File name stub
      data yr3filename/'yr3x+','yr3y+','yr3z+','yr3x-','yr3y-','yr3z-'/
      character*5 yr4filename(6) ! File name stub
      data yr4filename/'yr4x+','yr4y+','yr4z+','yr4x-','yr4y-','yr4z-'/
!     RESPONSE AMPLITUDES
      integer d_tr1(6)             ! SF handle for tr1 amplitudes
      integer l_tr1_offset(3)      ! Offset for tr1 file
      integer k_tr1_offset(3)      ! Offset for tr1 file
      integer size_tr1(3)          ! File size in doubles
      integer d_tr2(6)             ! SF handle for tr2 amplitudes
      integer l_tr2_offset(3)      ! Offset for tr2 file
      integer k_tr2_offset(3)      ! Offset for tr2 file
      integer size_tr2(3)          ! File size in doubles
      integer d_tr3(6)             ! SF handle for tr3 amplitudes
      integer l_tr3_offset(3)      ! Offset for tr3 file
      integer k_tr3_offset(3)      ! Offset for tr3 file
      integer size_tr3(3)          ! File size in doubles
      integer d_tr4(6)             ! SF handle for tr4 amplitudes
      integer l_tr4_offset(3)      ! Offset for tr4 file
      integer k_tr4_offset(3)      ! Offset for tr4 file
      integer size_tr4(3)          ! File size in doubles
!     RESPONSE RESIDUALS
      integer d_rr1(3)             ! SF handle for rr1 amplitudes
      integer d_rr2(3)             ! SF handle for rr2 amplitudes
      integer d_rr3(3)             ! SF handle for rr3 amplitudes
      integer d_rr4(3)             ! SF handle for rr4 amplitudes
!     LAMBDA RESPONSE AMPLITUDES
      integer d_yr1(6)             ! SF handle for yr1 amplitudes
      integer l_yr1_offset(3)      ! Offset for yr1 file
      integer k_yr1_offset(3)      ! Offset for yr1 file
      integer size_yr1(3)          ! File size in doubles
      integer d_yr2(6)             ! SF handle for yr2 amplitudes
      integer l_yr2_offset(3)      ! Offset for yr2 file
      integer k_yr2_offset(3)      ! Offset for yr2 file
      integer size_yr2(3)          ! File size in doubles
      integer d_yr3(6)             ! SF handle for yr3 amplitudes
      integer l_yr3_offset(3)      ! Offset for yr3 file
      integer k_yr3_offset(3)      ! Offset for yr3 file
      integer size_yr3(3)          ! File size in doubles
      integer d_yr4(6)             ! SF handle for yr4 amplitudes
      integer l_yr4_offset(3)      ! Offset for yr4 file
      integer k_yr4_offset(3)      ! Offset for yr4 file
      integer size_yr4(3)          ! File size in doubles
!     LAMBDA RESPONSE RESIDUALS
      integer d_yrr1(3)             ! SF handle for yrr1 amplitudes
      integer d_yrr2(3)             ! SF handle for yrr2 amplitudes
      integer d_yrr3(3)             ! SF handle for yrr3 amplitudes
      integer d_yrr4(3)             ! SF handle for yrr4 amplitudes
!     OTHER RESPONSE STUFF
      double precision omega       ! The frequency used in (A-Iw)x=b
      double precision lambda      ! The damping of Ab in (A^2-Iw^2)+lambda*Ab=0
      double precision omega1     ! The frequency used in (A-Iw)x=b
      double precision omega2     ! The frequency used in (A-Iw)x=b
      double precision omega3     ! The frequency used in (A-Iw)x=b
      integer omegacount
      integer omegasign
      integer dynfreq
      integer dynaxis
      integer axisA
      integer axisB
      integer axisC
      integer axisD
      integer d_a0                ! Dipole polarizability component handle
      integer l_a0_offset         ! Offset for a0 file
      integer k_a0_offset         ! Offset for a0 file
      integer size_a0             ! File size in doubles
      double precision alpha(3,3) ! Dipole polarizability tensor
      double precision alpha1     ! Dipole polarizability term 1
      double precision alpha2     ! Dipole polarizability term 2
      double precision alpha3     ! Dipole polarizability term 3
      double precision alpha4     ! Dipole polarizability term 4
      integer d_b0                ! Dipole hyperpolarizability component handle
      integer l_b0_offset         ! Offset for b0 file
      integer k_b0_offset         ! Offset for b0 file
      integer size_b0             ! File size in doubles
      double precision beta(3,3,3)! Dipole hyperpolarizability tensor
      double precision beta1      ! Dipole hyperpolarizability term 1
      double precision beta2      ! Dipole hyperpolarizability term 2
      double precision beta3      ! Dipole hyperpolarizability term 3
      double precision beta4      ! Dipole hyperpolarizability term 4
      integer d_g0                ! Dipole second hyperpolarizability component handle
      integer l_g0_offset         ! Offset for g0 file
      integer k_g0_offset         ! Offset for g0 file
      integer size_g0             ! File size in doubles
!      double precision gamm(3,3,3,3)! Dipole second hyperpolarizability tensor (gamma not safe variable name)
!      double precision gamm1        ! Dipole second hyperpolarizability term 1
!      double precision gamm2        ! Dipole second hyperpolarizability term 2
!      double precision gamm3        ! Dipole second hyperpolarizability term 3
!      double precision gamm4        ! Dipole second hyperpolarizability term 4
      double precision rr1,rr2,rr3,rr4
      double precision t1diagn ! T1 diagnostic
      double precision alphacopy(3,3)   ! Dipole polarizability tensor copy
      double precision alphaiso         ! Isotropic dipole polarizability
      double precision alphaani         ! Anisotropic dipole polarizability
      double precision alphaevr(3)      ! Dipole polarizability tensor eigenvalues (real)
      double precision alphaevi(3)      ! Dipole polarizability tensor eigenvalues (imag)
      double precision aeigvec(3,3)     ! dummy dgeev eigenvectors
!      integer iter             ! Iteration counter
!
!     Excited-state CI, CC, & MBPT
!
      integer irrep            ! Symmetry loop index
      integer irrep_g          ! Ground state symmetry
      integer d_rx1            ! RHS residual file
      integer d_rx2            ! RHS residual file
      integer d_rx3            ! RHS residual file
      integer d_rx4            ! RHS residual file
      integer d_ry1            ! LHS residual file
      integer d_ry2            ! LHS residual file
      integer d_ry3            ! LHS residual file
      integer d_ry4            ! LHS residual file
      integer d_x0             ! Offset for x0 file
      integer l_x0_offset      ! Offset for x0 file
      integer k_x0_offset      ! Offset for x0 file
      integer size_x0          ! File size in doubles
      integer l_x1_offset      ! Offset for x1 file
      integer k_x1_offset      ! Offset for x1 file
      integer size_x1          ! File size in doubles
      integer l_x2_offset      ! Offset for x2 file
      integer k_x2_offset      ! Offset for x2 file
      integer size_x2          ! File size in doubles
      integer l_x3_offset      ! Offset for x3 file
      integer k_x3_offset      ! Offset for x3 file
      integer size_x3          ! File size in doubles
      integer l_x4_offset      ! Offset for x4 file
      integer k_x4_offset      ! Offset for x4 file
      integer size_x4          ! File size in doubles
      integer d_y0             ! Offset for y0 file
      integer l_y0_offset      ! Offset for y0 file
      integer k_y0_offset      ! Offset for y0 file
      integer size_y0          ! File size in doubles
      integer l_y1_offset      ! Offset for y1 file
      integer k_y1_offset      ! Offset for y1 file
      integer size_y1          ! File size in doubles
      integer l_y2_offset      ! Offset for y2 file
      integer k_y2_offset      ! Offset for y2 file
      integer size_y2          ! File size in doubles
      integer l_y3_offset      ! Offset for y3 file
      integer k_y3_offset      ! Offset for y3 file
      integer size_y3          ! File size in doubles
      integer l_y4_offset      ! Offset for y4 file
      integer k_y4_offset      ! Offset for y4 file
      integer size_y4          ! File size in doubles
      integer ivec,jvec        ! Current trial vector
      integer l_omegax         ! Excitation energy
      integer k_omegax         ! Excitation energy
      integer l_omegay         ! Excitation energy
      integer k_omegay         ! Excitation energy
      integer l_residual       ! Residual
      integer k_residual       ! Residual
      integer max_ioprocs_in_org,max_ioprocs_in_new
      logical converged
      double precision pt3_1   ! [T] energy
      double precision pt3_2   ! (T)-[T] energy
! LR-CCSD(T) method
      double precision pt3_3,pt3_4,pt3_5,pt3_6
! --------------------------------------------
! LR-CCSD(TQ) method
      double precision pt4_lr
! --------------------------------------------
!
      double precision au2ev    ! Conversion factor from a.u. to eV
      double precision au2ang   ! Conversion factor from bohr to Angstrom
      double precision au2ang3  ! Conversion factor from bohr^3 to Angstrom^3
      double precision auXnm    ! Conversion factor from a.u. (frequency) to nm (wavelength)
      double precision numerator   ! <L|D exp(T)|R>
      double precision denominator ! <L|exp(T)|R>
      double precision transition(3) ! Transition moments
      parameter (au2ev=27.2113961d0)
      parameter (auXnm=45.563353d0)
      parameter (au2ang=5.29177249d-01)
      parameter (au2ang3=au2ang*au2ang*au2ang)
      double precision pi
      parameter(pi = 3.14159265358979323846264338327950288419D0)
      character*255 modelname
      integer fldgts
cducc
      integer l_transf,k_transf,
     1      l_act1,k_act1,
     1      l_prt1,k_prt1,
     1      l_act2,k_act2,
     1      l_prt2,k_prt2,
     1      l_h,k_h,
     1      l_horb,k_horb,
     1      l_v,k_v,
     1      l_t1,k_t1,
     1      l_t2,k_t2,
     1      l_ht,k_ht,
     1      l_vt,k_vt
      logical do_ccd, do_lccd, do_ccsd, do_lccsd, do_qcisd
      logical do_pt, do_lam_pt, do_cr_pt, do_lcr_pt, do_2t_pt
!
!     Fragment extension
!
      logical excitedfragment
!
      logical ccsdtleft
      logical tce_skip_scf
      nodezero=ga_nodeid().eq.0
      ccsdtleft=.false.
      do_mbpt=.false.
      filename(1:255) = ''
!
!     =========================
!     Ground-state HF/DFT first
!     =========================
!
      tce_energy = .false.
      if (.not.rtdb_get(rtdb,'tce:reference',mt_int,1,reference)) then
        reference = 1
        if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1    call errquit('tce_energy: failed to write reference',0,
     2    RTDB_ERR)
      endif
!
!     TCE restart precludes SCF calculation
!
      tce_skip_scf = .false.
      if (.not.rtdb_get(rtdb,'tce:read_integrals',mt_log,5,
     1         read_integrals)) then
        read_integrals(1)=.false. ! 1-e integrals
        read_integrals(2)=.false. ! 2-e integrals
        read_integrals(3)=.false. ! dipole integrals
        read_integrals(4)=.false. ! quadrupole integrals
        read_integrals(5)=.false. ! octopole integrals
        if (.not.rtdb_put(rtdb,'tce:read_integrals',mt_log,5,
     1           read_integrals))
     1  call errquit('tce_init: failed to write read_integrals',
     2               0,RTDB_ERR)
        tce_skip_scf = (read_integrals(1).and.read_integrals(2))
      endif
      if (.not.tce_skip_scf) then
        if (reference.eq.0) then
          if (.not.dft_energy(rtdb)) return
          if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,ref))
     1      call errquit('tce_energy: failed to get dft energy',0,
     2      RTDB_ERR)
        elseif (reference.eq.1) then
          if (.not.scf(rtdb)) return
          if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,ref))
     1      call errquit('tce_energy: failed to get scf energy',0,
     2      RTDB_ERR)
        else
          call errquit('tce_energy: illegal reference',reference,
     1      UNKNOWN_ERR)
        endif
      endif
!
!     ===================
!     Print utility start
!     ===================
!
      call util_print_push
      call util_print_rtdb_load(rtdb,'tce')
!
!     ==========
!     Initialize
!     ==========
!
      call tce_init(rtdb)
      if (.not.rtdb_get(rtdb,'tce:densmat',mt_log,1,idens))
     1   call errquit('tce_energy:rtdb_get densmat failed',0, RTDB_ERR)
      if (nodezero) call util_flush(LuOut)
!      if (lineresp) left=.true.
!
!      if (nodezero) then
!        write(LuOut,*) "lineresp = ",lineresp
!        write(LuOut,*) "leftresp = ",leftresp
!        write(LuOut,*) "quadresp = ",quadresp
!        write(LuOut,*) "anumfreq = ",anumfreq
!        do i=1,anumfreq
!          write(LuOut,*) "afreq(",i,") = ",afreq(i)
!        enddo
!        write(LuOut,*) "respaxis(1) = ",respaxis(1)
!        write(LuOut,*) "respaxis(2) = ",respaxis(2)
!        write(LuOut,*) "respaxis(3) = ",respaxis(3)
!        print*,"respaxis ",respaxis(1),respaxis(2),respaxis(3)
!        call util_flush(LuOut)
!      endif
!
!     =================
!     Totally symmetric  RESPONSE WILL CHANGE THESE - BE CAREFUL
!     =================
!
      irrep_x = 0
      irrep_y = 0
      irrep_c = 0
      irrep_g = 0
!
      idiv2e  = 2
!
       thresh_l=thresh
       if (.not.rtdb_get(rtdb,'tce:threshl',mt_dbl,1,aux_x)) then
         if(nodezero) then
          write(LuOut,9520)thresh_l
         end if
       else
          thresh_l=aux_x
          if(nodezero) then
           write(LuOut,9520)thresh_l
          end if
       end if
!
       thresh_eom=thresh
       if (.not.rtdb_get(rtdb,'tce:thresheom',mt_dbl,1,aux_x)) then
         if(nodezero) then
          write(LuOut,9521)thresh_eom
         end if
       else
          thresh_eom=aux_x
          if(nodezero) then
           write(LuOut,9521)thresh_eom
          end if
       end if
!
       eom_start=.false.
       if (.not.rtdb_get(rtdb,'tce:eomstart',mt_log,1,eom_start)) then
         if(nodezero) then
          write(LuOut,*)'no EOMCCSD initial starts read in'
         end if
       else
         if(nodezero) then
          write(LuOut,*)'EOMCCSD initial starts read in)'
         end if
       end if
! restart files
       if (.not.rtdb_get(rtdb,'tce:tceiop',mt_int,1,recforio)) then
cold        rec_mem=1
        rec_mem=2048
       else
        rec_mem=recforio
       end if

       if (.not.rtdb_get(rtdb,'tce:maxioprocs',mt_int,1,max_ioprocs_in)) 
     4      max_ioprocs_in=4
!QuantumComputing
       if (.not.rtdb_get(rtdb,'tce:ucc4',
     4  mt_log,1,ucc4)) ucc4=.false.
!QuantumComputing
       if (.not.rtdb_get(rtdb,'tce:print_integrals',
     4  mt_log,1,quantumc)) quantumc=.false.
!QuantumComputing
       if (.not.rtdb_get(rtdb,'tce:qducc',
     4  mt_log,1,qducc)) qducc=.false.
!QuantumComputing
       if(qducc) then 
       nactv=0
        if (.not.rtdb_get(rtdb,'tce:nactv',mt_int,1,nactv)) then
         call errquit('tce_energy:nactv=0 while qducc=true',0,
     2      RTDB_ERR)
        end if
       end if ! if qducc
!QuantumComputing qrotation
       if (.not.rtdb_get(rtdb,'tce:qrotation',
     4  mt_log,1,qrotation)) qrotation=.false.
       nacto=0
       nrotv=0
       nrot=0
      if(qrotation) then
        if (.not.rtdb_get(rtdb,'tce:nacto',mt_int,1,nacto))
     1    call errquit('tce_energy:nacto=0 while qrotation T',0,
     2      RTDB_ERR)
        if (.not.rtdb_get(rtdb,'tce:nrotv',mt_int,1,nrotv))
     1    call errquit('tce_energy:nrotv=0 while qrotation T',0,
     2      RTDB_ERR)
         nrot=nacto+nrotv
       endif !qrotation --> nacto/nrotv
c qrotation
c
       if(.not.rtdb_cget(rtdb,'scf:scftype',1,hftype))
     1   call errquit('tce_energy: failed to read scftype',0,RTDB_ERR)
c *** debug ***
c       if(nodezero) then
c         write(6,*)'hftype ',hftype
c         call util_flush(6)
c       endif
c *************
c
!QuantumComputing
       if (.not.rtdb_get(rtdb,'tce:pwtce',
     4  mt_log,1,pwtce)) pwtce=.false.
c
!QuantumComputing
!
! Sigma2
!!!       res_log=.false.
!!!       if (.not.rtdb_get(rtdb,'tce:sigma2',
!!!     4  mt_log,1,res_log)) then
!!!        sig2=.false.
!!!       else
!!!        sig2=res_log
!!!       end if
!
       res_log=.false.
       if (.not.rtdb_get(rtdb,'tce:readint',mt_log,1,res_log)) then
        read_int=.false.
       else
        read_int=res_log
       end if
!
       res_log=.false.
       if (.not.rtdb_get(rtdb,'tce:writeint',mt_log,1,res_log)) then
        write_int=.false.
       else
        write_int=res_log
       end if
!
       res_log=.false.
       if (.not.rtdb_get(rtdb,'tce:readt',mt_log,1,res_log)) then
        read_ta=.false.
       else
        read_ta=res_log
       end if
!
       res_log=.false.
       if (.not.rtdb_get(rtdb,'tce:writet',mt_log,1,res_log)) then
        write_ta=.false.
       else
        write_ta=res_log
       end if
!
       res_log=.false.
       if (.not.rtdb_get(rtdb,'tce:readx',mt_log,1,res_log)) then
        read_xa=.false.
       else
        read_xa=res_log
       end if
!
       res_log=.false.
       if (.not.rtdb_get(rtdb,'tce:writex',mt_log,1,res_log)) then
        write_xa=.false.
       else
        write_xa=res_log
       end if
!
       res_log=.false.
       if (.not.rtdb_get(rtdb,'tce:readi3',mt_log,1,res_log)) then
        read_in3=.false.
       else
        read_in3=res_log
       end if
!
       res_log=.false.
       if (.not.rtdb_get(rtdb,'tce:writei3',mt_log,1,res_log)) then
        write_in3=.false.
       else
        write_in3=res_log
       end if
!
       if (.not.rtdb_get(rtdb,'tce:nts',mt_log,1,res_log)) 
     c      res_log=.true.
       if(res_log) ccsd_var='ic'
!
       slice=.false.
       if (rtdb_get(rtdb,'tce:xmem',mt_int,1,xmem))slice=.true. 
!
       d4d5=.false.
       if (.not.rtdb_get(rtdb,'tce:d4d5',mt_log,1,d4d5)) then
         continue
       else
         d4d5=.true.
       end if
!
       if(nodezero) then
        write(LuOut,*)'TCE RESTART OPTIONS'
        write(LuOut,*)'READ_INT:  ',read_int
        write(LuOut,*)'WRITE_INT: ',write_int
        write(LuOut,*)'READ_TA:   ',read_ta
        write(LuOut,*)'WRITE_TA:  ',write_ta
        write(LuOut,*)'READ_XA:   ',read_xa
        write(LuOut,*)'WRITE_XA:  ',write_xa
        write(LuOut,*)'READ_IN3:  ',read_in3
        write(LuOut,*)'WRITE_IN3: ',write_in3
        write(LuOut,*)'SLICE:     ',slice
!cx        write(LuOut,*)'XMEM:      ',xmem
        write(LuOut,*)'D4D5:      ',d4d5
        write(LuOut,*)'quantumc:  ',quantumc
        write(LuOut,*)'qrotation: ',qrotation
        write(LuOut,*)'pwtce:     ',pwtce
        write(LuOut,*)'qducc:     ',qducc
        write(LuOut,*)'nacto:     ',nacto
        write(LuOut,*)'nrotv:     ',nrotv
        write(LuOut,*)'nrot:      ',nrot
        call util_flush(LuOut)
       end if

c     number of digits for parallel file prefix (util_file_name0 last arg)
       if (.not.rtdb_get(rtdb,'tce:filedgts',mt_int,1,fldgts)) 
     A      fldgts=5
       if(max_ioprocs_in.gt.10**fldgts) then
          fldgts=log10(dble(max_ioprocs_in))+1
          if(ga_nodeid().eq.0) write(6,*) 
     R         ' tce:filedgts redefined to ',fldgts
       endif


!
!     ==============
!     Create a mutex
!     ==============
!
      need_mutex=ioalg.ne.2
#if 0
       if(ga_nodeid().eq.0) then
          write(LuOut,*)' value of need_mutex ',need_mutex
          write(LuOut,*)' value of ioalg ',ioalg
       endif
#endif
      if(need_mutex) then
         if (.not.ga_create_mutexes(1))
     1        call errquit('tce_energy: GA problem',0,GA_ERR)
      endif
!
!     =========================
!     Reorder and tile orbitals
!     =========================
!
! quantum rotation
      if(qrotation) then
c  getting from rtdb q_ele_a q_ele_b q_orb
         if (.not.rtdb_get(rtdb,'tce:qela',mt_int,1,q_ele_a))
     1     call errquit('tce_energy: failed to read qela',0,
     2     rtdb_err)
         if (.not.rtdb_get(rtdb,'tce:qelb',mt_int,1,q_ele_b))
     1     call errquit('tce_energy: failed to read qelb',0,
     2     rtdb_err)
         if (.not.rtdb_get(rtdb,'tce:qorb',mt_int,1,q_orb))
     1     call errquit('tce_energy: failed to read qorb',0,
     2     rtdb_err)
c
c       if(nodezero) then
c        write(6,*)'nacto ',nacto
c        write(6,*)'q_ele_a ',q_ele_a
c        write(6,*)'nrotv ',nrotv
c        write(6,*)'nmo(1)-nocc(1)',nmo(1)-nocc(1)
c        write(6,*)'nmo(1)  ',nmo(1)
c        write(6,*)'nocc(1) ',nocc(1)
c        call util_flush(1)
c       endif
c sanity check
      if(nacto.ge.q_ele_a)
     &   call errquit('tce: xort: nacto problem',0,MA_ERR)
      if(nrotv.ge.nmo(1)-nocc(1))
     &   call errquit('tce: xort: nacto problem',0,MA_ERR)
c opening nrot x nrot xrot matrix
       if(.not.ma_push_get(mt_dbl,nrot*nrot,'xrot',l_xrot,k_xrot))
     1     call errquit('tce: xrot MA problem',0,MA_ERR)
c
       do i=1,nrot*nrot
        dbl_mb(k_xrot+i-1)=0.0d0
       enddo
c
      call xrot_ini(nrot,dbl_mb(k_xrot))
c
      call ga_sync()
c
      call rot_vir(rtdb,q_orb,q_ele_a,q_ele_b,nacto,nrotv,
     &                   nrot,dbl_mb(k_xrot))
c
      call ga_sync()
c closing xrot
      if (.not.ma_pop_stack(l_xrot))
     1     call errquit('l_xrot',15,MA_ERR)
      endif !qrotation
c
c
      call tce_tile(rtdb)
      call sf_test
!
!     ================================
!     Parallel integral transformation
!     ================================
!
      if (left) then
        cpu = - util_cpusec()
        wall = - util_wallsec()
        call tce_nud1(dipole_scf)
        if (nodezero.and.util_print('multipole',print_debug)) then
            write(LuOut,*)
            write(LuOut,9120) '==== tce_aod1 ===='
            call util_flush(LuOut)
        endif
        call tce_aod1(rtdb,g_aod1,dipole_scf)
        if (nodezero.and.util_print('multipole',print_debug)) then
            write(LuOut,9120) '=================='
            call util_flush(LuOut)
        endif
        if (multipole.gt.1) then
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,*)
              write(LuOut,9120) '==== tce_aoq1 ===='
              call util_flush(LuOut)
          endif
          call tce_aoq1(rtdb,g_aoq1,quadrupole_scf)
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,9120) '=================='
              call util_flush(LuOut)
          endif
        endif
        if (multipole.gt.2) then
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,*)
              write(LuOut,9120) '==== tce_aoo1 ===='
              call util_flush(LuOut)
          endif
          call tce_aoo1(rtdb,g_aoo1,octupole_scf)
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,9120) '=================='
              call util_flush(LuOut)
          endif
        endif
        cpu = cpu + util_cpusec()
        wall = wall + util_wallsec()
        dipole_scf_t = dsqrt( dipole_scf(1)*dipole_scf(1) +
     1                        dipole_scf(2)*dipole_scf(2) +
     2                        dipole_scf(3)*dipole_scf(3) )
        if (nodezero.and.util_print('mod1',print_default)) then
          write(LuOut,9430) "SCF",dipole_scf(1),dipole_scf(1)*debye,
     1                            dipole_scf(2),dipole_scf(2)*debye,
     2                            dipole_scf(3),dipole_scf(3)*debye,
     3                            dipole_scf_t ,dipole_scf_t*debye
          if (multipole.gt.1) then
            write(LuOut,9432) 'SCF','second',
     1             quadrupole_scf(1),quadrupole_scf(1)*buck,
     2             quadrupole_scf(2),quadrupole_scf(2)*buck,
     3             quadrupole_scf(3),quadrupole_scf(3)*buck,
     4             quadrupole_scf(4),quadrupole_scf(4)*buck,
     5             quadrupole_scf(5),quadrupole_scf(5)*buck,
     6             quadrupole_scf(6),quadrupole_scf(6)*buck
          endif
          if (multipole.gt.2) then
            write(LuOut,9433) 'SCF','third',
     1             octupole_scf(1), octupole_scf(1)*octo,
     2             octupole_scf(2), octupole_scf(2)*octo,
     3             octupole_scf(3), octupole_scf(3)*octo,
     4             octupole_scf(4), octupole_scf(4)*octo,
     5             octupole_scf(5), octupole_scf(5)*octo,
     6             octupole_scf(6), octupole_scf(6)*octo,
     7             octupole_scf(7), octupole_scf(7)*octo,
     8             octupole_scf(8), octupole_scf(8)*octo,
     9             octupole_scf(9), octupole_scf(9)*octo,
     1             octupole_scf(10),octupole_scf(10)*octo
          endif
          write(LuOut,*)
          write(LuOut,9020) cpu, wall
          call util_flush(LuOut)
        endif
        !
        ! dipole
        !
        do axis = 1, 3
          irrep_d=sym_abelian_axis(geom,axis)
          call sym_irrepname(geom,irrep_d+1,irrepname)
          call tce_mod1_offset(l_d1_offset(axis),k_d1_offset(axis),
     1      size_d1(axis))
          call tce_filename(dfilename(axis),filename)
          call createfile(filename,d_d1(axis),size_d1(axis))
          if (nodezero.and.util_print('mod1',print_default)) then
            write(LuOut,*)
            write(LuOut,9440) axisname(axis),irrepname
            write(LuOut,9000) 'dipole',size_d1(axis)
            write(LuOut,9010) 'dipole',trim(filename)
          endif
        enddo
        if (.not.read_integrals(3)) then
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,*)
              write(LuOut,9120) '==== tce_mod1 ===='
              call util_flush(LuOut)
          endif
          call tce_mod1(g_aod1,d_d1)
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,9120) '=================='
              call util_flush(LuOut)
          endif
          if (save_integrals(3)) then
            if(nodezero) then
              write(LuOut,*) 'Saving dipole integrals now...'
            endif
            do axis = 1, 3
              irrep_d=sym_abelian_axis(geom,axis)
              call d1_restart_save(d_d1(axis),k_d1_offset(axis),
     &             size_d1(axis),axis,handle_d1(axis),irrep_d)
            enddo
          endif
        elseif (read_integrals(3)) then
          if(nodezero) then
            write(LuOut,*) 'Reading in dipole integrals now...'
          endif
          do axis = 1, 3
            if (.not.ga_destroy( g_aod1(axis) ) )
     1        call errquit('tce_energy: GA problem',axis,GA_ERR)
            irrep_d=sym_abelian_axis(geom,axis)
            call d1_restart_read(d_d1(axis),k_d1_offset(axis),
     &           size_d1(axis),axis,handle_d1(axis),irrep_d)
          enddo
        endif
        do axis = 1, 3
          call gatoeaf(d_d1(axis))
        enddo
        call util_flush(LuOut)
        !
        ! quadrupole
        !
        if (multipole.gt.1) then
          do axis = 1, 6
            irrep_d=sym_abelian_axis(geom,axis+3)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            call tce_mod1_offset_new(l_q1_offset(axis),
     1                               k_q1_offset(axis),
     2                               size_q1(axis),irrep_d)
            call tce_filename(qfilename(axis),filename)
            call createfile(filename,d_q1(axis),size_q1(axis))
            if (nodezero.and.util_print('moq1',print_default)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis+3),irrepname
              write(LuOut,9000) 'quadrupole',size_q1(axis)
              write(LuOut,9010) 'quadrupole',trim(filename)
            endif
          enddo
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,*)
              write(LuOut,9120) '==== tce_moq1 ===='
              call util_flush(LuOut)
          endif
          call tce_moq1(g_aoq1,d_q1)
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,9120) '=================='
              call util_flush(LuOut)
          endif
          do axis = 1, 6
            call gatoeaf(d_q1(axis))
          enddo
        endif ! quadrupole
        call util_flush(LuOut)
        !
        ! octupole
        !
        if (multipole.gt.2) then
          do axis = 1, 10
            irrep_d=sym_abelian_axis(geom,axis+9)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            call tce_mod1_offset_new(l_o1_offset(axis),
     1                               k_o1_offset(axis),
     2                               size_o1(axis),irrep_d)
            call tce_filename(ofilename(axis),filename)
            call createfile(filename,d_o1(axis),size_o1(axis))
            if (nodezero.and.util_print('moo1',print_default)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis+9),irrepname
              write(LuOut,9000) 'octupole',size_o1(axis)
              write(LuOut,9010) 'octupole',trim(filename)
            endif
          enddo
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,*)
              write(LuOut,9120) '==== tce_moo1 ===='
              call util_flush(LuOut)
          endif
          call tce_moo1(g_aoo1,d_o1)
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,9120) '=================='
              call util_flush(LuOut)
          endif
          do axis = 1, 10
            call gatoeaf(d_o1(axis))
          enddo
        endif ! octupole
        call util_flush(LuOut)
!!
!!       new dipole integrals - never finished
!!       (save the hnd_mtpole2 reference in case useful later)
!!
!        if (multipole.ge.1) then
!          do i = 1, multipole
!            if (i.eq.1) then
!              call hnd_mtpole2(rtdb,ao_bas_han,geom,1,dipole_scf)
!            elseif (i.eq.2) then
!              call hnd_mtpole2(rtdb,ao_bas_han,geom,2,quadrupole_scf)
!            elseif (i.eq.3) then
!              call hnd_mtpole2(rtdb,ao_bas_han,geom,3,octupole_scf)
!            endif
!            j = (i+1)*(i+2)/2 ! 3, 6, 10
!            do axis = 1, j
!            enddo
!          enddo
!        endif
      endif
!
!     1-e integrals
!
! === jaguar ===
      if(.not.read_in3) then ! --- read_in3 ---
! ===============
      cpu = - util_cpusec()
      wall = - util_wallsec()
      if (.not.rtdb_get(rtdb,'tce:recompf',mt_log,1,recompf))
     1  call errquit('tce_energy: failed reading from rtdb',0,
     2  RTDB_ERR)
      if (recompf.and.(.not.read_integrals(1))) then
        call tce_ao1e_fock2e(rtdb,g_ao1e)
      endif
!
!->d3p975
!      do i = 1, ipol
!        if (.not.ga_destroy(g_movecs(i)))
!     1    call errquit('tce_energy: GA problem',0,GA_ERR)
!      enddo
!<-d3p975
!
      call tce_filename('f1',filename)
      call tce_mo1e_offset(l_f1_offset,k_f1_offset,size_1e)
      call createfile(filename,d_f1,size_1e)
      if (.not.read_integrals(1)) then
        if (recompf) then
! === jaguar ===
      if(read_int)  then
      call util_file_name0('f1int',.false.,.true.,filename,fldgts)
      unitn=77
      if(nodezero) write(LuOut,*)'size_1e',size_1e
       call read_tensor(filename,d_f1,size_1e,unitn)
!     re-split the int files?
       if (rtdb_get(rtdb,'tce:maxioprocs_new',
     M      mt_int,1,max_ioprocs_in_new))  then
          max_ioprocs_in_org=max_ioprocs_in
          max_ioprocs_in=max_ioprocs_in_new
          call write_tensor(filename,d_f1,size_1e,unitn)
          max_ioprocs_in=max_ioprocs_in_org
       endif
       call ga_sync()
       go to 2002
      end if
! ==============
          call tce_mo1e(g_ao1e,d_f1,k_f1_offset)
! === jaguar ===
      if(write_int)  then
      call util_file_name0('f1int',.false.,.true.,filename,fldgts)

      unitn=77
      if(nodezero) write(LuOut,*)'size_1e',size_1e
      call write_tensor(filename,d_f1,size_1e,unitn)
      call ga_sync()
      end if
! ==============
 2002     continue
        else
          call tce_mo1e_epsilon(d_f1)
        endif 
        if (save_integrals(1))then
          if(nodezero) then
            write(LuOut,*) 'Saving 1-electron integrals now...'
          endif
          call f1_restart_save(d_f1,k_f1_offset,size_1e,handle_f1)
        endif
      elseif (read_integrals(1)) then
        if(nodezero) then
          write(LuOut,*) 'Reading in 1-electron integrals now...'
        endif
        call f1_restart_read(d_f1,k_f1_offset,size_1e,handle_f1)
      endif
      call gatoeaf(d_f1)
      cpu = cpu + util_cpusec()
      wall = wall + util_wallsec()
      if (nodezero.and.util_print('mo1e',print_default)) then
        write(LuOut,*)
        if (recompf) then
          write(LuOut,9120) 'Fock matrix recomputed'
        else
          write(LuOut,9120) 'Epsilons used as Fock matrix'
        endif
        write(LuOut,9000) '1-e',size_1e
        write(LuOut,9010) '1-e',trim(filename)
        write(LuOut,9020) cpu, wall
        call util_flush(LuOut)
      endif
#ifdef MIRROR_FOCK
      if (nodezero) write(LuOut,*) 'Mirroring Fock matrix'
!
!     allocate space for local copy
!
      if (.not.ma_push_get(mt_dbl,size_1e,'f1temp',l_f1temp,k_f1temp))
     1    call errquit('tce_energy: MA problem',3737,MA_ERR)
!
!     make a mirror GA for Fock matrix
!
      call tce_filename('f1temp',filename)
      call createfile_mirror(filename,d_f1temp,size_1e)
!
!     copy the entire distributed GA to local buffer
!
      call ga_get(d_f1,1,size_1e,1,1,dbl_mb(k_f1temp),1)
      call ga_sync()
!
!     put it back in the mirrored GA
!
      call ga_put(d_f1temp,1,size_1e,1,1,dbl_mb(k_f1temp),1)
      call ga_sync() ! no need to merge_mirrored since copying entire array
!
!     delete the distributed version
!
      if (.not.ga_destroy(d_f1))
     1    call errquit('tce_energy: GA problem',3737,GA_ERR)
!
!     delete local temporary space
!
      if (.not.ma_pop_stack(l_f1temp))
     1  call errquit('tce_energy: MA problem',3737,MA_ERR)
!
!     reassign the handle so that it points to the mirrored GA
!
      d_f1 = d_f1temp
#endif
! === jaguar ===
      end if ! --- read_in3 ---
! ==============
!
      io_v2=IOGA
! ccsd_act/eomccsd_act
      if((model.eq.'ccsd_act').and.(.not.intorb)) then
       call errquit('tce_energy: ccsd_act with 2eorb only (RHF,ROHF)',
     1               911,GA_ERR)
      end if
      if((model.eq.'ccsd_act').and.
     1   (.not.((fast2e.eq.13).or.(fast2e.eq.14)))) then
       call errquit('tce_energy: ccsd_act 2emet 13 or 14 only',
     1               911,GA_ERR)
      end if
!
      if(.not.intorb) THEN  !--------------
        if (read_integrals(2)) then
          call errquit('tce_energy: cannot restart without 2eorb',
     1                  911,GA_ERR)
        endif
!
!       2-e integrals first half
!
        cpu = - util_cpusec()
        wall = - util_wallsec()
        if(fast2e.eq.1) then
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'tce_ao2e: fast2e=1'
            write(LuOut,*) 'half-transformed integrals in memory'
            call util_flush(LuOut)
          endif
c
        if(.not.read_int) then
          call tce_ao2e(rtdb,d_ao2e)
        end if
c
        elseif((fast2e.eq.2).or.(fast2e.eq.3)) then
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'tce_ao2e_disk: fast2e>1'
            write(LuOut,*) 'half-transformed integrals on disk'
            call util_flush(LuOut)
          endif
          call tce_ao2e_disk(rtdb,d_ao2e)
        else
          write(LuOut,9123)fast2e
          call errquit('tce_energy: invalid 2emet: ',fast2e,CALC_ERR)
        end if
!        elseif(fast2e.eq.4) then
!          if (nodezero) then
!            write(LuOut,*)
!            write(LuOut,*) 'tce_ao2e: fast2e=4'
!            write(LuOut,*) 'half-transformed integrals in memory'
!            call util_flush(LuOut)
!          endif
!!           call tce_ao2e_new(rtdb,geom,ao_bas_han,tol2e,ipol,d_ao2e)
!          call tce_ao2e_new(rtdb,ao_bas_han,d_ao2e)
!        endif
        cpu = cpu + util_cpusec()
        wall = wall + util_wallsec()
        if (nodezero.and.util_print('time',print_default))
     1    write(LuOut,9020) cpu, wall
        if (nodezero) call util_flush(LuOut)
!
!       2-e integrals second half
!
        cpu = - util_cpusec()
        wall = - util_wallsec()
        call tce_mo2e_offset(l_v2_offset,k_v2_offset,size_2e)
        if(fast2e.eq.1) then
          call tce_filename('v2',filename)
          call createfile(filename,d_v2,size_2e)
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'tce_mo2e: fast2e=1'
            write(LuOut,*) '2-e integrals stored in memory'
            call util_flush(LuOut)
          endif
! === jaguar === SPINORB
         if(read_int)  then
c          call util_file_name('v2int',.false.,.true.,filename)
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=78
          call read_tensor(filename,d_v2,size_2e,unitn)
!     re-split the int files?
          if (rtdb_get(rtdb,'tce:maxioprocs_new',
     M         mt_int,1,max_ioprocs_in_new))  then
             max_ioprocs_in=max_ioprocs_in_new
             call write_tensor(filename,d_v2,size_2e,unitn)
          endif
          call ga_sync()
          go to 9003
         end if
! ==============
          call tce_mo2e(rtdb,d_ao2e,d_v2,k_v2_offset)
          call reconcilefile(d_v2,size_2e)
! === jaguar ===
         if(write_int)  then
c          call util_file_name('v2int',.false.,.true.,filename)
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=78
          call write_tensor(filename,d_v2,size_2e,unitn)
          call ga_sync()
         end if
! ==============
 9003    continue
        elseif(fast2e.eq.2) then
          call tce_filename('v2',filename)
          call createfile(filename,d_v2,size_2e)
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'tce_mo2e_disk: fast2e=2'
            write(LuOut,*) '2-e integrals stored in memory'
            call util_flush(LuOut)
          endif
          call tce_mo2e_disk(rtdb,d_ao2e,d_v2,k_v2_offset)
        elseif(fast2e.eq.3) then
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'tce_mo2e_disk: fast2e=3'
            write(LuOut,*) '2-e integrals stored on disk'
            call util_flush(LuOut)
          endif
          call tce_filename_sf('v2',filename)
          call createfile_sf(filename,d_v2,size_2e)
          io_v2 = IOSF
          d_v2orb = d_v2
          call tce_zero_sf(d_v2,size_2e)
          call tce_mo2e_disk(rtdb,d_ao2e,d_v2,k_v2_offset)
        elseif(fast2e.eq.4) then
          call tce_filename('v2',filename)
          call createfile(filename,d_v2,size_2e)
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'tce_mo2e: fast2e=4'
            write(LuOut,*) '2-e integrals stored in memory'
            call util_flush(LuOut)
          endif
          call tce_mo2e_new(d_ao2e,d_v2,k_v2_offset)
          call reconcilefile(d_v2,size_2e)
        endif
        cpu = cpu + util_cpusec()
        wall = wall + util_wallsec()
        if (nodezero.and.util_print('mo2e',print_default)) then
          write(LuOut,*)
          write(LuOut,9000) '2-e',size_2e
          write(LuOut,9010) '2-e',trim(filename)
!          write(LuOut,9090) '2-e',d_v2
          write(LuOut,9020) cpu, wall
          call util_flush(LuOut)
        endif
      elseif(intorb) THEN !----------------------
! === jaguar ===
        if(.not.read_in3) then ! -- read_in3 ---
! ==============
        if(nodezero) then
          write(LuOut,*)'4-electron integrals stored in orbital form'
        end if
        cpu = - util_cpusec()
        wall = - util_wallsec()
        length = 0
        if(fast2e.eq.12) then
          call tce_mo2e_offset_intorb_hybrid(lax_v2_alpha_offset_ga,
     &         kax_v2_alpha_offset_ga,size_2e_ga)
          call tce_mo2e_offset_intorb(lax_v2_alpha_offset_sf,
     &         kax_v2_alpha_offset_sf,size_2e_sf)
        else
         if(model.ne.'ccsd_act') then
          call tce_mo2e_offset_intorb(lax_v2_alpha_offset,
     &         kax_v2_alpha_offset,size_2e)
         else
! ccsd_act/eomccsd_act
          call tce_mo2e_offset_intorb_act(lax_v2_alpha_offset,
     &         kax_v2_alpha_offset,size_2e)
! --------------------
         end if
        endif
        l_v2_alpha_offset=lax_v2_alpha_offset
        k_v2_alpha_offset=kax_v2_alpha_offset
        if (fast2e.eq.9) then
           call tce_filename_sf('v2',filename)
           call createfile_sf_suffix(filename,d_v2,size_2e,V2_SUFFIX)
           io_v2 = IOSF
           d_v2orb = d_v2
           call tce_zero_sf(d_v2,size_2e)
        elseif (fast2e.eq.12) then
           io_v2 = IOHYBRID
!          The GA copy of 2-e integrals
           call tce_filename('v2ga',filename)
           call createfile(filename,d_v2ga,size_2e_ga)
           call ga_zero(d_v2ga)
           d_v2 = d_v2ga
!          The SF copy of 2-e integrals
           call tce_filename_sf('v2sf',filename)
           call createfile_sf_suffix(filename,d_v2sf,size_2e_sf,
     1                               V2_SUFFIX)
           call tce_zero_sf(d_v2sf,size_2e_sf)
           print*,'=========================='
           print*,'d_v2    = ',d_v2
           print*,'d_v2orb = ',d_v2orb
           print*,'d_v2ga  = ',d_v2ga
           print*,'d_v2sf  = ',d_v2sf
           print*,'=========================='
        else
           call tce_filename('v2',filename)
           call createfile(filename,d_v2,size_2e)
           call ga_zero(d_v2)
        endif
        if (fast2e.eq.12) then
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9000) 'v2ga ',size_2e_ga
            write(LuOut,9000) 'v2sf ',size_2e_sf
            call util_flush(LuOut)
          endif
        else
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9000) 'v2   ',size_2e
            call util_flush(LuOut)
          endif
        endif
        d_v2orb = d_v2
        cpu = - util_cpusec()
        wall = - util_wallsec()
        if(nodezero) then
          write(LuOut,9123)fast2e
          write(LuOut,9124)imaxsize
          write(LuOut,9125)idisk
        end if
        if(fast2e.ge.2) call atpart_nalength
!
!       ============================================================
!
        if((.not.read_integrals(2)).or.(fast2e.eq.9)) then
!
        if(fast2e.eq.1) then
          call tce_mo2e_trans(rtdb,d_v2,kax_v2_alpha_offset,size_2e)
        elseif(fast2e.eq.2) then
          call tce_mo2e_zones_4a_disk_2s_new(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
        elseif(fast2e.eq.3) then
          call tce_mo2e_zones_4a_disk_ga(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
        elseif(fast2e.eq.4) then
          call  tce_mo2e_zones_4a_disk_ga_chop(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
        elseif(fast2e.eq.5) then
          idisk=1
          if(idisk.ne.1) then
            call errquit('tce_energy: 2emet 5 with idiskx=1 only',
     2                                             idisk,CALC_ERR)
          endif
          call tce_mo2e_zones_4a_disk_4s_new_patch(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
        elseif(fast2e.eq.6) then
          idisk=1
          if(idisk.ne.1) then
            call errquit('tce_energy: 2emet 6 with idiskx=1 only',
     2                                             idisk,CALC_ERR)
          endif
          call tce_mo2e_zones_4a_disk_4s_new_patch_m2(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
!        elseif(fast2e.eq.7) then
!          if(idisk.ne.1) then
!            call errquit('tce_energy: 2emet 7 with idiskx=1 only',
!     2                                             idisk,CALC_ERR)
!          endif
!          call tce_mo2e_zones_4a_disk_4s_new_patch_m3(rtdb,d_v2,
!     1                                  kax_v2_alpha_offset,size_2e)
!        elseif(fast2e.eq.8) then
!          if(idisk.ne.1) then
!            call errquit('tce_energy: 2emet 8 with idiskx=1 only',
!     2                                             idisk,CALC_ERR)
!          endif
!          call tce_mo2e_zones_4a_disk_4s_new_patch_m4(rtdb,d_v2,
!     1                                  kax_v2_alpha_offset,size_2e)
        elseif ((fast2e.eq.9).and.(.not.(read_integrals(2)))) then
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'tce_mo2e_disk_2eorb: fast2e=9'
            call util_flush(LuOut)
          endif
          call tce_mo2e_disk_2eorb(rtdb,d_v2,kax_v2_alpha_offset,
     1                             size_2e)
        elseif ((fast2e.eq.9).and.(read_integrals(2))) then
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'tce_mo2e_disk_2eorb: fast2e=9'
            write(LuOut,*) 'restarting with integral file on disk'
            call util_flush(LuOut)
          endif
        elseif(fast2e.eq.10) then
          idisk=0
          if(idisk.ne.0) then
            call errquit('tce_energy: 2emet 10 with idiskx=0 only',
     2                                             idisk,CALC_ERR)
          endif
          call tce_mo2e_incore_2eorb(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
        elseif(fast2e.eq.11) then
          idisk=0
          if(idisk.ne.0) then
            call errquit('tce_energy: 2emet 11 with idiskx=0 only',
     2                                             idisk,CALC_ERR)
          endif
          call tce_mo2e_incore_2eorb_split(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
        elseif(fast2e.eq.12) then
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'fast2e=12'
            write(LuOut,*) '2-e integrals stored in GA and SF'
            write(LuOut,*)
            call util_flush(LuOut)
          endif
!          if (nodezero) write(LuOut,*) 'tce_mo2e_incore_2eorb'
!          call tce_mo2e_incore_2eorb(rtdb,d_v2ga,kax_v2_alpha_offset,
!     1                               size_2e)
!          if (nodezero) then
!            write(LuOut,*) 'tce_mo2e_zones_4a_disk_4s_new_patch'
!          endif
!          call tce_mo2e_zones_4a_disk_4s_new_patch(rtdb,d_v2ga,
!     1                                  kax_v2_alpha_offset,size_2e)
!          if (nodezero) write(LuOut,*) 'tce_mo2e_disk_2eorb'
!          idisk=.true.
!          call tce_mo2e_disk_2eorb(rtdb,d_v2sf,kax_v2_alpha_offset,
!     1                             size_2e)
          if (nodezero) write(LuOut,*) 'tce_mo2e_hybrid_2eorb_split'
          idisk=1
          call tce_mo2e_hybrid_2eorb_split(d_v2ga,d_v2sf,
     1         kax_v2_alpha_offset)
        elseif(fast2e.eq.13) then
! === jaguar ===
         if(read_int)  then
c          call util_file_name('v2int',.false.,.true.,filename)
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=78
          call read_tensor(filename,d_v2,size_2e,unitn)
!     re-split the int files?
          if (rtdb_get(rtdb,'tce:maxioprocs_new',
     M         mt_int,1,max_ioprocs_in_new))  then
             max_ioprocs_in=max_ioprocs_in_new
             call write_tensor(filename,d_v2,size_2e,unitn)
          endif
          call ga_sync()
          go to 3003
         end if
! ==============
         if(model.ne.'ccsd_act') then
          call tce_mo2e_zones_4a_disk_ga_N5(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
!
!
!
!
!
!
      if(pwtce) then
          if (.not.(hftype.eq.'RHF'))
     1     call errquit('tce_energy:hnormal only with RHF reference',0,
     2     RTDB_ERR)
           do g1b=1,noab+nvab
            if(int_mb(k_range+g1b-1).ne.1) then
              call errquit('tce_energy:pwtce tilesize not equal 1',0,
     2        RTDB_ERR)
            endif
           enddo
           nos=int_mb(k_offset+noab-1)+int_mb(k_range+noab-1)
           nvs=int_mb(k_offset+noab+nvab-1)+int_mb(k_range+noab+nvab-1)
     &         -nos
           noas=nos/2
           nobs=nos/2
           nvas=nvs/2
           nvbs=nvs/2
           mov=noas+nvas
           call pw_tce(rtdb,d_v2,k_v2_alpha_offset,d_f1,k_f1_offset,
     &       dbl_mb(k_evl_sorted),nos,nvs,noas,nobs,nvas,nvbs,mov,ehfa)
c only nodezero has ref
           if (.not.rtdb_get(rtdb,'tce:pwref',mt_dbl,1,ref))
     &     call errquit('pwtce: failed to write ehfa',0,RTDB_ERR)
c *** debug ***
c        if(nodezero) then
c         write(6,*)'ehfa pwtce: ',ref
c         call util_flush(6)
c        endif
c *************
      endif
!
!
!
!
!
!QuantumComputing
!
           if(quantumc) then
!ehfa
           if (.not. geom_nuc_rep_energy(geom, enrep_tce))
     &      call errquit('TCE:unable to get nuclear rep ',geom, 
     &      GEOM_ERR)
c
c
           range_ov=int_mb(k_offset_alpha+noa+nva-1)
     1          +int_mb(k_range_alpha+noa+nva-1)
           range_occ=int_mb(k_offset_alpha+noa-1)
     1          +int_mb(k_range_alpha+noa-1)
c
           size_f2m=range_ov*range_ov
           size_v4m=range_ov*range_ov*range_ov*range_ov
c
       if(.not.ma_push_get(mt_dbl,size_f2m,'f2m',l_f2m,k_f2m))
     1     call errquit('tce_f2m: MA problem',0,MA_ERR)
       if(.not.ma_push_get(mt_dbl,size_v4m,'v4m',l_v4m,k_v4m))
     1     call errquit('tce_v4m: MA problem',0,MA_ERR)
c
           spust=.true.
c
      if(nodezero) then
       write(6,*)'begin_one_electron_integrals'
       call util_flush(6)
      endif
c
           call kinetic_hcore(rtdb,dbl_mb(k_f2m),range_ov,spust)
           call util_flush(6)
      if(nodezero) then
       write(6,*)'end_one_electron_integrals'
       call util_flush(6)
      endif
c
      if(nodezero) then  !nodezero 
c      write(6,*)'range_ov  = ',range_ov
c      write(6,*)'range_occ = ',range_occ
       write(6, *) 'begin_two_electron_integrals'
       call util_flush(6)
      offset_test=0
      DO g3b = 1,noa+nva   !j
      DO g4b = g3b,noa+nva !i
      DO g1b = 1,noa+nva   !l
      DO g2b = g1b,noa+nva !k
      IF (int_mb(k_spin_alpha+g3b-1)+int_mb(k_spin_alpha+g4b-1).eq.
     &int_mb(k_spin_alpha+g1b-1)+int_mb(k_spin_alpha+g2b-1)) THEN
      IF (ieor(int_mb(k_sym_alpha+g3b-1),ieor(int_mb(k_sym_alpha+g4b-1),
     &    ieor(int_mb(k_sym_alpha+g1b-1),int_mb(k_sym_alpha+g2b-1))))
     &    .eq. irrep_v) THEN
      IROW=INDEX_PAIR(g4b,g3b)
      ICOL=INDEX_PAIR(g2b,g1b)
      IF(IROW.GE.ICOL) THEN
      size_test=int_mb(k_range_alpha+g3b-1) *
     &int_mb(k_range_alpha+g4b-1) * int_mb(k_range_alpha+g1b-1) *
     &int_mb(k_range_alpha+g2b-1)
c
       if(.not.ma_push_get(mt_dbl,size_test,'test',l_aux,k_aux))
     1     call errquit('tce_aux2: MA problem',0,MA_ERR)
        call ga_get(d_v2,offset_test+1,offset_test+size_test,1,1,
     1    dbl_mb(k_aux),size_test)
      i=0
      do igi1=1,int_mb(k_range_alpha+g1b-1)
      do igi2=1,int_mb(k_range_alpha+g2b-1)
      do igi3=1,int_mb(k_range_alpha+g3b-1)
      do igi4=1,int_mb(k_range_alpha+g4b-1)
       i=i+1
       xxx=dabs(dbl_mb(k_aux+i-1))
       tg1=igi1+int_mb(k_offset_alpha+g1b-1)
       tg2=igi2+int_mb(k_offset_alpha+g2b-1)
       tg3=igi3+int_mb(k_offset_alpha+g3b-1)
       tg4=igi4+int_mb(k_offset_alpha+g4b-1)
       tirow=index_pair(tg4,tg3)
       ticol=index_pair(tg2,tg1)
c
       if((tg3.le.tg4).and.(tg1.le.tg2)) then
       if(tirow.ge.ticol) then
         if(xxx.gt.1.0d-10) then
c yaml  two_electron_integrals
          write(6,4881) tg4,
     &                  tg3,
     &                  tg2,
     &                  tg1,
     7              dbl_mb(k_aux+i-1)
         endif !xxx
        else
         if(irow.ne.icol) then
         if(xxx.gt.1.0d-10) then
c yaml  two_electron_integrals
          write(6,4881) tg2,
     &                  tg1,
     &                  tg4,
     &                  tg3,
     7              dbl_mb(k_aux+i-1)
c
c          write(6,*)'--------'
         endif !xxx
         endif
       endif !row-column
       endif !tg3<=tg4 tg1<=tg2
      enddo
      enddo
      enddo
      enddo
c
      call util_flush(6)
c
      if (.not.ma_pop_stack(l_aux))
     1     call errquit('l_aux2',15,MA_ERR)
c
      offset_test=offset_test+size_test

      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
c     ehfa
      if(nodezero) then
         write(6, *) 'end_two_electron_integrals'
         call util_flush(6)
      endif
      ehfa=0.0d0
      endif !nodezero
c
      call active_shift(rtdb,d_v2,dbl_mb(k_v4m),dbl_mb(k_f2m),range_occ,
     &                 range_ov,ehfa)
c
      if(nodezero) then
c yaml enerep_tce - coulomb_repulsion
       write(6,*)'enrep_tce = ',enrep_tce
c ehfa - HF energy for active space - not used in yaml 
       write(6,6887) ehfa
c yaml red - this is hf_energy
       write(6,6888) ref
c yaml scf_energy_offset  (ref-ehfa) same and energy_offset
       write(6,6889) ref-ehfa
c
      call util_flush(6)
c
c
c      l1norm=0.0d0
c      do igi1=1,range_ov
c      do igi2=1,range_ov
c       xxx=0.0d0
c       do igi3=1,range_ov
c       do igi4=1,range_ov
c        xxx=xxx+dabs(v4m(igi4,igi3,igi2,igi1))
c       enddo
c       enddo
c       if(xxx.gt.l1norm) l1norm=xxx
c      enddo
c      enddo
c
c      write(6,*) 'L1 norm of v2 = ',l1norm
c
      endif !nodezero
c
      if (.not.ma_pop_stack(l_v4m))
     1     call errquit('l_v4m',15,MA_ERR)
      if (.not.ma_pop_stack(l_f2m))
     1     call errquit('l_f2m',15,MA_ERR)
c
      endif ! quantumc
!
!QuantumComputing
!
!
!
         else
! ccsd_act/eomccsd_act
          call tce_mo2e_zones_4a_disk_ga_act(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
! --------------------
         end if
! === jaguar ===
         if(write_int)  then
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=78
          call write_tensor(filename,d_v2,size_2e,unitn)
          call ga_sync()
         end if
! ==============
 3003    continue
!
        elseif(fast2e.eq.14) then
! === jaguar ===
         if(read_int)  then
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=78
          call read_tensor(filename,d_v2,size_2e,unitn)
          call ga_sync()
          go to 3004
         end if
! ==============
         if(model.ne.'ccsd_act') then
            write(6,*) ' chop v2 '
          call  tce_mo2e_zones_4a_disk_ga_chop_N5(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
         else 
! ccsd_act/eomccsd_act
          call  tce_mo2e_zones_4a_disk_ga_chop_act(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
! --------------------
         end if
! === jaguar ===
         if(write_int)  then
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=78
          call write_tensor(filename,d_v2,size_2e,unitn)
          call ga_sync()
         end if
! ==============
 3004    continue
!
        elseif(fast2e.eq.15) then
! === jaguar ===
         if(read_int)  then
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=178
          call read_tensor(filename,d_v2,size_2e,unitn)
          call ga_sync()
!     re-split the int files?
          if (rtdb_get(rtdb,'tce:maxioprocs_new',
     M         mt_int,1,max_ioprocs_in_new))  then
             max_ioprocs_in=max_ioprocs_in_new
             call write_tensor(filename,d_v2,size_2e,unitn)
          endif
          go to 3095
         end if
! ==============
          call  new_ga4ind_N5(rtdb,d_v2,
     1                        kax_v2_alpha_offset,
     1                        size_2e)
! === jaguar ===
         if(write_int)  then
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=178
          call write_tensor(filename,d_v2,size_2e,unitn)
          call ga_sync()
         end if
! ==============
 3095    continue
c
        elseif(fast2e.eq.16) then
! === jaguar ===
         if(read_int)  then
c          call util_file_name('v2int',.false.,.true.,filename)
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=78
          call read_tensor(filename,d_v2,size_2e,unitn)
!     re-split the int files?
          if (rtdb_get(rtdb,'tce:maxioprocs_new',
     M         mt_int,1,max_ioprocs_in_new))  then
             max_ioprocs_in=max_ioprocs_in_new
             call write_tensor(filename,d_v2,size_2e,unitn)
          endif
          call ga_sync()
          go to 3743
         end if
! ==============
         if(model.ne.'ccsd_act') then
c          call tce_mo2e_zones_4a_disk_ga_N5(rtdb,d_v2,
c     1                                  kax_v2_alpha_offset,size_2e)

          call new4ind_N5_osu(rtdb,d_v2,
     1                                kax_v2_alpha_offset,
     1                                size_2e)   

         end if
! === jaguar ===
         if(write_int)  then
          call util_file_name0('v2int',.false.,.true.,filename,fldgts)
          unitn=78
          call write_tensor(filename,d_v2,size_2e,unitn)
          call ga_sync()
         end if
! ==============
 3743    continue
c
        else
          write(LuOut,9123)fast2e
          call errquit('tce_energy: invalid 2emet: ',fast2e,CALC_ERR)
        end if
        cpu = cpu + util_cpusec()
        wall = wall + util_wallsec()
        if (nodezero) then
          write(LuOut,9020) cpu, wall
          call util_flush(LuOut)
        endif
        if (save_integrals(2).and.(.not.(fast2e.eq.9))) then
          if(nodezero) then
            write(LuOut,*) 'Saving 2-electron integrals now...'
          endif
          call v2_restart_save(d_v2,k_v2_alpha_offset,size_2e,
     &                         handle_v2)
        elseif (save_integrals(2).and.(fast2e.eq.9)) then
          if(nodezero) then
            write(LuOut,*) 'Restarting with 2emet=9 is trivial'
            write(LuOut,*) 'assuming the file is still available'
          endif
        endif
!
!       ============================================================
!
        elseif (read_integrals(2).and.(.not.(fast2e.eq.9))) then
          if(nodezero) then
            write(LuOut,*) 'Reading in 2-electron integrals now...'
          endif
          call v2_restart_read(d_v2,k_v2_alpha_offset,size_2e,
     &                         handle_v2)
        endif
!
!       ============================================================
!
        call tce_mo2e_offset_size(length)
        call e2_offset_2_offset(lx_o2o_offset,kx_o2o_offset,
     &       lx_o2ot_offset,kx_o2ot_offset,
     &       lx_v2spin_offset,kx_v2spin_offset,length)
        l_o2o_offset    = lx_o2o_offset
        k_o2o_offset    = kx_o2o_offset
        l_o2ot_offset   = lx_o2ot_offset
        k_o2ot_offset   = kx_o2ot_offset
        l_v2spin_offset = lx_v2spin_offset
        k_v2spin_offset = kx_v2spin_offset
        k_v2_offset = k_v2spin_offset
! === jaguar ===
        end if  ! read_in3
! ==============
      end if ! intorb -----------------------
c
c
c
c
c  Tsources
c  for UCC(4) model='ccsdt'
c
      if(intorb.and.((model.eq.'ccsdt').or.(model.eq.'ccsdtq'))) then
        if(nodezero) then
         write(6,*)'IN MAPPING 1 for UCC(4),CCSDT,CCSDTQ'
         call util_flush(6)
        end if
        call tce_mo2e_offset(l_v2aux_offset,k_v2aux_offset,size_2eaux)
        call tce_filename('v2aux',filename)
        call createfile(filename,d_v2aux,size_2eaux)
        call mapping_v2(rtdb,d_v2aux,k_v2aux_offset,
     &             d_v2,k_v2_offset)
      endif
c end of Tsource
c
c
c
c
cc =====   Sigma2 ======
!!!      if(sig2) then
!!!          omega = 0.010d0
!!!          eta   = 0.001d0
!!!          size_s2=(nmo(1)+nmo(ipol))*(nmo(1)+nmo(ipol))
!!!          if(.not.MA_PUSH_GET(mt_dbl,size_s2,'s2_re',l_s2_re,k_s2_re))
!!!     1      call errquit('k_s2_re',1,MA_ERR)
!!!          if(.not.MA_PUSH_GET(mt_dbl,size_s2,'s2_im',l_s2_im,k_s2_im))
!!!     1       call errquit('k_s2_im',1,MA_ERR)
!!!          if (.not.MA_PUSH_GET(mt_dbl,size_s2,'s2_aux',l_fs2_a,k_fs2_a))
!!!     1       call errquit('k_s2_aux',1,MA_ERR)
!!!c
!!!       call  sigma2(rtdb,d_v2,k_v2_offset,omega,eta,
!!!     &          dbl_mb(k_s2_re),dbl_mb(k_s2_im))
!!!c
!!!         if(.not.MA_POP_STACK(l_fs2_a))
!!!     &      call errquit('l_fs2_a',4,MA_ERR)
!!!         if(.not.MA_POP_STACK(l_s2_im))
!!!     &      call errquit('l_s2_im',4,MA_ERR)
!!!         if(.not.MA_POP_STACK(l_s2_re))
!!!     &      call errquit('l_fs2_a',4,MA_ERR)
!!!      endif ! sig2
cc === end of Sigma2 ===
c
c
c
c
c
c
!
!     ====================
!     Initial t amplitudes
!     ====================
!
      needt1 = .false.
      needt2 = .false.
      needt3 = .false.
      needt3a= .false.
      needt4 = .false.
      needu1 = .false.
      needu2 = .false.
      needw1 = .false.
      needw2 = .false.
      needw3 = .false.
      needw4 = .false.
      neede2 = .false.
      if ((model.eq.'ccd').or.(model.eq.'lccd')) then
        needt2 = .true.
      elseif (((model.eq.'ccsd').or.(model.eq.'lccsd')).or.
     1        ((model.eq.'qcisd').or.(model.eq.'cisd'))) then
        needt1 = .true.
        needt2 = .true.
! ccsd_act/eomccsd_act
      elseif(model.eq.'ccsd_act') then
        needt1 = .true.
        needt2 = .true.
      elseif ((model.eq.'ccsdt').or.(model.eq.'cisdt')) then
#if !defined(CCSDT)
      CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
        needt1 = .true.
        needt2 = .true.
        needt3 = .true.
#endif
      elseif (model.eq.'ccsdta') then
#if !defined(CCSDT)
      CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
        needt1  = .true.
        needt2  = .true.
        needt3a = .true.
#endif
      elseif ((model.eq.'ccsdtq').or.(model.eq.'cisdtq')) then
#if !defined(CCSDTQ)
      CALL ERRQUIT('CCSDTQ not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
        needt1 = .true.
        needt2 = .true.
        needt3 = .true.
        needt4 = .true.
#endif
      elseif ((model.eq.'mbpt2').or.
     1         (model.eq.'mbpt3').or.
     2         (model.eq.'mbpt4').or.
     2         (model.eq.'mbpt4sdq_t').or.
     3         (model.eq.'mbpt4sdq')) then
        do_mbpt = .true.
        needu1 = .true.
        needu2 = .true.
        needt1 = .true.
        needt2 = .true.
        if ((model.eq.'mbpt3').or.(model.eq.'mbpt4').or.
     1      (model.eq.'mbpt4sdq_t').or.(model.eq.'mbpt4sdq')) then
          needw1 = .true.
          needw2 = .true.
          if (model.eq.'mbpt4') then
            neede2 = .true.
            needw3 = .true.
            needt3 = .true.
          endif ! mbpt4
          if ((model.eq.'mbpt4sdq').or.(model.eq.'mbpt4sdq_t')) then
            neede2 = .true.
          endif ! mbpt4sdq
        endif ! mbpt3
      elseif (model.eq.'multi') then
        needt1 = .true.
        needt2 = .true.
      endif ! model
      if (.not.rtdb_get(rtdb,'tce:do_ccd',mt_log,1,do_ccd)) then
        do_ccd = .false.
      endif
      if (.not.rtdb_get(rtdb,'tce:do_lccd',mt_log,1,do_lccd)) then
        do_lccd = .false.
      endif
      if (.not.rtdb_get(rtdb,'tce:do_ccsd',mt_log,1,do_ccsd)) then
        do_ccsd = .false.
      endif
      if (.not.rtdb_get(rtdb,'tce:do_lccsd',mt_log,1,do_lccsd)) then
        do_lccsd = .false.
      endif
      if (.not.rtdb_get(rtdb,'tce:do_qcisd',mt_log,1,do_qcisd)) then
        do_qcisd = .false.
      endif
      if (.not.rtdb_get(rtdb,'tce:do_pt',mt_log,1,do_pt)) then
        do_pt = .false.
      endif
      if (.not.rtdb_get(rtdb,'tce:do_lam_pt',mt_log,1,do_lam_pt)) then
        do_lam_pt = .false.
      endif
      if (.not.rtdb_get(rtdb,'tce:do_cr_pt',mt_log,1,do_cr_pt)) then
        do_cr_pt = .false.
      endif
      if (.not.rtdb_get(rtdb,'tce:do_lcr_pt',mt_log,1,do_lcr_pt)) then
        do_lcr_pt = .false.
      endif
      if (.not.rtdb_get(rtdb,'tce:do_2t_pt',mt_log,1,do_2t_pt)) then
        do_2t_pt = .false.
      endif
!
!     t1 amplitudes
!
      if (needt1) then
        cpu = - util_cpusec()
        call tce_filename('t1',filename)
        call tce_t1_offset(l_t1_offset,k_t1_offset,size_t1)
! -------
        if(nodezero) then
         write(LuOut,*)'T1-number-of-tasks',int_mb(k_t1_offset)
         call util_flush(LuOut)
        end if
! -------------
#ifdef MIRROR_SING
        call createfile_mirror(filename,d_t1,size_t1)
#else
        call createfile(filename,d_t1,size_t1)
#endif
! ccsd_act/eomccsd_act
        if(model.eq.'ccsd_act') then
          call ga_zero(d_t1)
          go to 8711
        end if
! === jaguar ===
        if(.not.read_in3) then !---
! ==============
        if (.not.read_t(1)) then
          call tce_guess_t1(d_t1,k_t1_offset)
! === jaguar ===
        end if !---
! ==============
#ifdef MIRROR_SING
          call reconcilefile(d_t1)
#endif
        elseif (read_t(1)) then
          call x1_restart_read(d_t1,k_t1_offset,size_t1,0,
     1                         handle_t1,irrep_t)
        endif
!c        call tce_guess_t1_zero(d_t1,k_t1_offset)
! ccsd_act/eomccsd_act
 8711   continue
        call reconcilefile(d_t1,size_t1)
        cpu = cpu + util_cpusec()
        if (nodezero.and.util_print('t1',print_default)) then
          write(LuOut,*)
          write(LuOut,9000) 't1',size_t1
          write(LuOut,9010) 't1',trim(filename)
          write(LuOut,9090) 't1',d_t1
!         write(LuOut,9020) cpu
          call util_flush(LuOut)
        endif
        if (left) then
          irrep_y=0
          call tce_y1_offset(l_l1_offset,k_l1_offset,size_l1)
          call tce_filename('lambda1',filename)
          call createfile(filename,d_lambda1,size_l1)
          if (read_l(1)) then
            call y1_restart_read(d_lambda1,k_l1_offset,size_l1,0,
     1                           handle_l1,irrep_y)
          endif
          call gatoeaf(d_lambda1)
          if (nodezero.and.util_print('l1',print_debug)) then
            write(LuOut,*)
            write(LuOut,9000) 'l1',size_l1
            write(LuOut,9090) 'l1',d_lambda1
            call util_flush(LuOut)
          endif
        endif
        if (lineresp) then
          do axis = 1, 3
          if (respaxis(axis)) then
            irrep_x=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_x+1,irrepname)
            call tce_x1_offset(l_tr1_offset(axis),k_tr1_offset(axis),
     &           size_tr1(axis))
            call tce_filename(tr1filename(axis),filename)
            call createfile(filename,d_tr1(axis),size_tr1(axis))
            if (read_tr(1)) then
              call tr1_restart_read(d_tr1(axis),k_tr1_offset(axis),
     1             size_tr1(axis),axis,handle_tr1(axis),irrep_x)
            endif
            call gatoeaf(d_tr1(axis))
            call tce_filename(tr1filename(axis+3),filename)
            call createfile(filename,d_tr1(axis+3),
     &           size_tr1(axis))
            if (read_tr(4+1)) then
              call tr1_restart_read(d_tr1(axis+3),k_tr1_offset(axis),
     1          size_tr1(axis),axis+3,handle_tr1(axis+3),irrep_x)
            elseif (read_tr(1)) then
              call tr1_restart_read(d_tr1(axis+3),k_tr1_offset(axis),
     1          size_tr1(axis),axis,handle_tr1(axis+3),irrep_x)
            endif
            call gatoeaf(d_tr1(axis+3))
            if (nodezero.and.util_print('tr1',print_debug)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
              write(LuOut,9000) 'tr1',size_tr1(axis)
              write(LuOut,9090) 'tr1',d_tr1(axis)
              write(LuOut,9090) 'tr1',d_tr1(axis+3)
              call util_flush(LuOut)
            endif
          endif
          enddo ! axis
        endif ! lineresp
        if (leftresp) then
          do axis = 1, 3
          if (respaxis(axis)) then
            irrep_y=sym_abelian_axis(geom,axis)
            call tce_y1_offset(l_yr1_offset(axis),k_yr1_offset(axis),
     &           size_yr1(axis))
            call tce_filename(yr1filename(axis),filename)
            call createfile(filename,d_yr1(axis),size_yr1(axis))
            if (read_lr(1)) then
              call yr1_restart_read(d_yr1(axis),k_yr1_offset(axis),
     1          size_yr1(axis),axis,handle_yr1(axis),irrep_y)
            endif
            call gatoeaf(d_yr1(axis))
            call tce_filename(yr1filename(axis+3),filename)
            call createfile(filename,d_yr1(axis+3),
     &           size_yr1(axis))
            if (read_lr(1)) then
              call yr1_restart_read(d_yr1(axis+3),k_yr1_offset(axis),
     1          size_yr1(axis),axis,handle_yr1(axis+3),irrep_y)
            endif
            call gatoeaf(d_yr1(axis+3))
          endif
          enddo ! axis
        endif ! leftresp
      endif
!
!     t2 amplitudes
!
      if (needt2) then
!--       if(.not.intorb) then ! ----------
        cpu = - util_cpusec()
        call tce_filename('t2',filename)
        if(model.ne.'ccsd_act') then
        call tce_t2_offset(l_t2_offset,k_t2_offset,size_t2)
        else
! ccsd_act/eomccsd_act
         call tce_t2a_offset(l_t2_offset,k_t2_offset,size_t2)
! --------------------
        end if
! ------------
        if(nodezero) then
         write(LuOut,*)'T2-number-of-boxes',int_mb(k_t2_offset)
         call util_flush(LuOut)
        end if
! -------------
        call createfile(filename,d_t2,size_t2)
        if(nodezero) then
         call util_flush(LuOut)
        end if
        
! ccsd_act/eomccsd_act
        if(model.eq.'ccsd_act') go to 8712
        if (.not.read_t(2)) then
! === jaguar ===
        if(.not.read_in3) then ! --- read_in3
! ==============
          call tce_guess_t2(d_v2,k_v2_offset,d_t2,k_t2_offset)
          call reconcilefile(d_t2,size_t2)
! starts for t1 amplitudes t1(in the second order of MBPT)
          if(needt1.and.((model.eq.'ccsd').or.(model.eq.'multi'))) then
            if ((ioalg.eq.2).and.(io_v2.eq.IOGA)) then
              call t1mp2(d_f1,d_t1,d_t2,d_v2,k_f1_offset,k_t1_offset,
     &                                       k_t2_offset,k_v2_offset)
              call reconcilefile(d_t1,size_t1)
              call tce_guess_t1(d_t1,k_t1_offset)
              call reconcilefile(d_t1,size_t1)
            end if
          end if
! === jaguar ===
        end if ! --- read_in3
! ==============
        elseif (read_t(2)) then
          call x2_restart_read(d_t2,k_t2_offset,size_t2,0,
     1                         handle_t2,irrep_t)
        endif
        if (util_print('mbpt2',print_debug)) then
          call tce_mbpt2(d_v2,k_v2_offset,d_t2,k_t2_offset,corr)
          if (nodezero) then
            write(LuOut,9030) corr
            write(LuOut,9040) ref + corr
          endif
        endif
 8712   continue
        cpu = cpu + util_cpusec()
        if (nodezero.and.util_print('t2',print_default)) then
          write(LuOut,*)
          write(LuOut,9000) 't2',size_t2
          write(LuOut,9010) 't2',trim(filename)
          write(LuOut,9090) 't2',d_t2
!         write(LuOut,9020) cpu
          call util_flush(LuOut)
        endif
        if (left) then
          irrep_y=0
          call tce_y2_offset(l_l2_offset,k_l2_offset,size_l2)
          call tce_filename('lambda2',filename)
          call createfile(filename,d_lambda2,size_l2)
          if (read_l(2)) then
            call y2_restart_read(d_lambda2,k_l2_offset,size_l2,0,
     1                           handle_l2,irrep_y)
          endif
          call gatoeaf(d_lambda2)
          if (nodezero.and.util_print('l2',print_debug)) then
            write(LuOut,*)
            write(LuOut,9000) 'l2',size_l2
            write(LuOut,9090) 'l2',d_lambda2
            call util_flush(LuOut)
          endif
        endif
!--       else  ! ----
!--        cpu = - util_cpusec()
!--        call tce_filename('t2',filename)
!--        call tce_t2_offset(l_t2_offset,k_t2_offset,size_t2)
!--        call createfile(filename,d_t2,size_t2)
!--        if(nodezero) then
!--          write(LuOut,*)'before tce_guess_t2'
!--          call util_flush(LuOut)
!--        end if
!--        call ga_zero(d_t2)
!--       if(nodezero) then
!--         write(LuOut,*)'after tce_guess_t2'
!--         call util_flush(LuOut)
!--       end if
!--        cpu = cpu + util_cpusec()
!--        if (nodezero.and.util_print('t2',print_default)) then
!--          write(LuOut,*)
!--          write(LuOut,9000) 't2',size_t2
!--          write(LuOut,9010) 't2',trim(filename)
!--c         write(LuOut,9090) 't2',d_t2
!--c         write(LuOut,9020) cpu
!--          call util_flush(LuOut)
!--        endif
!--        if (left) then
!--          call tce_y2_offset(l_l2_offset,k_l2_offset,size_l2)
!--          call tce_filename('lambda2',filename)
!--          call createfile(filename,d_lambda2,size_l2)
!--          call gatoeaf(d_lambda2)
!--        endif
!--       endif !---
        if (lineresp) then
          do axis = 1, 3
          if (respaxis(axis)) then
            irrep_x=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_x+1,irrepname)
            call tce_x2_offset(l_tr2_offset(axis),k_tr2_offset(axis),
     &           size_tr2(axis))
            call tce_filename(tr2filename(axis),filename)
            call createfile(filename,d_tr2(axis),size_tr2(axis))
            if (read_tr(2)) then
              call tr2_restart_read(d_tr2(axis),k_tr2_offset(axis),
     1          size_tr2(axis),axis,handle_tr2(axis),irrep_x)
            endif
            call gatoeaf(d_tr2(axis))
            call tce_filename(tr2filename(axis+3),filename)
            call createfile(filename,d_tr2(axis+3),
     &           size_tr2(axis))
            if (read_tr(4+2)) then
              call tr2_restart_read(d_tr2(axis+3),k_tr2_offset(axis),
     1          size_tr2(axis),axis+3,handle_tr2(axis+3),irrep_x)
            elseif (read_tr(2)) then
              call tr2_restart_read(d_tr2(axis+3),k_tr2_offset(axis),
     1          size_tr2(axis),axis,handle_tr2(axis+3),irrep_x)
            endif
            call gatoeaf(d_tr2(axis+3))
            if (nodezero.and.util_print('tr2',print_debug)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
              write(LuOut,9000) 'tr2',size_tr1(axis)
              write(LuOut,9090) 'tr2',d_tr1(axis)
              write(LuOut,9090) 'tr2',d_tr1(axis+3)
              call util_flush(LuOut)
            endif
          endif
          enddo ! axis
        endif ! lineresp
        if (leftresp) then
          do axis = 1, 3
          if (respaxis(axis)) then
            irrep_y=sym_abelian_axis(geom,axis)
            call tce_y2_offset(l_yr2_offset(axis),k_yr2_offset(axis),
     &           size_yr2(axis))
            call tce_filename(yr2filename(axis),filename)
            call createfile(filename,d_yr2(axis),size_yr2(axis))
            if (read_lr(2)) then
              call yr2_restart_read(d_yr2(axis),k_yr2_offset(axis),
     1          size_yr2(axis),axis,handle_yr2(axis),irrep_y)
            endif
            call gatoeaf(d_yr2(axis))
            call tce_filename(yr2filename(axis+3),filename)
            call createfile(filename,d_yr2(axis+3),
     &           size_yr2(axis))
            if (read_lr(2)) then
              call yr2_restart_read(d_yr2(axis+3),k_yr2_offset(axis),
     1          size_yr2(axis),axis,handle_yr2(axis+3),irrep_y)
            endif
            call gatoeaf(d_yr2(axis+3))
          endif
          enddo ! axis
        endif ! leftresp
      endif
!
!     t3 amplitudes
!
      if (needt3) then
#if !defined(CCSDT)
        CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
        cpu = - util_cpusec()
        call tce_filename('t3',filename)
        call tce_t3_offset(l_t3_offset,k_t3_offset,size_t3)
        call createfile(filename,d_t3,size_t3)
        if (read_t(3)) then
          call x3_restart_read(d_t3,k_t3_offset,size_t3,0,
     1                         handle_t3,irrep_t)
        endif
        call gatoeaf(d_t3)
        cpu = cpu + util_cpusec()
        if (nodezero.and.util_print('t3',print_default)) then
          write(LuOut,*)
          write(LuOut,9000) 't3',size_t3
          write(LuOut,9010) 't3',trim(filename)
!         write(LuOut,9090) 't3',d_t3
!         write(LuOut,9020) cpu
          call util_flush(LuOut)
        endif
        if (left) then
          irrep_y=0
          call tce_y3_offset(l_l3_offset,k_l3_offset,size_l3)
          call tce_filename('lambda3',filename)
          call createfile(filename,d_lambda3,size_l3)
          if (read_l(3)) then
            call y3_restart_read(d_lambda3,k_l3_offset,size_l3,0,
     1                           handle_l3,irrep_y)
          endif
          call gatoeaf(d_lambda3)
        endif
        if (lineresp) then
          do axis = 1, 3
          if (respaxis(axis)) then
            irrep_x=sym_abelian_axis(geom,axis)
            call tce_x3_offset(l_tr3_offset(axis),k_tr3_offset(axis),
     &           size_tr3(axis))
            call tce_filename(tr3filename(axis),filename)
            call createfile(filename,d_tr3(axis),size_tr3(axis))
            if (read_tr(3)) then
              call tr3_restart_read(d_tr3(axis),k_tr3_offset(axis),
     1          size_tr3(axis),axis,handle_tr3(axis),irrep_x)
            endif
            call gatoeaf(d_tr3(axis))
            call tce_filename(tr3filename(axis+3),filename)
            call createfile(filename,d_tr3(axis+3),
     &           size_tr3(axis))
            if (read_tr(4+3)) then
              call tr3_restart_read(d_tr3(axis+3),k_tr3_offset(axis),
     1          size_tr3(axis),axis+3,handle_tr3(axis+3),irrep_x)
            elseif (read_tr(3)) then
              call tr3_restart_read(d_tr3(axis+3),k_tr3_offset(axis),
     1          size_tr3(axis),axis,handle_tr3(axis+3),irrep_x)
            endif
            call gatoeaf(d_tr3(axis+3))
          endif
          enddo ! axis
        endif ! lineresp
        if (leftresp) then
          do axis = 1, 3
          if (respaxis(axis)) then
            irrep_y=sym_abelian_axis(geom,axis)
            call tce_y3_offset(l_yr3_offset(axis),k_yr3_offset(axis),
     &           size_yr3(axis))
            call tce_filename(yr3filename(axis),filename)
            call createfile(filename,d_yr3(axis),size_yr3(axis))
            call gatoeaf(d_yr3(axis))
            call tce_filename(yr3filename(axis+3),filename)
            call createfile(filename,d_yr3(axis+3),
     &           size_yr3(axis))
            call gatoeaf(d_yr3(axis+3))
          endif
          enddo ! axis
        endif ! leftresp
#endif
      endif
!
!     t3a amplitudes
!
      if (needt3a) then
#if !defined(CCSDT)
        CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
        cpu = - util_cpusec()
        call tce_filename('t3',filename)
        call tce_t3a_offset(l_t3_offset,k_t3_offset,size_t3)
        call createfile(filename,d_t3,size_t3)
        if (read_t(3)) then
          call x3a_restart_read(d_t3,k_t3_offset,size_t3,0,
     1                         handle_t3,irrep_t)
        endif
        call gatoeaf(d_t3)
        cpu = cpu + util_cpusec()
        if (nodezero.and.util_print('t3',print_default)) then
          write(LuOut,*)
          write(LuOut,9000) 't3a',size_t3
          write(LuOut,9010) 't3a',trim(filename)
!         write(LuOut,9090) 't3a',d_t3
!         write(LuOut,9020) cpu
          call util_flush(LuOut)
        endif
#endif
      endif
!
!     t4 amplitudes
!
      if (needt4) then
#if !defined(CCSDTQ)
      CALL ERRQUIT('CCSDTQ not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
        cpu = - util_cpusec()
        call tce_filename('t4',filename)
        call tce_t4_offset(l_t4_offset,k_t4_offset,size_t4)
        call createfile(filename,d_t4,size_t4)
        if (read_t(4)) then
          call x4_restart_read(d_t4,k_t4_offset,size_t4,0,
     1                         handle_t4,irrep_t)
        endif
        call gatoeaf(d_t4)
        cpu = cpu + util_cpusec()
        if (nodezero.and.util_print('t4',print_default)) then
          write(LuOut,*)
          write(LuOut,9000) 't4',size_t4
          write(LuOut,9010) 't4',trim(filename)
!         write(LuOut,9090) 't4',d_t4
!         write(LuOut,9020) cpu
          call util_flush(LuOut)
        endif
        if (left) then
          irrep_y=0
          call tce_y4_offset(l_l4_offset,k_l4_offset,size_l4)
          call tce_filename('lambda4',filename)
          call createfile(filename,d_lambda4,size_l4)
          if (read_l(4)) then
            call y4_restart_read(d_lambda4,k_l4_offset,size_l4,0,
     1                           handle_l4,irrep_y)
          endif
          call gatoeaf(d_lambda4)
        endif
        if (lineresp) then
          do axis = 1, 3
          if (respaxis(axis)) then
            irrep_x=sym_abelian_axis(geom,axis)
            call tce_x4_offset(l_tr4_offset(axis),k_tr4_offset(axis),
     &           size_tr4(axis))
            call tce_filename(tr4filename(axis),filename)
            call createfile(filename,d_tr4(axis),size_tr4(axis))
            if (read_tr(4)) then
              call tr4_restart_read(d_tr4(axis),k_tr4_offset(axis),
     1          size_tr4(axis),axis,handle_tr4(axis),irrep_x)
            endif
            call gatoeaf(d_tr4(axis))
            call tce_filename(tr4filename(axis+3),filename)
            call createfile(filename,d_tr4(axis+3),
     &           size_tr4(axis))
            if (read_tr(4+4)) then
              call tr4_restart_read(d_tr4(axis+3),k_tr4_offset(axis),
     1          size_tr4(axis),axis+3,handle_tr4(axis+3),irrep_x)
            elseif (read_tr(4)) then
              call tr4_restart_read(d_tr4(axis+3),k_tr4_offset(axis),
     1          size_tr4(axis),axis,handle_tr4(axis+3),irrep_x)
            endif
            call gatoeaf(d_tr4(axis+3))
          endif
          enddo ! axis
        endif ! lineresp
        if (leftresp) then
          do axis = 1, 3
          if (respaxis(axis)) then
            irrep_y=sym_abelian_axis(geom,axis)
            call tce_y4_offset(l_yr4_offset(axis),k_yr4_offset(axis),
     &           size_yr4(axis))
            call tce_filename(yr4filename(axis),filename)
            call createfile(filename,d_yr4(axis),size_yr4(axis))
            call gatoeaf(d_yr4(axis))
            call tce_filename(yr4filename(axis+3),filename)
            call createfile(filename,d_yr4(axis+3),
     &           size_yr4(axis))
            call gatoeaf(d_yr4(axis+3))
          endif
          enddo ! axis
        endif ! leftresp
#endif
      endif
!
!     =======================
!     Create residual offsets
!     =======================
!
      call tce_e_offset(l_e_offset,k_e_offset,size_e)
      if (left) then
        call tce_e_offset(l_d0_offset,k_d0_offset,size_d0)
        if (multipole.gt.1) then
          call tce_e_offset(l_q0_offset,k_q0_offset,size_q0)
        endif
        if (multipole.gt.2) then
          call tce_e_offset(l_o0_offset,k_o0_offset,size_o0)
        endif
      endif
      if (lineresp) then
        call tce_e_offset(l_a0_offset,k_a0_offset,size_a0)
      endif
      if (leftresp) then
        call tce_e_offset(l_b0_offset,k_b0_offset,size_b0)
      endif
      if (quadresp) then
        call tce_e_offset(l_g0_offset,k_g0_offset,size_g0)
      endif
      if (do_mbpt) then
       if (needt1) call tce_t1_offset(l_r1_offset,k_r1_offset,size_r1)
       if (needt2) then
        if(model.ne.'ccsd_act') then
         call tce_t2_offset(l_r2_offset,k_r2_offset,size_r2)
        else
         call tce_t2a_offset(l_r2_offset,k_r2_offset,size_r2) 
        end if
       end if
       if (needt3) call tce_t3_offset(l_r3_offset,k_r3_offset,size_r3)
       if (needt3a) call tce_t3a_offset(l_r3_offset,k_r3_offset,size_r3)
       if (needt4) call tce_t4_offset(l_r4_offset,k_r4_offset,size_r4)
      endif
!
!     ======================
!     Create auxiliary files
!     ======================
!
      if (neede2) then
        call tce_filename('e2',filename)
        call tce_e_offset(l_e2_offset,k_e2_offset,size_e2)
        call createfile(filename,d_e2,size_e2)
        call gatoeaf(d_e2)
!       write(LuOut,9090) 'e2',d_e2
      endif
      if (needu1) then
        call tce_filename('u1',filename)
        call tce_t1_offset(l_u1_offset,k_u1_offset,size_u1)
        call createfile(filename,d_u1,size_u1)
        call tce_guess_t1(d_u1,k_u1_offset)
        call reconcilefile(d_u1,size_u1)
!       write(LuOut,9090) 'u2',d_u1
      endif
      if (needu2) then
!--       if(.not.intorb) then
        call tce_filename('u2',filename)
        call tce_t2_offset(l_u2_offset,k_u2_offset,size_u2)
        call createfile(filename,d_u2,size_u2)
        call tce_guess_t2(d_v2,k_v2_offset,d_u2,k_u2_offset)
        call reconcilefile(d_u2,size_u2)
!       write(LuOut,9090) 'u2',d_u2
!--       else
!--        call tce_filename('u2',filename)
!--        call tce_t2_offset(l_u2_offset,k_u2_offset,size_u2)
!--        call createfile(filename,d_u2,size_u2)
!--c        call tce_guess_t2(d_v2,k_v2_offset,d_u2,k_u2_offset)
!--        call ga_zero(d_u2)
!--        call reconcilefile(d_u2,size_u2)
!--c       write(LuOut,9090) 'u2',d_u2
!--       endif
      endif
      if (needw1) then
        call tce_filename('w1',filename)
        call tce_t1_offset(l_w1_offset,k_w1_offset,size_w1)
        call createfile(filename,d_w1,size_w1)
        call gatoeaf(d_w1)
!       write(LuOut,9090) 'w1',d_w1
      endif
      if (needw2) then
        call tce_filename('w2',filename)
        call tce_t2_offset(l_w2_offset,k_w2_offset,size_w2)
        call createfile(filename,d_w2,size_w2)
        call gatoeaf(d_w2)
!       write(LuOut,9090) 'w2',d_w2
      endif
      if (needw3) then
        call tce_filename('w3',filename)
        call tce_t3_offset(l_w3_offset,k_w3_offset,size_w3)
        call createfile(filename,d_w3,size_w3)
        call gatoeaf(d_w3)
!       write(LuOut,9090) 'w3',d_w3
      endif
      if (needw4) then
        call tce_filename('w4',filename)
        call tce_t4_offset(l_w4_offset,k_w4_offset,size_w4)
        call createfile(filename,d_w4,size_w4)
        call gatoeaf(d_w4)
!       write(LuOut,9090) 'w4',d_w4
      endif
!
!     =========================
!     CC / CI / MBPT iterations
!     =========================
!
      if (model.eq.'eionly') then
        if (nodezero) then
          write(LuOut,9120) 'BEGIN PRINT'
          write(LuOut,9122) 'Number of alpha electrons = ',nocc(1)
          write(LuOut,9122) 'Number of beta electrons  = ',nocc(2)
          write(LuOut,9122) 'Number of alpha MOs = ',nmo(1)
          write(LuOut,9122) 'Number of beta MOs = ',nmo(2)
        endif
        if (nodezero) write(LuOut,9120) 'BEGIN 1-E'
        call tce_print_g1(rtdb,d_f1,k_f1_offset,thresh,irrep_f)
        if (nodezero) write(LuOut,9120) 'END 1-E'
        if (nodezero) write(LuOut,9120) 'BEGIN 2-E'
        call tce_print_g2(rtdb,d_v2,k_v2_offset,thresh,irrep_v)
        if (nodezero) write(LuOut,9120) 'END 2-E'
        if (left) then
          do axis = 1, 3
            irrep_d=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            if (nodezero) write(LuOut,9440) axisname(axis),irrepname
            if (nodezero) write(LuOut,9120) 'BEGIN DIPOLE'
            call tce_print_g1(rtdb,d_d1(axis),k_d1_offset(axis),thresh,
     &                        irrep_d)
            if (nodezero) write(LuOut,9120) 'END DIPOLE'
          enddo ! axis
        endif ! left
        if (nodezero) write(LuOut,9120) 'END PRINT'
        if (nodezero) call util_flush(LuOut)
      endif
      if ((model.eq.'ccd').or.do_ccd) then
! -------------
! CCD     right
! -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CCD"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call ccd_e(d_e,d_t2,d_v2,
     1               k_e_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call ccd_t2(d_f1,d_r2,d_t2,d_v2,
     1                k_f1_offset,k_t2_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_r2,size_t2)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          residual = r2
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "CCD",corr
              write(LuOut,9080) "CCD",ref + corr
            endif
            call deletefile(d_r2)
            call deletefile(d_e)
            if (ampnorms) then
              call tce_residual_t2(d_t2,k_t2_offset,r2)
              if (nodezero) then
                write(LuOut,9082) "T doubles",r2
              endif
            endif
            call tce_print_x2(d_t2,k_t2_offset,printtol,irrep_t)
            call tce_diis_tidy()
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
            goto 101
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          call tce_diis(.false.,iter,.false.,.true.,.false.,.false.,
     1                  dummy,dummy,dummy,dummy,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r2)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
  101 continue
      endif
      if ((model.eq.'lccd').or.do_lccd) then
! -------------
! LCCD    right
! -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "LCCD"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call lccd_e(d_e,d_t2,d_v2,
     1                k_e_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call lccd_t2(d_f1,d_r2,d_t2,d_v2,
     1                 k_f1_offset,k_t2_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_r2,size_t2)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          residual = r2
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "LCCD",corr
              write(LuOut,9080) "LCCD",ref + corr
            endif
            call deletefile(d_r2)
            call deletefile(d_e)
            if (ampnorms) then
              call tce_residual_t2(d_t2,k_t2_offset,r2)
              if (nodezero) then
                write(LuOut,9082) "T doubles",r2
              endif
            endif
            call tce_print_x2(d_t2,k_t2_offset,printtol,irrep_t)
            call tce_diis_tidy()
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
            goto 102
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          call tce_diis(.false.,iter,.false.,.true.,.false.,.false.,
     1                  dummy,dummy,dummy,dummy,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r2)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
  102 continue
      endif
!
!
       if (model.eq.'ccsd_act') then
! ------------------------
!  ccsd_act/eomccsd_act
! ------------------------
! === jaguar ===
      if(read_ta)  then
       call tce_zero(d_t1,size_t1)
       call tce_zero(d_t2,size_t2)
       call util_file_name0('t1amp',.false.,.true.,filename,fldgts)
       unitn=79
       call read_tensor(filename,d_t1,size_t1,unitn)
       call util_file_name0('t2amp',.false.,.true.,filename,fldgts)
       unitn=80
       call read_tensor(filename,d_t2,size_t2,unitn)
       call ga_sync()
       go to 3203
      end if
! ==============
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CCSD-act"
! - T1/X1 LOCALIZATION ----------
       if(ccsd_var.ne.'cc2') then
!    opening l_t1_local
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
!    copy d_t1 ==> l_t1_local
        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
       end if
! -------------------------------
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
! - T1/X1 LOCALIZATION --
          call ccsd_ea(d_f1,d_e,k_t1_local,d_t2,d_v2,
     1                k_f1_offset,k_e_offset,
     2                k_t1_offset,k_t2_offset,k_v2_offset)
! -----------------------
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          size_r1=size_t1
          call createfile(filename,d_r1,size_r1)
! - T1/X1 LOCALIZATION --
          call ccsd_t1a(d_f1,d_r1,k_t1_local,d_t2,d_v2,
     1                 k_f1_offset,k_t1_offset,
     2                 k_t1_offset,k_t2_offset,k_v2_offset)
! -----------------------
          call reconcilefile(d_r1,size_r1)
          size_r2=size_t2
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_r2)
! - T1/X1 LOCALIZATION --
          call ccsd_t2a(d_f1,d_r2,k_t1_local,d_t2,d_v2,
     1                 k_f1_offset,k_t2_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset)!,size_t2)
! -----------------------
          call reconcilefile(d_r2,size_r2)
          call tce_residual_t1(d_r1,k_t1_offset,r1)
          call tce_residual_t2a(d_r2,k_t2_offset,r2)
          residual = max(r1,r2)
          call get_block(d_e,corr,1,0)
! - T1/X1 LOCALIZATION --
       if(ccsd_var.ne.'cc2') then
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',2,MA_ERR)
       end if
! -----------------------
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
                write(LuOut,9070) "CCSD",corr
                write(LuOut,9080) "CCSD",ref + corr
            endif
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            irrep_x = irrep_t ! Caution
            call tce_print_x1(d_t1,k_t1_offset,printtol,irrep_t)
            call tce_print_x2a(d_t2,k_t2_offset,printtol,irrep_t)
            call tce_diis_tidy()
! === jaguar ===
      if(write_ta)  then
       call util_file_name0('t1amp',.false.,.true.,filename,fldgts)
       unitn=79
       call write_tensor(filename,d_t1,size_t1,unitn)
       call util_file_name0('t2amp',.false.,.true.,filename,fldgts)
        unitn=80
       call write_tensor(filename,d_t2,size_t2,unitn)
       call ga_sync()
      end if
! ==============
            goto 3203
          endif
          call tce_diis_act(.false.,iter,.true.,.true.,.false.,.false.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 3203   continue
!
        end if ! ccsd_act/eomccsd_act
!
      if ((model.eq.'ccsd').or.do_ccsd) then
! -------------
! CCSD    right
! -------------
! === jaguar ===
          if(read_ta)  then
           call tce_zero(d_t1,size_t1)
           call tce_zero(d_t2,size_t2)
           call util_file_name0('t1amp',.false.,.true.,filename,fldgts)
           unitn=79
           call read_tensor(filename,d_t1,size_t1,unitn)
           call util_file_name0('t2amp',.false.,.true.,filename,fldgts)
           unitn=80
           call read_tensor(filename,d_t2,size_t2,unitn)
           call ga_sync()
!cx           go to 320
          end if
! ==============
#ifdef MIRROR_SING
          call ccsd_energy_noloc(d_e,d_f1,d_v2,d_t1,d_t2,
     1         k_e_offset,k_f1_offset,k_v2_offset,
     2         k_t1_offset,k_t2_offset,
     3         size_t1,size_t2,ref,corr)
#else
          call ccsd_energy_loc(d_e,d_f1,d_v2,d_t1,d_t2,
     1         k_e_offset,k_f1_offset,k_v2_offset,
     2         k_t1_offset,k_t2_offset,
     3         size_t1,size_t2,ref,corr,rtdb)
!QuantumComputing
          if(quantumc) then
          if(nodezero) then
           write(6,*)'Ground state specification:'
           write(6,*)'exp(T)|Ref>~(1+T)|Ref>'
           call util_flush(6)
          endif
          if(nodezero) write(6,*)'------- a(p)+/a(q) strings ----------'
           xnorm1 = 0.0d0
           xnorm2 = 0.0d0
           xnorm  = 0.0d0
c
           call tce_string_x1_norm(d_t1,k_t1_offset,0.10d0,irrep_t,
     &                             xnorm1)
           call tce_string_x2_norm(d_t2,k_t2_offset,0.10d0,irrep_t,
     &                             xnorm2)
c normalization only on nodezero / reference included
           if(nodezero)  then
            xnorm=(1.0d0)/(dsqrt(1.0d0+xnorm1+xnorm2))
           endif
c
          call tce_string_ref(xnorm)
          call tce_string_x1(d_t1,k_t1_offset,0.10d0,irrep_t,xnorm)
          call tce_string_x2(d_t2,k_t2_offset,0.10d0,irrep_t,xnorm)
          if(nodezero) write(6,*)'-------------------------------------'
c
          endif !quantumc
c
c qducc starts here
c
          if(qducc) then
             if (.not.(hftype.eq.'RHF'))
     1        call errquit('tce_energy:qducc only with RHF reference',0,
     2        RTDB_ERR)
           nos=int_mb(k_offset+noab-1)+int_mb(k_range+noab-1)
           nvs=int_mb(k_offset+noab+nvab-1)+int_mb(k_range+noab+nvab-1)
     &         -nos
           noas=nos/2
           nobs=nos/2
           nvas=nvs/2
           nvbs=nvs/2
c -------------------
c DUCC SECTION
c
           if (.not. geom_nuc_rep_energy(geom, enrep_tce))
     &      call errquit('TCE:unable to get nuclear rep ',geom, 
     &      GEOM_ERR)
c
          if(nodezero) then
           write(6,*)'Ground state specification:'
           write(6,*)'exp(T)|Ref>~(1+T)|Ref>'
           call util_flush(6)
          endif
          if(nodezero) write(6,*)'------- a(p)+/a(q) strings ----------'
           xnorm1 = 0.0d0
           xnorm2 = 0.0d0
           xnorm  = 0.0d0
c
           call tce_string_x1_norm(d_t1,k_t1_offset,0.10d0,irrep_t,
     &                             xnorm1)
           call tce_string_x2_norm(d_t2,k_t2_offset,0.10d0,irrep_t,
     &                             xnorm2)
c normalization only on nodezero / reference included
           if(nodezero)  then
            xnorm=(1.0d0)/(dsqrt(1.0d0+xnorm1+xnorm2))
           endif
c
          call tce_string_ref(xnorm)
          call tce_string_x1(d_t1,k_t1_offset,0.10d0,irrep_t,xnorm)
          call tce_string_x2(d_t2,k_t2_offset,0.10d0,irrep_t,xnorm)
          if(nodezero) write(6,*)'-------------------------------------'
c
      if(nodezero) then
c yaml enerep_tce - coulomb_repulsion
       write(6,*)'enrep_tce = ',enrep_tce
c yaml red - this is hf_energy
       write(6,6888) ref
c yaml scf_energy_offset
       write(6,6889) 0.0d0
       write(6,*)'ALL OCCUPIED ORBITALS ARE ACTIVE IN DUCC --> NO SHIFT'
       write(6,*)'Number of active orbitals',nactv+(nos/2)
       write(6,*)'Number of active alpha electrons',noas
       write(6,*)'Number of active beta electrons',nobs
c
      call util_flush(6)
      endif 

         if(nodezero) then
          write(6,*)' nos ',nos
          write(6,*)' nvs ',nvs
          write(6,*)'noas ',noas
          write(6,*)'nobs ',nobs
          write(6,*)'nvas ',nvas
          write(6,*)'nvbs ',nvbs
          write(6,*)'nactv ',nactv
          call util_flush(6)
         endif
c -------------------
        if (.not.MA_PUSH_GET(mt_int,nos+nvs,'transform',
     1      l_transf,k_transf))
     1      call errquit('transform',100,MA_ERR)
        if (.not.MA_PUSH_GET(mt_int,(nos+nvs)**2,'act1',
     1      l_act1,k_act1))
     1      call errquit('act1',101,MA_ERR)
        if (.not.MA_PUSH_GET(mt_int,(nos+nvs)**2,'prt1',
     1      l_prt1,k_prt1))
     1      call errquit('transform',102,MA_ERR)
        if (.not.MA_PUSH_GET(mt_byte,(nos+nvs)**4,'act2',
     1      l_act2,k_act2))
     1      call errquit('act2',103,MA_ERR)
        if (.not.MA_PUSH_GET(mt_byte,(nos+nvs)**4,'prt2',
     1      l_prt2,k_prt2))
     1      call errquit('prt2',104,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos+nvs)**2,'h',
     1      l_h,k_h))
     1      call errquit('h',200,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,((nos+nvs)/2)**2,'horb',
     1      l_horb,k_horb))
     1      call errquit('horb',201,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos+nvs)**4,'v',
     1      l_v,k_v))
     1      call errquit('v',202,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,nos*(nos+nvs),'t1',
     1      l_t1,k_t1))
     1      call errquit('t1',203,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos*nvs)**2,'t2',
     1      l_t2,k_t2))
     1      call errquit('t2',204,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos+nvs)**2,'ht',
     1      l_ht,k_ht))
     1      call errquit('ht',205,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos+nvs)**4,'vt',
     1      l_vt,k_vt))
     1      call errquit('vt',206,MA_ERR)
           call ducc(rtdb,d_t1,k_t1_offset,d_t2,k_t2_offset,
     &    d_v2,k_v2_offset,dbl_mb(k_evl_sorted),
     &    nos,nvs,noas,nobs,nvas,nvbs,nactv,
     D        int_mb(k_transf),int_mb(k_act1),int_mb(k_prt1),
     D        byte_mb(k_act2),byte_mb(k_prt2),
     D        dbl_mb(k_h),dbl_mb(k_horb),dbl_mb(k_v),
     D        dbl_mb(k_t1),dbl_mb(k_t2),
     D        dbl_mb(k_ht),dbl_mb(k_vt))
              if (.not.ma_chop_stack(l_transf))
     1        call errquit('tce_energy: chopstack ducc',4,MA_ERR)
          endif ! qducc
c
c
!QuantumComputing
#endif
! === jaguar ===
         if(write_ta)  then
          call util_file_name0('t1amp',.false.,.true.,filename,fldgts)
          unitn=79
          call write_tensor(filename,d_t1,size_t1,unitn)
          call util_file_name0('t2amp',.false.,.true.,filename,fldgts)
          unitn=80
          call write_tensor(filename,d_t2,size_t2,unitn)
          call ga_sync()
         end if
! ==============
  320    continue
! -------------
! CCSD    left
! -------------
        if (left) then
          call ccsd_lambda(
     1           d_f1,d_v2,d_d1,d_q1,d_o1,
     2           d_t1,d_t2,d_lambda1,d_lambda2,
     3           k_f1_offset,k_v2_offset,
     4           k_d1_offset,k_q1_offset,k_o1_offset,
     5           k_t1_offset,k_t2_offset,
     6           k_l1_offset,k_l2_offset,
     7           size_l1,size_l2,
     8           dipole_scf,    dipole_cor,
     9           quadrupole_scf,quadrupole_cor,
     1           octupole_scf,  octupole_cor)
!
!->d3p975
!
! calculated ccsd 1-particle reduced density matrix/ces
!
            if (nroots.eq.0.and.idens) then
               if (.not.intorb)
     1            call ccsd_1prdm(rtdb,d_t1,d_t2,d_lambda1,d_lambda2,
     1                            k_t1_offset,k_t2_offset,
     2                            k_l1_offset,k_l2_offset)
            endif
!
!<-d3p975
!
        endif ! left
 2237   continue
! -------------
! CCSD-LR
! -------------
        if (lineresp) then
#if !defined(CCSDLR)
        call errquit('CCSDLR not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
          ! no need to pass d_a0 as argument here
          call ccsd_lr_alpha(rtdb,d_d0,d_a0,d_f1,d_v2,d_d1,
     1         d_t1,d_t2,d_lambda1,d_lambda2,d_tr1,d_tr2,
     2         k_d0_offset,k_a0_offset,
     3         k_f1_offset,k_v2_offset,k_d1_offset,
     4         k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset,
     5         k_tr1_offset,k_tr2_offset,
     6         size_tr1,size_tr2,dipole_scf,dipole_cor)
#endif
        endif ! lineresp
! -------------
! CCSD-YR
! -------------
        if (leftresp) then
#if !defined(CCSDLR)
        call errquit('CCSDLR not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
          ! no need to pass d_a0 as argument here
          call ccsd_qr_beta(rtdb,d_a0,d_b0,d_f1,d_v2,d_d1,
     1         d_t1,d_t2,d_lambda1,d_lambda2,d_yr1,d_yr2,d_tr1,d_tr2,
     2         k_a0_offset,k_b0_offset,
     3         k_f1_offset,k_v2_offset,k_d1_offset,
     4         k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset,
     5         k_tr1_offset,k_tr2_offset,k_yr1_offset,k_yr2_offset,
     6         size_yr1,size_yr2)
#endif
        endif ! leftresp
      endif ! (model.eq.'ccsd').or.do_ccsd
! -------------------
! CCSD(T) and CCSD[T]
! -------------------
        if (do_pt.or.((perturbative.eq.'(t)').or.
     1                (perturbative.eq.'[t]'))) then
!
          if (nodezero) write(LuOut,*) 'CCSD(T)'
          if (do_pt) perturbative='(t)'
!
          cpu = - util_cpusec()
          wall = - util_wallsec()
#ifdef USE_PSTAT
c
c     profiling
c
          call ccsd_t_pstat_init(rtdb)
          call ccsd_t_fdist_init()
#endif
          if (restart_ccsd_t.eq.0) then
#if defined(TCE_CUDA) || defined(TCE_HIP)
              if(.not.rtdb_get(rtdb,'tce:cuda',mt_int,1,icuda)) then
               if (nodezero) write(LuOut,*) 'Using plain CCSD(T) code'
               if (nodezero) call util_flush(LuOut)
               call ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                    d_v2,k_v2_offset,pt3_1,pt3_2,size_t1)
              else
               if (nodezero) write(LuOut,*) 'Using CUDA CCSD(T) code'
!kbn           call errquit('tce_energy: reading cuda',icuda,CALC_ERR) 
               if (nodezero) call util_flush(LuOut)
               call ccsd_t_gpu(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                    d_v2,k_v2_offset,pt3_1,pt3_2,size_t1,icuda)
              end if 
#else
              if (nodezero) write(LuOut,*) 'Using plain CCSD(T) code'
              if (nodezero) call util_flush(LuOut)
              if(.not.slice) then
              call ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                    d_v2,k_v2_offset,pt3_1,pt3_2,size_t1)
              else
              if (nodezero) write(LuOut,*) 'Using sliced CCSD(T) code'
              call ccsd_t_6dts(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                    d_v2,k_v2_offset,pt3_1,pt3_2,size_t1,xmem)
              end if
#endif
#ifdef USE_PSTAT
c     
c     cleanup profiling
c
              call ccsd_t_pstat_print
#endif

          elseif (restart_ccsd_t.gt.0) then
            pt3_1 = 0.0d0
            call ccsd_t_restart(rtdb,d_t1,k_t1_offset,
     1                               d_t2,k_t2_offset,
     2                               d_v2,k_v2_offset,pt3_2)
          endif
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
          if (nodezero) then
            write(LuOut,*)
            if (restart_ccsd_t.eq.0) then
              write(LuOut,9075) "CCSD[T]",pt3_1
              write(LuOut,9070) "CCSD[T]",corr+pt3_1
              write(LuOut,9080) "CCSD[T]",ref+corr+pt3_1
            endif
            write(LuOut,9075) "CCSD(T)",pt3_2
            write(LuOut,9070) "CCSD(T)",corr+pt3_2
            write(LuOut,9080) "CCSD(T)",ref+corr+pt3_2
            write(LuOut,9020) cpu, wall
          endif
          if (perturbative.eq.'[t]') then
            corr = corr + pt3_1
          elseif (perturbative.eq.'(t)') then
            corr = corr + pt3_2
          endif
        endif
! ---------------------------------
! Lambda-CCSD(T) and Lambda-CCSD[T]
! ---------------------------------
        if (do_lam_pt.or.((perturbative.eq.'lambda(t)').or.
     &                    (perturbative.eq.'lambda[t]'))) then
!
          if (nodezero) write(LuOut,*) 'Lambda-CCSD(T)'
          if (do_lam_pt) perturbative='lambda(t)'
!
          cpu = - util_cpusec()
          wall = - util_wallsec()
          call lambda_ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                d_lambda1,k_l1_offset,d_lambda2,k_l2_offset,
     2                d_f1,k_f1_offset,d_v2,k_v2_offset,pt3_1,pt3_2)
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9075) "Lambda-CCSD[T]",pt3_1
            write(LuOut,9070) "Lambda-CCSD[T]",corr+pt3_1
            write(LuOut,9080) "Lambda-CCSD[T]",ref+corr+pt3_1
            write(LuOut,9075) "Lambda-CCSD(T)",pt3_2
            write(LuOut,9070) "Lambda-CCSD(T)",corr+pt3_2
            write(LuOut,9080) "Lambda-CCSD(T)",ref+corr+pt3_2
            write(LuOut,9020) cpu, wall
          endif
          if (perturbative.eq.'lambda[t]') then
            corr = corr + pt3_1
          elseif (perturbative.eq.'lambda(t)') then
            corr = corr + pt3_2
          endif
        endif
! -------------------------
! CR-CCSD(T) and CR-CCSD[T]
! -------------------------
        if(do_cr_pt.or.((perturbative.eq.'cr_(t)')
     1             .or. (perturbative.eq.'cr_[t]'))) then
!
          if (nodezero) write(LuOut,*) 'CR-CCSD(T)'
          if (do_cr_pt) perturbative='cr_(t)'
!
          if ((model.eq.'ccsd').or.(model.eq.'multi')) then
            cpu = - util_cpusec()
            wall = - util_wallsec()
            call tce_filename('e',filename)
            call createfile(filename,d_e,size_e)
            call cr_ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                     d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                     d_e,k_e_offset,pt3_1,pt3_2,size_t1)
            call deletefile(d_e)
            cpu = cpu + util_cpusec()
            wall = wall + util_wallsec()
            if (nodezero) then
              write(LuOut,*)
              write(LuOut,9075) "CR-CCSD[T]",pt3_1
              write(LuOut,9070) "CR-CCSD[T]",corr+pt3_1
              write(LuOut,9080) "CR-CCSD[T]",ref+corr+pt3_1
              write(LuOut,9075) "CR-CCSD(T)",pt3_2
              write(LuOut,9070) "CR-CCSD(T)",corr+pt3_2
              write(LuOut,9080) "CR-CCSD(T)",ref+corr+pt3_2
              write(LuOut,9020) cpu, wall
            endif
            if (perturbative.eq.'cr_[t]') then
              corr = corr + pt3_1
            elseif (perturbative.eq.'cr_(t)') then
              corr = corr + pt3_2
            endif
          endif
        endif
! ----------------------------------------------------
! LR-CCSD(T) and/or LR-CCSD(T)-1 and/or LR-CCSD(TQ)-1p
! ----------------------------------------------------
        if (do_lcr_pt.or.((perturbative.eq.'lr_(t)').or.
     &                    (perturbative.eq.'lr_(tq1p)').or.
     &                    (perturbative.eq.'lr_(tq1)'))) then
!
          if (nodezero) write(LuOut,*) 'LR-CCSD(T)'
          if (do_lcr_pt) perturbative='lr_(t)'
!
          cpu = - util_cpusec()
          wall = - util_wallsec()
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call lr_ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                   d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                   d_e,k_e_offset,
     3                   pt3_1,pt3_2,pt3_3,pt3_4,pt3_5,pt3_6,
     4                   size_t1)
          call deletefile(d_e)
!
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9070) "LR-CCSD(T)IA",corr+pt3_1
            write(LuOut,9080) "LR-CCSD(T)IA",ref+corr+pt3_1
            write(LuOut,*)
            write(LuOut,9070) "LR-CCSD(T)IB",corr+pt3_2
            write(LuOut,9080) "LR-CCSD(T)IB",ref+corr+pt3_2
            write(LuOut,*)
            write(LuOut,9070) "LR-CCSD(T)IIA",corr+pt3_3
            write(LuOut,9080) "LR-CCSD(T)IIA",ref+corr+pt3_3
            write(LuOut,*)
            write(LuOut,9070) "LR-CCSD(T)IIB",corr+pt3_4
            write(LuOut,9080) "LR-CCSD(T)IIB",ref+corr+pt3_4
            write(LuOut,*)
            write(LuOut,9070) "LR-CCSD(T)IIIA",corr+pt3_5
            write(LuOut,9080) "LR-CCSD(T)IIIA",ref+corr+pt3_5
            write(LuOut,*)
            write(LuOut,9070) "LR-CCSD(T)IIIB",corr+pt3_6
            write(LuOut,9080) "LR-CCSD(T)IIIB",ref+corr+pt3_6
            write(LuOut,*)
            write(LuOut,9020) cpu, wall
          endif
          if (perturbative.eq.'lr_(t)') then
            corr = corr + pt3_2
          endif
! --------------------------------------
! LR-CCSD(TQ)-1p
! --------------------------------------
          if(perturbative.eq.'lr_(tq1p)') then
            cpu = - util_cpusec()
            wall = - util_wallsec()
            call tce_filename('e',filename)
            call createfile(filename,d_e,size_e)
            call  lr_ccsd_tq1p(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                         d_v2,k_v2_offset,
     2                         d_e,k_e_offset,pt4_lr)
            call deletefile(d_e)
            cpu = cpu + util_cpusec()
            wall = wall + util_wallsec()
            if (nodezero) then
              write(LuOut,*)
              write(LuOut,9070) "LR-CCSD(TQ)-1p,B",corr+pt3_2+pt4_lr
              write(LuOut,9080) "LR-CCSD(TQ)-1p,B",ref+corr+pt3_2+pt4_lr
              write(LuOut,*)
              write(LuOut,9020) cpu, wall
            end if
            if (perturbative.eq.'lr_(tq1p)') then
              corr = corr + pt3_2 + pt4_lr
            end if
          end if
! --------------------------------------------------------
! LR-CCSD(TQ)-1
! --------------------------------------------------------
          if(perturbative.eq.'lr_(tq1)') then
            cpu = - util_cpusec()
            wall = - util_wallsec()
            call tce_filename('e',filename)
            call createfile(filename,d_e,size_e)
            call  lr_ccsd_tq_1(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                         d_v2,k_v2_offset,
     2                         d_e,k_e_offset,pt4_lr)
            call deletefile(d_e)
            cpu = cpu + util_cpusec()
            wall = wall + util_wallsec()
            if (nodezero) then
              write(LuOut,*)
              write(LuOut,9070) "LR-CCSD(TQ)-1,B",corr+pt3_2+pt4_lr
              write(LuOut,9080) "LR-CCSD(TQ)-1,B",ref+corr+pt3_2+pt4_lr
              write(LuOut,*)
              write(LuOut,9020) cpu, wall
            end if
            if (perturbative.eq.'lr_(tq1)') then
              corr = corr + pt3_2 + pt4_lr
            end if
          end if
        end if
! ---------------------
! CCSD(2)_T and CCSD(2)
! ---------------------
        if (do_2t_pt.or.((perturbative.eq.'2_t').or.
     1                   (perturbative.eq.'2_tq').or.
     2                   (perturbative.eq.'2_tqf'))) then
!
          if (nodezero) write(LuOut,*) 'CCSD(2)_T'
          if (do_2t_pt) perturbative='2_t'
!
          cpu = - util_cpusec()
          wall = - util_wallsec()
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call ccsd2_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                 d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                 d_lambda1,k_l1_offset,d_lambda2,k_l2_offset,
     3                 d_e,k_e_offset,pt3_1,size_t1)
          if (perturbative.eq.'2_tq') then
            call ccsd2_q(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                   d_v2,k_v2_offset,d_lambda2,k_l2_offset,
     3                   d_e,k_e_offset,pt3_2)
          endif
          call deletefile(d_e)
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
          if (nodezero) then
            write(LuOut,*)
            if (perturbative.eq.'2_t') then
              write(LuOut,9070) "CCSD(2)_T",corr+pt3_1
              write(LuOut,9080) "CCSD(2)_T",ref+corr+pt3_1
            elseif (perturbative.eq.'2_tq') then
              write(LuOut,9070) "CCSD(2)",corr+pt3_1+pt3_2
              write(LuOut,9080) "CCSD(2)",ref+corr+pt3_1+pt3_2
            endif
            write(LuOut,9020) cpu, wall
          endif
          corr = corr + pt3_1
        endif
!cx       endif
      !JEFFendif
      if ((model.eq.'lccsd').or.do_lccsd) then
! -------------
! LCCSD   right
! -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "LCCSD"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call lccsd_e(d_f1,d_e,d_t1,d_t2,d_v2,
     1                 k_f1_offset,k_e_offset,
     2                 k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_t1)
          call lccsd_t1(d_f1,d_r1,d_t1,d_t2,d_v2,
     1                  k_f1_offset,k_t1_offset,
     2                  k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_r1,size_t1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call lccsd_t2(d_f1,d_r2,d_t1,d_t2,d_v2,
     1                  k_f1_offset,k_t2_offset,
     2                  k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_r2,size_t2)
          call tce_residual_t1(d_r1,k_t1_offset,r1)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          residual = max(r1,r2)
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "LCCSD",corr
              write(LuOut,9080) "LCCSD",ref + corr
            endif
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            if (ampnorms) then
              call tce_residual_t1(d_t1,k_t1_offset,r1)
              call tce_residual_t2(d_t2,k_t2_offset,r2)
              if (nodezero) then
                write(LuOut,9082) "T singles",r1
                write(LuOut,9082) "T doubles",r2
              endif
            endif
            call tce_print_x1(d_t1,k_t1_offset,printtol,irrep_t)
            call tce_print_x2(d_t2,k_t2_offset,printtol,irrep_t)
            call tce_diis_tidy()
            if (save_t(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1 now...'
              endif
              call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                             handle_t1,irrep_t)
            endif
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
            goto 104
          endif
          if (save_t(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1 now...'
            endif
            call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                           handle_t1,irrep_t)
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          call tce_diis(.false.,iter,.true.,.true.,.false.,.false.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
  104 continue
      endif
      if (model.eq.'ccsdt') then
#if !defined(CCSDT)
        CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
! -------------
! CCSDT   right
! -------------
        call tce_diis_init()
c  Tsource
      if(intorb) then
         d_v2copy=d_v2
         k_v2copy_offset=k_v2_offset
         d_v2=d_v2aux
         k_v2_offset=k_v2aux_offset
      endif !intorb
c end of Tsource
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if(.not.ucc4) then
           if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CCSDT"
          else !ucc4
           if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "UCC4"
          endif
c
         if(ucc4) then
          if(nodezero) then
           write(6,*)'ucc4: iteration initlzation'
            call util_flush(6)
          endif
          irrep_y=0
          call tce_y2_offset(l_l2_offset,k_l2_offset,size_l2)
          call tce_filename('lambda2',filename)
          call createfile(filename,d_lambda2,size_l2)
c call transition procedure
          call t2_2_y2(d_t2,d_lambda2,k_t2_offset,k_l2_offset)
         endif
c
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          if(.not.ucc4) then
          call ccsdt_e(d_f1,d_e,d_t1,d_t2,d_v2,
     1                 k_f1_offset,k_e_offset,
     2                 k_t1_offset,k_t2_offset,k_v2_offset)
          else ! ucc4
           if(nodezero) then
            write(6,*)'UCC4 energy'
            call util_flush(6)
          endif
           call ucc4_e_new(d_e,d_t2,d_v2,d_lambda2,k_e_offset,
     &          k_t2_offset,
     &          k_v2_offset,k_l2_offset)
          endif
          call reconcilefile(d_e,size_e)
c
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_t1)
          if(.not.ucc4) then
          call ccsdt_t1(d_f1,d_r1,d_t1,d_t2,d_t3,d_v2,
     1                  k_f1_offset,k_t1_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset)
          else  ! ucc4
          if(nodezero) then
            write(6,*)'UCC4 singles'
            call util_flush(6)
          endif
       call ucc4_t1_new(d_f1,d_r1,d_t1,d_t2,d_v2,k_f1_offset,
     &         k_t1_offset,
     &         k_t1_offset,k_t2_offset,k_v2_offset)
          endif
          call reconcilefile(d_r1,size_t1)
c
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          if(.not.ucc4) then
          call ccsdt_t2(d_f1,d_r2,d_t1,d_t2,d_t3,d_v2,
     1                  k_f1_offset,k_t2_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset,size_t2)
          else !ucc4
          if(nodezero) then
           write(6,*)'UCC4 doubles'
           call util_flush(6)
          endif
           call ucc4_t2(d_f1,d_r2,d_t1,d_t2,d_t3,d_v2,k_f1_offset,
     &           k_t2_offset,
     &           k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset)
           call ucc4_t2_t2_v_t(d_r2,d_t2,d_v2,d_lambda2,k_t2_offset,
     &           k_t2_offset,k_v2_offset,k_l2_offset)
          endif
          call reconcilefile(d_r2,size_t2)
c
          call tce_filename('r3',filename)
          call createfile(filename,d_r3,size_t3)
          if(.not.ucc4) then
          call ccsdt_t3(d_f1,d_r3,d_t1,d_t2,d_t3,d_v2,
     1                  k_f1_offset,k_t3_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset)
          else  ! ucc4
          if(nodezero) then
           write(6,*)'UCC4 triples'
           call util_flush(6)
          endif
           call ucc4_t3(d_f1,d_r3,d_t2,d_t3,d_v2,k_f1_offset,
     &       k_t3_offset,k_t2_offset,k_t3_offset,k_v2_offset)
          endif
          call reconcilefile(d_r3,size_t3)
          call tce_residual_t1(d_r1,k_t1_offset,r1)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          call tce_residual_t3(d_r3,k_t3_offset,r3)
          residual = max(r1,r2,r3)
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
             if(.not.ucc4) then
              write(LuOut,9060)
              write(LuOut,9070) "CCSDT",corr
              write(LuOut,9080) "CCSDT",ref + corr
             else 
              write(LuOut,9060)
              write(LuOut,9070) "UCC(4)",corr
              write(LuOut,9080) "UCC(4)",ref + corr
             endif 
            endif
            if (diagnt1) call tce_diagnose_t1(d_t1,k_t1_offset,t1diagn)
            if (ampnorms) then
              call tce_residual_t1(d_t1,k_t1_offset,r1)
              call tce_residual_t2(d_t2,k_t2_offset,r2)
              call tce_residual_t3(d_t3,k_t3_offset,r3)
              if (nodezero) then
                write(LuOut,9082) "T singles",r1
                write(LuOut,9082) "T doubles",r2
                write(LuOut,9082) "T triples",r3
              endif
            endif
            call deletefile(d_r3)
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
c
c delete d_lambda related objects: ucc4
c
            if(ucc4) then
              call deletefile(d_lambda2)
              if (.not.ma_pop_stack(l_l2_offset))
     1        call errquit("tce_energy: l2_offset ucc4",4,MA_ERR)
            endif
c
            call tce_print_x1(d_t1,k_t1_offset,printtol,irrep_t)
            call tce_print_x2(d_t2,k_t2_offset,printtol,irrep_t)
c
c
c
c Tsource
          if(qducc) then
           if(nodezero) then
            if(ucc4) then
             write(LuOut,*)' UCC(4)  '
            else
             write(LuOut,*)' CCSDT  '
            endif
           endif
             if (.not.(hftype.eq.'RHF'))
     1        call errquit('tce_energy:qducc only with RHF reference',0,
     2        RTDB_ERR)
           nos=int_mb(k_offset+noab-1)+int_mb(k_range+noab-1)
           nvs=int_mb(k_offset+noab+nvab-1)+int_mb(k_range+noab+nvab-1)
     &         -nos
           noas=nos/2
           nobs=nos/2
           nvas=nvs/2
           nvbs=nvs/2

CNPB           if (.not. geom_nuc_rep_energy(geom, enrep_tce))
CNPB     &      call errquit('TCE:unable to get nuclear rep ',geom, 
CNPB     &      GEOM_ERR)
CNPBc
CNPB          if(nodezero) then
CNPB           write(6,*)'Ground state specification:'
CNPB           write(6,*)'exp(T)|Ref>~(1+T)|Ref>'
CNPB           call util_flush(6)
CNPB          endif
CNPB          if(nodezero) write(6,*)'------- a(p)+/a(q) strings ----------'
CNPB           xnorm1 = 0.0d0
CNPB           xnorm2 = 0.0d0
CNPB           xnorm  = 0.0d0
CNPBc
CNPB           call tce_string_x1_norm(d_t1,k_t1_offset,0.10d0,irrep_t,
CNPB     &                             xnorm1)
CNPB           call tce_string_x2_norm(d_t2,k_t2_offset,0.10d0,irrep_t,
CNPB     &                             xnorm2)
CNPBc normalization only on nodezero / reference included
CNPB           if(nodezero)  then
CNPB            xnorm=(1.0d0)/(dsqrt(1.0d0+xnorm1+xnorm2))
CNPB           endif
CNPBc
CNPB          call tce_string_ref(xnorm)
CNPB          call tce_string_x1(d_t1,k_t1_offset,0.10d0,irrep_t,xnorm)
CNPB          call tce_string_x2(d_t2,k_t2_offset,0.10d0,irrep_t,xnorm)
CNPB          if(nodezero) write(6,*)'-------------------------------------'
CNPBc
      if(nodezero) then
c yaml enerep_tce - coulomb_repulsion
       write(6,*)'enrep_tce = ',enrep_tce
c yaml red - this is hf_energy
       write(6,6888) ref
c yaml scf_energy_offset
       write(6,6889) 0.0d0
       write(6,*)'ALL OCCUPIED ORBITALS ARE ACTIVE IN DUCC --> NO SHIFT'
       write(6,*)'Number of active orbitals',nactv+(nos/2)
       write(6,*)'Number of active alpha electrons',noas
       write(6,*)'Number of active beta electrons',nobs
c
      call util_flush(6)
      endif 
c -------------------
         if(nodezero) then
          write(6,*)' nos ',nos
          write(6,*)' nvs ',nvs
          write(6,*)'noas ',noas
          write(6,*)'nobs ',nobs
          write(6,*)'nvas ',nvas
          write(6,*)'nvbs ',nvbs
          write(6,*)'nactv ',nactv
          call util_flush(6)
         endif
c -------------------
        if (.not.MA_PUSH_GET(mt_int,nos+nvs,'transform',
     1      l_transf,k_transf))
     1      call errquit('transform',100,MA_ERR)
        if (.not.MA_PUSH_GET(mt_int,(nos+nvs)**2,'act1',
     1      l_act1,k_act1))
     1      call errquit('act1',101,MA_ERR)
        if (.not.MA_PUSH_GET(mt_int,(nos+nvs)**2,'prt1',
     1      l_prt1,k_prt1))
     1      call errquit('transform',102,MA_ERR)
        if (.not.MA_PUSH_GET(mt_byte,(nos+nvs)**4,'act2',
     1      l_act2,k_act2))
     1      call errquit('act2',103,MA_ERR)
        if (.not.MA_PUSH_GET(mt_byte,(nos+nvs)**4,'prt2',
     1      l_prt2,k_prt2))
     1      call errquit('prt2',104,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos+nvs)**2,'h',
     1      l_h,k_h))
     1      call errquit('h',200,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,((nos+nvs)/2)**2,'horb',
     1      l_horb,k_horb))
     1      call errquit('horb',201,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos+nvs)**4,'v',
     1      l_v,k_v))
     1      call errquit('v',202,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,nos*(nos+nvs),'t1',
     1      l_t1,k_t1))
     1      call errquit('t1',203,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos*nvs)**2,'t2',
     1      l_t2,k_t2))
     1      call errquit('t2',204,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos+nvs)**2,'ht',
     1      l_ht,k_ht))
     1      call errquit('ht',205,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,(nos+nvs)**4,'vt',
     1      l_vt,k_vt))
     1      call errquit('vt',206,MA_ERR)
c -------------------
           call ducc(rtdb,d_t1,k_t1_offset,d_t2,k_t2_offset,
     &    d_v2copy,k_v2copy_offset,dbl_mb(k_evl_sorted),
     &    nos,nvs,noas,nobs,nvas,nvbs,nactv,
     D        int_mb(k_transf),int_mb(k_act1),int_mb(k_prt1),
     D        byte_mb(k_act2),byte_mb(k_prt2),
     D        dbl_mb(k_h),dbl_mb(k_horb),dbl_mb(k_v),
     D        dbl_mb(k_t1),dbl_mb(k_t2),
     D        dbl_mb(k_ht),dbl_mb(k_vt))
              if (.not.ma_chop_stack(l_transf))
     1        call errquit('tce_energy: chopstack ducc',4,MA_ERR)
          endif ! qducc
c end of Tsource
c
c
c
c
c
            call tce_diis_tidy()
            if (save_t(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1 now...'
              endif
              call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                             handle_t1,irrep_t)
            endif
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
            if (save_t(3)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T3 now...'
              endif
              call x3_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                             handle_t3,irrep_t)
            endif
            goto 2327
!ccc            ccsdtleft=.true.
!ccc            goto 100
!edo            goto 330
          endif
          if (save_t(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1 now...'
            endif
            call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                           handle_t1,irrep_t)
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          if (save_t(3).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T3 now...'
            endif
            call x3_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                           handle_t3,irrep_t)
          endif
          call tce_diis(.false.,iter,.true.,.true.,.true.,.false.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  d_r3,d_t3,k_t3_offset,size_t3,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r3)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
c
c delete d_lambda related objects: ucc4
c
            if(ucc4) then
              call deletefile(d_lambda2)
              if (.not.ma_pop_stack(l_l2_offset))
     1        call errquit("tce_energy: l2_offset ucc4",4,MA_ERR)
            endif
c
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)

 2327 continue
c Tsource
      if(intorb) then
         d_v2=d_v2copy
         k_v2_offset=k_v2copy_offset
      endif !qducc
c Tsource 
! -------------
! CCSDT   left
! -------------
        if (left) then
        irrep_y=0
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9400) "CCSDT Lambda"
          call tce_filename('lr1',filename)
          call createfile(filename,d_lr1,size_l1)
          call ccsdt_lambda1(d_f1,d_lr1,d_t1,d_t2,d_t3,d_v2,d_lambda1,
     1       d_lambda2,d_lambda3,k_f1_offset,k_l1_offset,k_t1_offset,
     2       k_t2_offset,k_t3_offset,k_v2_offset,k_l1_offset,
     3       k_l2_offset,k_l3_offset)
          call reconcilefile(d_lr1,size_l1)
          call tce_filename('lr2',filename)
          call createfile(filename,d_lr2,size_l2)
          call ccsdt_lambda2(d_f1,d_lr2,d_t1,d_t2,d_t3,d_v2,d_lambda1,
     1       d_lambda2,d_lambda3,k_f1_offset,k_l2_offset,k_t1_offset,
     2       k_t2_offset,k_t3_offset,k_v2_offset,k_l1_offset,
     3       k_l2_offset,k_l3_offset)
          call reconcilefile(d_lr2,size_l2)
          call tce_filename('lr3',filename)
          call createfile(filename,d_lr3,size_l3)
          call ccsdt_lambda3(d_f1,d_lr3,d_t1,d_t2,d_v2,d_lambda1,
     1       d_lambda2,d_lambda3,k_f1_offset,k_l3_offset,k_t1_offset,
     2       k_t2_offset,k_v2_offset,k_l1_offset,
     3       k_l2_offset,k_l3_offset)
          call reconcilefile(d_lr3,size_l3)
          r1 = ddotfile(d_lr1,d_lr1,size_l1)
          r2 = ddotfile(d_lr2,d_lr2,size_l2)
          r3 = ddotfile(d_lr3,d_lr3,size_l3)
          residual = max(r1,r2,r3)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9410)
            endif
            if (ampnorms) then
              r1 = ddotfile(d_lambda1,d_lambda1,size_l1)
              r2 = ddotfile(d_lambda2,d_lambda2,size_l2)
              r3 = ddotfile(d_lambda3,d_lambda3,size_l3)
              if (nodezero) then
                write(LuOut,9082) "Lambda singles",r1
                write(LuOut,9082) "Lambda doubles",r2
                write(LuOut,9082) "Lambda triples",r3
              endif
            endif
            call deletefile(d_lr3)
            call deletefile(d_lr2)
            call deletefile(d_lr1)
            call tce_diis_tidy()
            if (save_l(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving L1 now...'
              endif
              call y1_restart_save(d_lambda1,k_l1_offset,size_l1,0,
     1                             handle_l1,irrep_y)
            endif
            if (save_l(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving L2 now...'
              endif
              call y2_restart_save(d_lambda2,k_l2_offset,size_l2,0,
     1                             handle_l2,irrep_y)
            endif
            if (save_l(3)) then
              if(nodezero) then
                write(LuOut,*) 'Saving L3 now...'
              endif
              call y3_restart_save(d_lambda3,k_l3_offset,size_l3,0,
     1                             handle_l3,irrep_y)
            endif
            do axis=1,3
              irrep_d=sym_abelian_axis(geom,axis)
              call tce_filename('d0',filename)
              call createfile(filename,d_d0,size_d0)
              call ccsdt_density1(d_d1(axis),d_d0,d_t1,d_t2,d_t3,
     1          d_lambda1,d_lambda2,d_lambda3,k_d1_offset(axis),
     2          k_d0_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     3          k_l1_offset,k_l2_offset,k_l3_offset)
              call reconcilefile(d_d0,size_d0)
              call get_block(d_d0,cdipole,1,0)
              dipole_cor(axis)=dipole_scf(axis)+cdipole
              call deletefile(d_d0)
            enddo
            dipole_cor_t = dsqrt( dipole_cor(1)*dipole_cor(1) +
     1                            dipole_cor(2)*dipole_cor(2) +
     2                            dipole_cor(3)*dipole_cor(3) )
            if (nodezero) write(LuOut,9430) "CCSDT",
     1                    dipole_cor(1),dipole_cor(1)*debye,
     2                    dipole_cor(2),dipole_cor(2)*debye,
     3                    dipole_cor(3),dipole_cor(3)*debye,
     4                    dipole_cor_t ,dipole_cor_t*debye
!cc            goto 2010
            goto 2337
          endif
          if (save_l(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L1 now...'
            endif
            call y1_restart_save(d_lambda1,k_l1_offset,size_l1,0,
     1                           handle_l1,irrep_y)
          endif
          if (save_l(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L2 now...'
            endif
            call y2_restart_save(d_lambda2,k_l2_offset,size_l2,0,
     1                           handle_l2,irrep_y)
          endif
          if (save_l(3).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L3 now...'
            endif
            call y3_restart_save(d_lambda3,k_l3_offset,size_l3,0,
     1                           handle_l3,irrep_y)
          endif
          call tce_diis(.true.,iter,.true.,.true.,.true.,.false.,
     1                  d_lr1,d_lambda1,k_l1_offset,size_l1,
     2                  d_lr2,d_lambda2,k_l2_offset,size_l2,
     3                  d_lr3,d_lambda3,k_l3_offset,size_l3,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_lr3)
          call deletefile(d_lr2)
          call deletefile(d_lr1)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
        endif ! (left)
 2337   continue
! ----------
! CCSDT(2)_Q
! ----------
        if (perturbative.eq.'2_q') then
          cpu = - util_cpusec()
          wall = - util_wallsec()
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call ccsdt2_q(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                  d_t3,k_t3_offset,d_f1,k_f1_offset,
     1                  d_v2,k_v2_offset,d_lambda1,k_l1_offset,
     2                  d_lambda2,k_l2_offset,d_lambda3,k_l3_offset,
     3                  d_e,k_e_offset,pt3_1)
          call deletefile(d_e)
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9070) "CCSDT(2)_Q",corr+pt3_1
            write(LuOut,9080) "CCSDT(2)_Q",ref+corr+pt3_1
            write(LuOut,9020) cpu, wall
          endif
          corr = corr + pt3_1
        endif
! -------------
! CCSDT-LR
! -------------
        if (lineresp) then
#if !defined(CCSDTLR)
      CALL ERRQUIT('CCSDTLR not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
!          if (nodezero) write(LuOut,*) "lineresp is go"
          do omegacount=1,anumfreq
!            if (nodezero) write(LuOut,*) "omegacount is go"
            omega = afreq(omegacount)
            if (omega.eq.(0.0d0)) then
              dynfreq = 1
            else
              dynfreq = 2
            endif
!
          do axis = 1, 3
!            if (nodezero) write(LuOut,*) "axis is go"
          if (respaxis(axis)) then
!            if (nodezero) write(LuOut,*) "respaxis(axis) is go"
            irrep_d=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            if (nodezero.and.util_print('mod1',print_default)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
            endif
            irrep_o=irrep_d
            irrep_x=irrep_d
!
          do omegasign=1,dynfreq
!            if (nodezero) write(LuOut,*) "omegasign is go"
            if (omegasign.eq.1) then
              omega = (1.0d0)*afreq(omegacount)
              dynaxis = 0
            elseif (omegasign.eq.2) then
              omega = (-1.0d0)*afreq(omegacount)
              dynaxis = 3
            endif
            if (nodezero) write(LuOut,9431) omega
!            if (nodezero) write(LuOut,*) "omega = ",omega
!            if (nodezero) write(LuOut,*) "dynaxis = ",dynaxis
!
            call tce_diis_init()
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSDT-LR"
!
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
              call tce_zero(d_rr1(axis),size_tr1(axis))
!        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega,d_tr1(axis+dynaxis),
     1             d_rr1(axis),size_tr1(axis))
!        write(LuOut,*) "eomccsdt_x1"
              call eomccsdt_x1(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
!        write(LuOut,*) "ccsdt_o1"
              call ccsdt_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     &             k_tr1_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset)
!
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
              call tce_zero(d_rr2(axis),size_tr2(axis))
!        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega,d_tr2(axis+dynaxis),
     1             d_rr2(axis),size_tr2(axis))
!        write(LuOut,*) "eomccsdt_x2"
              call eomccsdt_x2(d_f1,d_rr2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
!        write(LuOut,*) "ccsdt_o2"
              call ccsdt_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,d_t3,
     &             k_tr2_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset)
!
              call tce_filename(rr3filename(axis),filename)
              call createfile(filename,d_rr3(axis),size_tr3(axis))
              call tce_zero(d_rr3(axis),size_tr3(axis))
!        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega,d_tr3(axis+dynaxis),
     1             d_rr3(axis),size_tr3(axis))
!        write(LuOut,*) "eomccsdt_x3"
              call eomccsdt_x3(d_f1,d_rr3(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr3_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
!        write(LuOut,*) "ccsdt_o3"
              call ccsdt_o3(d_rr3(axis),d_d1(axis),d_t1,d_t2,d_t3,
     &             k_tr3_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset)
!
!        write(LuOut,*) "reconcilefile"
              call reconcilefile(d_rr1(axis),size_tr1(axis))
              call reconcilefile(d_rr2(axis),size_tr2(axis))
              call reconcilefile(d_rr3(axis),size_tr3(axis))
!        write(LuOut,*) "tce_residual_tr1"
              call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
!        write(LuOut,*) "tce_residual_tr2"
              call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
!        write(LuOut,*) "tce_residual_tr3"
              call tce_residual_tr3(d_rr3(axis),k_tr3_offset(axis),rr3)
              residual = max(rr1,rr2,rr3)
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
              if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
              if (residual .lt. thresh) then
                if (nodezero) then
                  write(LuOut,9410)
                endif
                if (ampnorms) then
                  call tce_residual_tr1(d_tr1(axis+dynaxis),
     1                                  k_tr1_offset(axis),rr1)
                  call tce_residual_tr2(d_tr2(axis+dynaxis),
     1                                  k_tr2_offset(axis),rr2)
                  call tce_residual_tr3(d_tr3(axis+dynaxis),
     1                                  k_tr3_offset(axis),rr3)
                  if (nodezero) then
                    write(LuOut,9082) "T(1) singles",rr1
                    write(LuOut,9082) "T(1) doubles",rr2
                    write(LuOut,9082) "T(1) triples",rr3
                  endif
                endif
                call deletefile(d_rr3(axis))
                call deletefile(d_rr2(axis))
                call deletefile(d_rr1(axis))
                call tce_diis_tidy()
                if (save_tr(1)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving T1(1) now...'
                  endif
                  call tr1_restart_save(d_tr1(axis+dynaxis),
     1                 k_tr1_offset(axis),size_tr1(axis),
     2                 axis+dynaxis,handle_tr1(axis),irrep_x)
                endif
                if (save_tr(2)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving T2(1) now...'
                  endif
                  call tr2_restart_save(d_tr2(axis+dynaxis),
     1                 k_tr2_offset(axis),size_tr2(axis),
     2                 axis+dynaxis,handle_tr2(axis),irrep_x)
                endif
                if (save_tr(3)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving T3(1) now...'
                  endif
                  call tr3_restart_save(d_tr3(axis),k_tr3_offset(axis),
     1                 size_tr3(axis),axis,handle_tr3(axis),irrep_x)
                endif
                goto 2338
              endif
              if (save_tr(1).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving T1(1) now...'
                endif
                call tr1_restart_save(d_tr1(axis+dynaxis),
     1               k_tr1_offset(axis),size_tr1(axis),
     2               axis+dynaxis,handle_tr1(axis),irrep_x)
              endif
              if (save_tr(2).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving T2(1) now...'
                endif
                call tr2_restart_save(d_tr2(axis+dynaxis),
     1               k_tr2_offset(axis),size_tr2(axis),
     2               axis+dynaxis,handle_tr2(axis),irrep_x)
              endif
              if (save_tr(3).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving T3(1) now...'
                endif
                call tr3_restart_save(d_tr3(axis),k_tr3_offset(axis),
     1               size_tr3(axis),axis,handle_tr3(axis),irrep_x)
              endif
              call tce_diis2(.false.,iter,.true.,.true.,.true.,.false.,
     1             d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2             size_tr1(axis),
     3             d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4             size_tr2(axis),
     5             d_rr3(axis),d_tr3(axis+dynaxis),k_tr3_offset(axis),
     6             size_tr3(axis),
     6             dummy,dummy,dummy,dummy,omega)
              call deletefile(d_rr3(axis))
              call deletefile(d_rr2(axis))
              call deletefile(d_rr1(axis))
              if (nodezero) call util_flush(LuOut)
            enddo ! iter loop
          call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 2338 continue
!
          enddo ! omegasign loop
          endif ! respaxis(axis)
          enddo ! axis loop
!
!         Evaluate Dipole Polarizability
!
          do axisA = 1, 3
          do axisB = 1, axisA
          alpha(axisA,axisB)=0.0d0
          if (respaxis(axisA).and.respaxis(axisB)) then
            irrep_a=sym_abelian_axis(geom,axisA)
            irrep_b=sym_abelian_axis(geom,axisB)
            irrep_y=irrep_g
!
!            if (nodezero) write(LuOut,*) "axisA = ",axisA
!            if (nodezero) write(LuOut,*) "axisB = ",axisB
!
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,size_a0)
!
              alpha1=0.0d0
              alpha2=0.0d0
!
              if (omega.ne.0.0d0) then
                dynaxis = 3
              else
                dynaxis = 0
              endif
!
              irrep_c=irrep_a
              irrep_oa=irrep_a
              irrep_ob=irrep_b
!              if (nodezero) write(LuOut,*) "ccsdt_lr_alpha_offdiag #1"
              call ccsdt_lr_alpha_offdiag(d_tr1(axisA),
     1             d_tr2(axisA),d_tr3(axisA),
     2             d_tr1(axisB+dynaxis),d_tr2(axisB+dynaxis),
     3             d_tr3(axisB+dynaxis),d_tr1(axisA),
     4             d_tr2(axisA),d_tr3(axisA),
     5             d_f1,d_a0,d_d1(axisA),d_d1(axisB),d_t1,d_t2,d_t3,
     6             d_v2,d_lambda1,d_lambda2,d_lambda3,
     7             k_tr1_offset(axisA),k_tr2_offset(axisA),
     8             k_tr3_offset(axisA),k_tr1_offset(axisB),
     9             k_tr2_offset(axisB),k_tr3_offset(axisB),
     1             k_tr1_offset(axisA),k_tr2_offset(axisA),
     2             k_tr3_offset(axisA),k_f1_offset,k_a0_offset,
     3             k_d1_offset(axisA),k_d1_offset(axisB),
     4             k_t1_offset,k_t2_offset,k_t3_offset,
     5             k_v2_offset,k_l1_offset,k_l2_offset,k_l3_offset)
!
              call reconcilefile(d_a0,size_a0)
              call get_block(d_a0,alpha1,1,0)
              call tce_zero(d_a0,size_a0)
!              write(LuOut,*) "alpha1 = ",alpha1
!
              if ((omega.eq.0.0d0).or.(axisA.eq.axisB)) then
                alpha2=alpha1
              else
!                if (nodezero) write(LuOut,*) "ccsdt_lr_alpha_offdiag #2"
                call ccsdt_lr_alpha_offdiag(d_tr1(axisA+dynaxis),
     1               d_tr2(axisA+dynaxis),d_tr3(axisA+dynaxis),
     2               d_tr1(axisB),d_tr2(axisB),
     3               d_tr3(axisB),d_tr1(axisA+dynaxis),
     4               d_tr2(axisA+dynaxis),d_tr3(axisA+dynaxis),
     5               d_f1,d_a0,d_d1(axisA),d_d1(axisB),d_t1,d_t2,d_t3,
     6               d_v2,d_lambda1,d_lambda2,d_lambda3,
     7               k_tr1_offset(axisA),k_tr2_offset(axisA),
     8               k_tr3_offset(axisA),k_tr1_offset(axisB),
     9               k_tr2_offset(axisB),k_tr3_offset(axisB),
     1               k_tr1_offset(axisA),k_tr2_offset(axisA),
     2               k_tr3_offset(axisA),k_f1_offset,k_a0_offset,
     3               k_d1_offset(axisA),k_d1_offset(axisB),
     4               k_t1_offset,k_t2_offset,k_t3_offset,
     5               k_v2_offset,k_l1_offset,k_l2_offset,k_l3_offset)
!
                call reconcilefile(d_a0,size_a0)
                call get_block(d_a0,alpha2,1,0)
                call tce_zero(d_a0,size_a0)
              endif ! axisA.eq.axisB
!              write(LuOut,*) "alpha2 = ",alpha2
!
              alpha(axisA,axisB)=(-0.5d0)*(alpha1+alpha2)
!
            call deletefile(d_a0)
!
!            write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
!            if (nodezero) write(LuOut,9020) cpu, wall
!
          endif ! respaxis(axis)
          enddo ! axisB loop
          enddo ! axisA loop
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
!
          do i = 1, 3
            do j = 1, i
              alphacopy(i,j)=alpha(i,j)
              alphacopy(j,i)=alpha(i,j)
            enddo
          enddo
!
          call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
!
          alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
          alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
          alphaani = dsqrt(0.5d0*alphaani)
!
          if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1        .and.nodezero) write(LuOut,9911)
!
          if (nodezero) write(LuOut,9434) "CCSDT Linear Response",
     1      afreq(omegacount),auXnm/afreq(omegacount),
     2      alpha(1,1),alpha(2,1),alpha(3,1),
     3      au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4      alpha(2,1),alpha(2,2),alpha(3,2),
     5      au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6      alpha(3,1),alpha(3,2),alpha(3,3),
     7      au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8      alphaevr(1),alphaevr(2),alphaevr(3),
     9      au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1      alphaiso,au2ang3*alphaiso,
     2      alphaani,au2ang3*alphaani
          if (nodezero) write(LuOut,9020) cpu, wall
          call util_flush(LuOut)
!
          enddo ! omegacount loop
!
#endif
        endif ! lineresp
!
! -------------
! CCSDT-YR
! -------------
        if (leftresp) then
#if !defined(CCSDTLR)
      CALL ERRQUIT('CCSDTLR not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
!
          do omegacount=1,anumfreq
            omega1 = afreq(omegacount)
            if (omega1.eq.(0.0d0)) then
              dynfreq = 1
            else
              dynfreq = 2
            endif
!
          do axis = 1, 3
          if (respaxis(axis)) then
            irrep_d=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            if (nodezero.and.util_print('mod1',print_default)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
            endif
            irrep_o=irrep_d
            irrep_x=irrep_d
!
          do omegasign=1,dynfreq
            if (omegasign.eq.1) then
              omega1 = (1.0d0)*afreq(omegacount)
              dynaxis = 0
            elseif (omegasign.eq.2) then
              omega1 = (-1.0d0)*afreq(omegacount)
              dynaxis = 3
            endif
            if (nodezero) write(LuOut,9431) omega1
!            if (nodezero) write(LuOut,*) "omega1 = ",omega1
!            if (nodezero) write(LuOut,*) "dynaxis = ",dynaxis
!
            call tce_diis_init()
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSDT-YR (Static)"
!        write(LuOut,*) "tce_filename"
              call tce_filename(yrr1filename(axis),filename)
!        write(LuOut,*) "createfile"
              call createfile(filename,d_yrr1(axis),size_yr1(axis))
!        write(LuOut,*) "tce_zero"
              call tce_zero(d_yrr1(axis),size_yr1(axis))
              call daxpyfile(1,(1.0d0)*omega1,d_yr1(axis+dynaxis),
     1             d_yrr1(axis),size_yr1(axis))
!        write(LuOut,*) "ccsdt_y_o1"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_lambda_o1(d_yrr1(axis),d_d1(axis),d_t1,d_t2,
     1             d_t3,d_lambda1,d_lambda2,d_lambda3,
     2             k_yr1_offset(axis),k_d1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,
     3             k_l1_offset,k_l2_offset,k_l3_offset)
!              call ga_print(d_yrr1(axis))
!        write(LuOut,*) "ccsdt_y_tr1"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_y_tr1(d_f1,d_yrr1(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_lambda1,d_lambda2,d_lambda3,
     2             k_f1_offset,k_yr1_offset(axis),k_t1_offset,
     3             k_t2_offset,k_v2_offset,k_tr1_offset(axis),
     4             k_tr2_offset(axis),k_tr3_offset(axis),
     5             k_l1_offset,k_l2_offset,k_l3_offset)
!              call ga_print(d_yrr1(axis))
!        write(LuOut,*) "eomccsdt_y1"
              irrep_y=irrep_d ! using Lambda(1) amplitudes
              call eomccsdt_y1(d_f1,d_yrr1(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_yr1(axis+dynaxis),d_yr2(axis+dynaxis),
     2             d_yr3(axis+dynaxis),k_f1_offset,k_yr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_yr1_offset(axis),k_yr2_offset(axis),
     5             k_yr3_offset(axis))
!              call ga_print(d_yrr1(axis))
!        write(LuOut,*) "tce_filename"
              call tce_filename(yrr2filename(axis),filename)
!        write(LuOut,*) "createfile"
              call createfile(filename,d_yrr2(axis),size_yr2(axis))
!        write(LuOut,*) "tce_zero"
              call tce_zero(d_yrr2(axis),size_yr2(axis))
              call daxpyfile(1,(1.0d0)*omega1,d_yr2(axis+dynaxis),
     1             d_yrr2(axis),size_yr2(axis))
!        write(LuOut,*) "ccsdt_y_o2"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_lambda_o2(d_yrr2(axis),d_d1(axis),d_t1,d_t2,
     1             d_lambda1,d_lambda2,d_lambda3,k_yr2_offset(axis),
     2             k_d1_offset(axis),k_t1_offset,k_t2_offset,
     3             k_l1_offset,k_l2_offset,k_l3_offset)
!        write(LuOut,*) "ccsdt_y_tr2"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_y_tr2(d_f1,d_yrr2(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_lambda1,d_lambda2,
     2             d_lambda3,k_f1_offset,k_yr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis),k_l1_offset,k_l2_offset,
     6             k_l3_offset)
!        write(LuOut,*) "eomccsdt_y2"
              irrep_y=irrep_d ! using Lambda(1) amplitudes
              call eomccsdt_y2(d_f1,d_yrr2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_yr1(axis+dynaxis),d_yr2(axis+dynaxis),
     2             d_yr3(axis+dynaxis),k_f1_offset,k_yr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_yr1_offset(axis),k_yr2_offset(axis),
     5             k_yr3_offset(axis))
!        write(LuOut,*) "tce_filename"
              call tce_filename(yrr3filename(axis),filename)
!        write(LuOut,*) "createfile"
              call createfile(filename,d_yrr3(axis),size_yr3(axis))
!        write(LuOut,*) "tce_zero"
              call tce_zero(d_yrr2(axis),size_yr2(axis))
              call daxpyfile(1,(1.0d0)*omega1,d_yr2(axis+dynaxis),
     1             d_yrr2(axis),size_yr2(axis))
!        write(LuOut,*) "ccsdt_y_o3"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_lambda_o3(d_yrr3(axis),d_d1(axis),d_t1,
     1             d_lambda2,d_lambda3,k_yr3_offset(axis),
     2             k_d1_offset(axis),k_t1_offset,
     3             k_l2_offset,k_l3_offset)
!        write(LuOut,*) "ccsdt_y_tr3"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_y_tr3(d_f1,d_yrr3(axis),d_t1,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),d_lambda2,
     2             d_lambda3,k_f1_offset,k_yr3_offset(axis),
     3             k_t1_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_l2_offset,k_l3_offset)
!        write(LuOut,*) "eomccsdt_y3"
              irrep_y=irrep_d ! using Lambda(1) amplitudes
              call eomccsdt_y3(d_f1,d_yrr3(axis),d_t1,d_t2,d_v2,
     1             d_yr1(axis+dynaxis),d_yr2(axis+dynaxis),
     2             d_yr3(axis+dynaxis),k_f1_offset,k_yr3_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_yr1_offset(axis),k_yr2_offset(axis),
     5             k_yr3_offset(axis))
!        write(LuOut,*) "reconcilefile"
              call reconcilefile(d_yrr1(axis),size_yr1(axis))
              call reconcilefile(d_yrr2(axis),size_yr2(axis))
              call reconcilefile(d_yrr3(axis),size_yr3(axis))
!        write(LuOut,*) "ddotfile"
              rr1 = ddotfile(d_yrr1(axis),d_yrr1(axis),size_yr1(axis))
              rr2 = ddotfile(d_yrr2(axis),d_yrr2(axis),size_yr2(axis))
              rr3 = ddotfile(d_yrr3(axis),d_yrr3(axis),size_yr3(axis))
              residual = max(rr1,rr2,rr3)
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
              if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
              if (residual .lt. thresh) then
                if (nodezero) then
                  write(LuOut,9410)
                endif
                if (ampnorms) then
                  r1 = ddotfile(d_yr1(axis+dynaxis),
     1                          d_yr1(axis+dynaxis),size_yr1(axis))
                  r2 = ddotfile(d_yr2(axis+dynaxis),
     1                          d_yr2(axis+dynaxis),size_yr2(axis))
                  r3 = ddotfile(d_yr3(axis+dynaxis),
     1                          d_yr3(axis+dynaxis),size_yr3(axis))
                  if (nodezero) then
                    write(LuOut,9082) "Lambda(1) singles",r1
                    write(LuOut,9082) "Lambda(1) doubles",r2
                    write(LuOut,9082) "Lambda(1) triples",r3
                  endif
                endif
!
!                call ga_print(d_yrr1(axis))
!                call ga_print(d_yrr2(axis))
!                call ga_print(d_yrr3(axis))
!
                call deletefile(d_yrr3(axis))
                call deletefile(d_yrr2(axis))
                call deletefile(d_yrr1(axis))
                call tce_diis_tidy()
                if (save_lr(1)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving L1(1) now...'
                  endif
                  call yr1_restart_save(d_yr1(axis),k_yr1_offset(axis),
     1                 size_yr1(axis),axis,handle_yr1(axis),irrep_d)
                endif
                if (save_lr(2)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving L2(1) now...'
                  endif
                  call yr2_restart_save(d_yr2(axis),k_yr2_offset(axis),
     1                 size_yr2(axis),axis,handle_yr2(axis),irrep_d)
                endif
                if (save_lr(3)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving L3(1) now...'
                  endif
                  call yr3_restart_save(d_yr3(axis),k_yr3_offset(axis),
     1                 size_yr3(axis),axis,handle_yr3(axis),irrep_d)
                endif
                goto 2339
              endif
              if (save_lr(1).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving L1(1) now...'
                endif
                call yr1_restart_save(d_yr1(axis),k_yr1_offset(axis),
     1               size_yr1(axis),axis,handle_yr1(axis),irrep_d)
              endif
              if (save_lr(2).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving L2(1) now...'
                endif
                call yr2_restart_save(d_yr2(axis),k_yr2_offset(axis),
     1               size_yr2(axis),axis,handle_yr2(axis),irrep_d)
              endif
              if (save_lr(3).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving L3(1) now...'
                endif
                call yr3_restart_save(d_yr3(axis),k_yr3_offset(axis),
     1               size_yr3(axis),axis,handle_yr3(axis),irrep_d)
              endif
!       write(LuOut,*) "before tce_diis2"
              irrep_y=irrep_d
              call tce_diis2(.true.,iter,.true.,.true.,.false.,.false.,
     1             d_yrr1(axis),d_yr1(axis+dynaxis),k_yr1_offset(axis),
     2             size_yr1(axis),
     3             d_yrr2(axis),d_yr2(axis+dynaxis),k_yr2_offset(axis),
     4             size_yr2(axis),
     5             d_yrr3(axis),d_yr3(axis+dynaxis),k_yr3_offset(axis),
     6             size_yr3(axis),
     6             dummy,dummy,dummy,dummy,omega)
!       write(LuOut,*) "after tce_diis2"
              call deletefile(d_yrr3(axis))
              call deletefile(d_yrr2(axis))
              call deletefile(d_yrr1(axis))
              if (nodezero) call util_flush(LuOut)
            enddo ! iter loop
          call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 2339 continue
!
          enddo ! omegasign loop
!
          endif ! respaxis(axis)
          enddo ! axis loop
          call errquit('tce_energy: manual stop',911,CALC_ERR)
!
!         EVALUATE POLARIZABILITY
!
          cpu=-util_cpusec()
          wall=-util_wallsec()
          do axisA = 1, 3
          irrep_a=sym_abelian_axis(geom,axisA)
          do axisB = 1, 3
          irrep_b=sym_abelian_axis(geom,axisB)
          alpha(axisA,axisB)=0.0d0
          if ( respaxis(axisA).and.respaxis(axisB) ) then
!
            alpha1=0.0d0
            alpha2=0.0d0
            alpha3=0.0d0
            alpha4=0.0d0
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,size_a0)
!
            irrep_y=irrep_g
            irrep_d=irrep_a
            irrep_tr=irrep_b
!
            call alpha_1(d_d1(axisA),d_a0,d_t1,d_t2,
     1           d_tr1(axisB),d_tr2(axisB),
     2           d_lambda1,d_lambda2,k_d1_offset(axisA),k_a0_offset,
     4           k_t1_offset,k_t2_offset,k_tr1_offset(axisB),
     5           k_tr2_offset(axisB),k_l1_offset,k_l2_offset)
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha1,1,0)
            call tce_zero(d_a0,size_a0)
!
            if (omega1.eq.0.0d0) then
              alpha2=alpha1
            else
              call alpha_1(d_d1(axisA),d_a0,d_t1,d_t2,
     1             d_tr1(axisB+3),d_tr2(axisB+3),
     2             d_lambda1,d_lambda2,k_d1_offset(axisA),k_a0_offset,
     4             k_t1_offset,k_t2_offset,k_tr1_offset(axisB),
     5             k_tr2_offset(axisB),k_l1_offset,k_l2_offset)
              call reconcilefile(d_a0,size_a0)
              call get_block(d_a0,alpha2,1,0)
              call tce_zero(d_a0,size_a0)
            endif
!
            irrep_o=irrep_a
            irrep_y=irrep_b
!
            call alpha_2_yr(d_a0,d_d1(axisA),d_t1,d_t2,
     1           d_yr1(axisB),d_yr2(axisB),
     2           k_a0_offset,k_d1_offset(axisA),k_t1_offset,
     3           k_t2_offset,k_yr1_offset(axisB),k_yr2_offset(axisB))
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha3,1,0)
            call tce_zero(d_a0,size_a0)
!
            if (omega1.eq.0.0d0) then
              alpha4=alpha3
            else
              call alpha_2_yr(d_a0,d_d1(axisA),d_t1,d_t2,
     1             d_yr1(axisB+3),d_yr2(axisB+3),
     2             k_a0_offset,k_d1_offset(axisA),k_t1_offset,
     3             k_t2_offset,k_yr1_offset(axisB),k_yr2_offset(axisB))
              call reconcilefile(d_a0,size_a0)
              call get_block(d_a0,alpha4,1,0)
              call tce_zero(d_a0,size_a0)
            endif
!
!            write(LuOut,*) "alpha1 = ",alpha1
!            write(LuOut,*) "alpha2 = ",alpha2
!            write(LuOut,*) "alpha3 = ",alpha3
!            write(LuOut,*) "alpha4 = ",alpha4
            alpha(axisA,axisB)=(-0.5d0)*(alpha1+alpha2+alpha3+alpha4)
            call deletefile(d_a0)
!
          endif ! respaxis(axisA,axisB)
          enddo ! axisB loop
          enddo ! axisA loop
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
!
          do j = 1, 3
            do i = 1, 3
              alphacopy(i,j)=alpha(i,j)
            enddo
          enddo
!
          call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
!
          alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
          alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
          alphaani = dsqrt(0.5d0*alphaani)
!
          if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1        .and.nodezero) write(LuOut,9911)
!
          if (nodezero) write(LuOut,9434) "CCSDT Lambda Response",
     1      afreq(omegacount),auXnm/afreq(omegacount),
     2      alpha(1,1),alpha(2,1),alpha(3,1),
     3      au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4      alpha(2,1),alpha(2,2),alpha(3,2),
     5      au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6      alpha(3,1),alpha(3,2),alpha(3,3),
     7      au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8      alphaevr(1),alphaevr(2),alphaevr(3),
     9      au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1      alphaiso,au2ang3*alphaiso,
     2      alphaani,au2ang3*alphaani
          if (nodezero) write(LuOut,9020) cpu, wall
          call util_flush(LuOut)
!
!
!         EVALUATE HYPERPOLARIZABILITY
!
          cpu=-util_cpusec()
          wall=-util_wallsec()
          do axisA = 1, 3
          do axisB = 1, 3
          do axisC = 1, 3
            beta(axisA,axisB,axisC)=0.0d0
          enddo
          enddo
          enddo
!
          do axisA = 1, 3
            irrep_a=sym_abelian_axis(geom,axisA)
          do axisB = 1, 3
            irrep_b=sym_abelian_axis(geom,axisB)
          do axisC = 1, 3
            irrep_c=sym_abelian_axis(geom,axisC)
          if ( ( respaxis(axisA).and.respaxis(axisB) )
     &                          .and.respaxis(axisC) ) then
!
!            write(LuOut,*) "axisA = ",axisA
!            write(LuOut,*) "axisB = ",axisB
!            write(LuOut,*) "axisC = ",axisC
!            call util_flush(LuOut)
!
            beta1=0.0d0
            beta2=0.0d0
            beta3=0.0d0
            beta4=0.0d0
            call tce_filename('b0',filename)
            call createfile(filename,d_b0,size_b0)
!
            irrep_y=irrep_g
            irrep_tra=irrep_a
            irrep_trb=irrep_b
            irrep_trc=irrep_c
!            write(LuOut,*) "beta_1"
            call beta_1(d_b0,d_t1,d_tr1(axisA),d_tr2(axisA),
     1           d_tr1(axisB),d_tr2(axisB),d_tr1(axisC),d_tr2(axisC),
     2           d_v2,d_lambda1,d_lambda2,k_b0_offset,k_t1_offset,
     3           k_tr1_offset(axisA),k_tr2_offset(axisA),
     4           k_tr1_offset(axisB),k_tr2_offset(axisB),
     5           k_tr1_offset(axisC),k_tr2_offset(axisC),
     6           k_v2_offset,k_l1_offset,k_l2_offset)
            call reconcilefile(d_b0,size_b0)
            call get_block(d_b0,beta1,1,0)
            call tce_zero(d_b0,size_b0)
!
            irrep_y=irrep_g
            irrep_d=irrep_a
            irrep_tra=irrep_b
            irrep_trb=irrep_c
!            write(LuOut,*) "beta_2"
            call beta_2(d_d1(axisA),d_b0,d_tr1(axisB),d_tr2(axisB),
     1           d_tr1(axisC),d_tr2(axisC),d_lambda1,d_lambda2,
     2           k_d1_offset(axisA),k_b0_offset,
     3           k_tr1_offset(axisB),k_tr2_offset(axisB),
     4           k_tr1_offset(axisC),k_tr2_offset(axisC),
     5           k_l1_offset,k_l2_offset)
            call reconcilefile(d_b0,size_b0)
            call get_block(d_b0,beta2,1,0)
            call tce_zero(d_b0,size_b0)
!
            irrep_y=irrep_a
            irrep_tra=irrep_b
            irrep_trb=irrep_c
!            write(LuOut,*) "beta_3"
            call beta_3(d_f1,d_b0,d_t1,d_t2,d_tr1(axisB),d_tr2(axisB),
     1           d_tr1(axisC),d_tr2(axisC),d_v2,
     2           d_yr1(axisA),d_yr2(axisA),
     3           k_f1_offset,k_b0_offset,k_t1_offset,k_t2_offset,
     4           k_tr1_offset(axisB),k_tr2_offset(axisB),
     5           k_tr1_offset(axisC),k_tr2_offset(axisC),k_v2_offset,
     6           k_yr1_offset(axisA),k_yr2_offset(axisA))
            call reconcilefile(d_b0,size_b0)
            call get_block(d_b0,beta3,1,0)
            call tce_zero(d_b0,size_b0)
!
            irrep_y=irrep_a
            irrep_d=irrep_b
            irrep_tr=irrep_c
!            write(LuOut,*) "beta_4"
            call beta_4(d_d1(axisB),d_b0,d_t1,d_t2,
     1           d_tr1(axisC),d_tr2(axisC),d_yr1(axisA),d_yr2(axisA),
     2           k_d1_offset(axisB),k_b0_offset,k_t1_offset,k_t2_offset,
     3           k_tr1_offset(axisC),k_tr2_offset(axisC),
     4           k_yr1_offset(axisA),k_yr2_offset(axisA))
            call reconcilefile(d_b0,size_b0)
            call get_block(d_b0,beta4,1,0)
            call tce_zero(d_b0,size_b0)
!
!            write(LuOut,*) "beta1 = ",beta1
!            write(LuOut,*) "beta2 = ",beta2
!            write(LuOut,*) "beta3 = ",beta3
!            write(LuOut,*) "beta4 = ",beta4
!
!           This term contributes to only one component
!
            beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta1
!
!           These terms both contribute to three components
!
            beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta2
            beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta2
            beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta2
!
            beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta3
            beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta3
            beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta3
!
!           This term contributes to six components
!
            beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta4
            beta(axisA,axisC,axisB)=beta(axisA,axisC,axisB)-beta4
            beta(axisB,axisA,axisC)=beta(axisB,axisA,axisC)-beta4
            beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta4
            beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta4
            beta(axisC,axisB,axisA)=beta(axisC,axisB,axisA)-beta4
!
            call deletefile(d_b0)
!
          endif ! respaxis(axisA,axisB,axisC)
          enddo ! axisB loop
          enddo ! axisB loop
          enddo ! axisA loop
!
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9120) "CCSDT Quadratic Response"
            write(LuOut,9504)
            write(LuOut,9502)
            do axisA = 1, 3
            do axisB = 1, 3
            do axisC = 1, 3
              if (dabs(beta(axisA,axisB,axisC)).gt.1.0d-8)
     1        write(LuOut,9503) axisname(axisA),axisname(axisB),
     2                          axisname(axisC),beta(axisA,axisB,axisC),
     3                          ' / au'
            enddo
            enddo
            enddo
            write(LuOut,9502)
          endif
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9020) cpu, wall
          call util_flush(LuOut)
!
          enddo ! omegacount loop
! CCSDTLR
#endif
        endif ! leftresp
! CCSDT
#endif
      endif
      if (model.eq.'ccsdta') then
#if !defined(CCSDT)
        CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
! -----------------------------
! CCSDt (CCSDT-active)   right
! -----------------------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CCSDt"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call ccsdt_e(d_f1,d_e,d_t1,d_t2,d_v2,
     1                 k_f1_offset,k_e_offset,
     2                 k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_t1)
          call ccsdt_t1a(d_f1,d_r1,d_t1,d_t2,d_t3,d_v2,
     1                  k_f1_offset,k_t1_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset)
          call reconcilefile(d_r1,size_t1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call ccsdt_t2a(d_f1,d_r2,d_t1,d_t2,d_t3,d_v2,
     1                  k_f1_offset,k_t2_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset)
          call reconcilefile(d_r2,size_t2)
          call tce_filename('r3',filename)
          call createfile(filename,d_r3,size_t3)
          call ccsdt_t3a(d_f1,d_r3,d_t1,d_t2,d_t3,d_v2,
     1                  k_f1_offset,k_t3_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset)
          call reconcilefile(d_r3,size_t3)
          call tce_residual_t1(d_r1,k_t1_offset,r1)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          call tce_residual_t3a(d_r3,k_t3_offset,r3)
          residual = max(r1,r2,r3)
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "CCSDt",corr
              write(LuOut,9080) "CCSDt",ref + corr
            endif
            call deletefile(d_r3)
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            call tce_diis_tidy()
            if (save_t(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1 now...'
              endif
              call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                             handle_t1,irrep_t)
            endif
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
            if (save_t(3)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T3a now...'
              endif
              call x3a_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                              handle_t3,irrep_t)
            endif
            ccsdtleft=.true.
            goto 100
!edo            goto 330
          endif
          if (save_t(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1 now...'
            endif
            call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                           handle_t1,irrep_t)
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          if (save_t(3).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T3a now...'
            endif
            call x3a_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                            handle_t3,irrep_t)
          endif
          call tce_diis(.false.,iter,.true.,.true.,.true.,.false.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  d_r3,d_t3,k_t3_offset,size_t3,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r3)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
#endif
      endif
      if (model.eq.'ccsdtq') then
#if !defined(CCSDTQ)
      CALL ERRQUIT('CCSDTQ not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
! -------------
! CCSDTQ  right
! -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CCSDTQ"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call ccsdtq_e(d_f1,d_e,d_t1,d_t2,d_v2,
     1                  k_f1_offset,k_e_offset,
     2                  k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_t1)
          call ccsdtq_t1(d_f1,d_r1,d_t1,d_t2,d_t3,d_v2,
     1                   k_f1_offset,k_t1_offset,k_t1_offset,
     2                   k_t2_offset,k_t3_offset,k_v2_offset)
          call reconcilefile(d_r1,size_t1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call ccsdtq_t2(d_f1,d_r2,d_t1,d_t2,d_t3,d_t4,d_v2,
     1                   k_f1_offset,k_t2_offset,k_t1_offset,
     2                   k_t2_offset,k_t3_offset,k_t4_offset,
     3                   k_v2_offset)
          call reconcilefile(d_r2,size_t2)
          call tce_filename('r3',filename)
          call createfile(filename,d_r3,size_t3)
          call ccsdtq_t3(d_f1,d_r3,d_t1,d_t2,d_t3,d_t4,d_v2,
     1                   k_f1_offset,k_t3_offset,k_t1_offset,
     2                   k_t2_offset,k_t3_offset,k_t4_offset,
     3                   k_v2_offset)
          call reconcilefile(d_r3,size_t3)
          call tce_filename('r4',filename)
          call createfile(filename,d_r4,size_t4)
          call ccsdtq_t4(d_f1,d_r4,d_t1,d_t2,d_t3,d_t4,d_v2,
     1                   k_f1_offset,k_t4_offset,k_t1_offset,
     2                   k_t2_offset,k_t3_offset,k_t4_offset,
     3                   k_v2_offset)
          call reconcilefile(d_r4,size_t4)
          call tce_residual_t1(d_r1,k_t1_offset,r1)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          call tce_residual_t3(d_r3,k_t3_offset,r3)
          call tce_residual_t4(d_r4,k_t4_offset,r4)
          residual = max(r1,r2,r3,r4)
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "CCSDTQ",corr
              write(LuOut,9080) "CCSDTQ",ref + corr
            endif
            if (diagnt1) call tce_diagnose_t1(d_t1,k_t1_offset,t1diagn)
            if (ampnorms) then
              call tce_residual_t1(d_t1,k_t1_offset,r1)
              call tce_residual_t2(d_t2,k_t2_offset,r2)
              call tce_residual_t3(d_t3,k_t3_offset,r3)
              call tce_residual_t4(d_t4,k_t4_offset,r4)
              if (nodezero) then
                write(LuOut,9082) "T singles",r1
                write(LuOut,9082) "T doubles",r2
                write(LuOut,9082) "T triples",r3
                write(LuOut,9082) "T quadruples",r4
              endif
            endif
            call deletefile(d_r4)
            call deletefile(d_r3)
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            call tce_print_x1(d_t1,k_t1_offset,printtol,irrep_t)
            call tce_print_x2(d_t2,k_t2_offset,printtol,irrep_t)
            call tce_diis_tidy()
            if (save_t(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1 now...'
              endif
              call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                             handle_t1,irrep_t)
            endif
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
            if (save_t(3)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T3 now...'
              endif
              call x3_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                             handle_t3,irrep_t)
            endif
            if (save_t(4)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T4 now...'
              endif
              call x4_restart_save(d_t4,k_t4_offset,size_t4,0,
     1                             handle_t4,irrep_t)
            endif
            goto 340
          endif
          if (save_t(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1 now...'
            endif
            call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                           handle_t1,irrep_t)
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          if (save_t(3).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T3 now...'
            endif
            call x3_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                           handle_t3,irrep_t)
          endif
          if (save_t(4).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T4 now...'
            endif
            call x4_restart_save(d_t4,k_t4_offset,size_t4,0,
     1                           handle_t4,irrep_t)
          endif
          call tce_diis(.false.,iter,.true.,.true.,.true.,.true.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  d_r3,d_t3,k_t3_offset,size_t3,
     4                  d_r4,d_t4,k_t4_offset,size_t4)
          call deletefile(d_r4)
          call deletefile(d_r3)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
! -------------
! CCSDTQ  left
! -------------
  340   continue
        if (left) then
        irrep_y=0
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1))
     1    write(LuOut,9400) "CCSDTQ Lambda"
          call tce_filename('lr1',filename)
          call createfile(filename,d_lr1,size_l1)
          call ccsdtq_lambda1(d_f1,d_lr1,d_t1,d_t2,d_t3,d_t4,d_v2,
     1      d_lambda1,d_lambda2,d_lambda3,d_lambda4,k_f1_offset,
     2      k_l1_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     3      k_t4_offset,k_v2_offset,k_l1_offset,k_l2_offset,
     4      k_l3_offset,k_l4_offset)
          call reconcilefile(d_lr1,size_l1)
          call tce_filename('lr2',filename)
          call createfile(filename,d_lr2,size_l2)
          call ccsdtq_lambda2(d_f1,d_lr2,d_t1,d_t2,d_t3,d_t4,d_v2,
     1      d_lambda1,d_lambda2,d_lambda3,d_lambda4,k_f1_offset,
     2      k_l2_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     3      k_t4_offset,k_v2_offset,k_l1_offset,k_l2_offset,
     4      k_l3_offset,k_l4_offset)
          call reconcilefile(d_lr2,size_l2)
          call tce_filename('lr3',filename)
          call createfile(filename,d_lr3,size_l3)
          call ccsdtq_lambda3(d_f1,d_lr3,d_t1,d_t2,d_t3,d_v2,
     1      d_lambda1,d_lambda2,d_lambda3,d_lambda4,k_f1_offset,
     2      k_l3_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     3      k_v2_offset,k_l1_offset,k_l2_offset,
     4      k_l3_offset,k_l4_offset)
          call reconcilefile(d_lr3,size_l3)
          call tce_filename('lr4',filename)
          call createfile(filename,d_lr4,size_l4)
          call ccsdtq_lambda4(d_f1,d_lr4,d_t1,d_t2,d_v2,
     1      d_lambda2,d_lambda3,d_lambda4,k_f1_offset,
     2      k_l4_offset,k_t1_offset,k_t2_offset,
     3      k_v2_offset,k_l2_offset,k_l3_offset,k_l4_offset)
          call reconcilefile(d_lr4,size_l4)
          r1 = ddotfile(d_lr1,d_lr1,size_l1)
          r2 = ddotfile(d_lr2,d_lr2,size_l2)
          r3 = ddotfile(d_lr3,d_lr3,size_l3)
          r4 = ddotfile(d_lr4,d_lr4,size_l4)
          residual = max(r1,r2,r3,r4)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9410)
            endif
            if (ampnorms) then
              r1 = ddotfile(d_lambda1,d_lambda1,size_l1)
              r2 = ddotfile(d_lambda2,d_lambda2,size_l2)
              r3 = ddotfile(d_lambda3,d_lambda3,size_l3)
              r4 = ddotfile(d_lambda4,d_lambda4,size_l4)
              if (nodezero) then
                write(LuOut,9082) "Lambda singles",r1
                write(LuOut,9082) "Lambda doubles",r2
                write(LuOut,9082) "Lambda triples",r3
                write(LuOut,9082) "Lambda quadruples",r4
              endif
            endif
            call deletefile(d_lr4)
            call deletefile(d_lr3)
            call deletefile(d_lr2)
            call deletefile(d_lr1)
            call tce_diis_tidy()
            if (save_l(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving L1 now...'
              endif
              call y1_restart_save(d_lambda1,k_l1_offset,size_l1,0,
     1                             handle_l1,irrep_y)
            endif
            if (save_l(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving L2 now...'
              endif
              call y2_restart_save(d_lambda2,k_l2_offset,size_l2,0,
     1                             handle_l2,irrep_y)
            endif
            if (save_l(3)) then
              if(nodezero) then
                write(LuOut,*) 'Saving L3 now...'
              endif
              call y3_restart_save(d_lambda3,k_l3_offset,size_l3,0,
     1                             handle_l3,irrep_y)
            endif
            if (save_l(4)) then
              if(nodezero) then
                write(LuOut,*) 'Saving L4 now...'
              endif
              call y4_restart_save(d_lambda4,k_l4_offset,size_l4,0,
     1                             handle_l4,irrep_y)
            endif
            do axis=1,3
              irrep_d=sym_abelian_axis(geom,axis)
              call tce_filename('d0',filename)
              call createfile(filename,d_d0,size_d0)
              call ccsdtq_density1(d_d1(axis),d_d0,d_t1,d_t2,d_t3,d_t4,
     1          d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     2          k_d1_offset(axis),k_d0_offset,k_t1_offset,k_t2_offset,
     3          k_t3_offset,k_t4_offset,k_l1_offset,k_l2_offset,
     4          k_l3_offset,k_l4_offset)
              call reconcilefile(d_d0,size_d0)
              call get_block(d_d0,cdipole,1,0)
              dipole_cor(axis)=dipole_scf(axis)+cdipole
              call deletefile(d_d0)
            enddo
            dipole_cor_t = dsqrt( dipole_cor(1)*dipole_cor(1) +
     1                            dipole_cor(2)*dipole_cor(2) +
     2                            dipole_cor(3)*dipole_cor(3) )
            if (nodezero) write(LuOut,9430) "CCSDTQ",
     1                    dipole_cor(1),dipole_cor(1)*debye,
     2                    dipole_cor(2),dipole_cor(2)*debye,
     3                    dipole_cor(3),dipole_cor(3)*debye,
     4                    dipole_cor_t ,dipole_cor_t*debye
!            goto 100
            goto 2437
          endif
          if (save_l(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L1 now...'
            endif
            call y1_restart_save(d_lambda1,k_l1_offset,size_l1,0,
     1                           handle_l1,irrep_y)
          endif
          if (save_l(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L2 now...'
            endif
            call y2_restart_save(d_lambda2,k_l2_offset,size_l2,0,
     1                           handle_l2,irrep_y)
          endif
          if (save_l(3).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L3 now...'
            endif
            call y3_restart_save(d_lambda3,k_l3_offset,size_l3,0,
     1                           handle_l3,irrep_y)
          endif
          if (save_l(4).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L4 now...'
            endif
            call y4_restart_save(d_lambda4,k_l4_offset,size_l4,0,
     1                           handle_l4,irrep_y)
          endif
          call tce_diis(.true.,iter,.true.,.true.,.true.,.true.,
     1                  d_lr1,d_lambda1,k_l1_offset,size_l1,
     2                  d_lr2,d_lambda2,k_l2_offset,size_l2,
     3                  d_lr3,d_lambda3,k_l3_offset,size_l3,
     4                  d_lr4,d_lambda4,k_l4_offset,size_l4)
          call deletefile(d_lr4)
          call deletefile(d_lr3)
          call deletefile(d_lr2)
          call deletefile(d_lr1)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
        endif
! -------------
! CCSDTQ-LR
! -------------
 2437   continue
        if (lineresp) then
!
!          if (nodezero) write(LuOut,*) "lineresp is go"
          do omegacount=1,anumfreq
!            if (nodezero) write(LuOut,*) "omegacount is go"
            omega1 = afreq(omegacount)
            if (omega1.eq.(0.0d0)) then
              dynfreq = 1
            else
              dynfreq = 2
            endif
!
          do axis = 1, 3
!            if (nodezero) write(LuOut,*) "axis is go"
          if (respaxis(axis)) then
!            if (nodezero) write(LuOut,*) "respaxis(axis) is go"
            irrep_d=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            if (nodezero.and.util_print('mod1',print_default)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
            endif
            irrep_o=irrep_d
            irrep_x=irrep_d
!
          do omegasign=1,dynfreq
!            if (nodezero) write(LuOut,*) "omegasign is go"
            if (omegasign.eq.1) then
              omega1 = (1.0d0)*afreq(omegacount)
              dynaxis = 0
            elseif (omegasign.eq.2) then
              omega1 = (-1.0d0)*afreq(omegacount)
              dynaxis = 3
            endif
            if (nodezero) write(LuOut,9431) omega1
!            if (nodezero) write(LuOut,*) "omega1 = ",omega1
!            if (nodezero) write(LuOut,*) "dynaxis = ",dynaxis
!
            call tce_diis_init()
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSDTQ-LR"
!
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
              call tce_zero(d_rr1(axis),size_tr1(axis))
!        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega1,d_tr1(axis+dynaxis),
     1             d_rr1(axis),size_tr1(axis))
!        write(LuOut,*) "eomccsdtq_x1"
              call eomccsdtq_x1(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,
     3             k_tr1_offset(axis),
     4             k_t1_offset,k_t2_offset,k_v2_offset,
     5             k_tr1_offset(axis),k_tr2_offset(axis),
     6             k_tr3_offset(axis))
!        write(LuOut,*) "ccsdtq1_o1"
              call ccsdtq_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     &             k_tr1_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset)
!
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
              call tce_zero(d_rr2(axis),size_tr2(axis))
!        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega1,d_tr2(axis+dynaxis),
     1             d_rr2(axis),size_tr2(axis))
!        write(LuOut,*) "eomccsdtq_x2"
              call eomccsdtq_x2(d_f1,d_rr2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_tr4(axis+dynaxis),
     3             k_f1_offset,k_tr2_offset(axis),
     4             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     5             k_tr1_offset(axis),k_tr2_offset(axis),
     6             k_tr3_offset(axis),k_tr4_offset(axis))
!        write(LuOut,*) "ccsdtq_o2"
              call ccsdtq_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,d_t3,
     &             k_tr2_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset)
!
              call tce_filename(rr3filename(axis),filename)
              call createfile(filename,d_rr3(axis),size_tr3(axis))
              call tce_zero(d_rr3(axis),size_tr3(axis))
!        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega1,d_tr3(axis+dynaxis),
     1             d_rr3(axis),size_tr3(axis))
!        write(LuOut,*) "eomccsdtq_x3"
              call eomccsdtq_x3(d_f1,d_rr3(axis),d_t1,d_t2,d_t3,d_t4,
     1             d_v2,d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_tr4(axis+dynaxis),
     3             k_f1_offset,k_tr3_offset(axis),
     4             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5             k_v2_offset,
     6             k_tr1_offset(axis),k_tr2_offset(axis),
     7             k_tr3_offset(axis),k_tr4_offset(axis))
!        write(LuOut,*) "ccsdtq_o3"
              call ccsdtq_o3(d_rr3(axis),d_d1(axis),d_t1,d_t2,d_t3,d_t4,
     &             k_tr3_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset)
!
              call tce_filename(rr4filename(axis),filename)
              call createfile(filename,d_rr4(axis),size_tr4(axis))
              call tce_zero(d_rr4(axis),size_tr4(axis))
!        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega1,d_tr4(axis+dynaxis),
     1             d_rr4(axis),size_tr4(axis))
!        write(LuOut,*) "eomccsdtq_x4"
              call eomccsdtq_x4(d_f1,d_rr4(axis),d_t1,d_t2,d_t3,d_t4,
     1             d_v2,d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_tr4(axis+dynaxis),
     3             k_f1_offset,k_tr4_offset(axis),
     4             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5             k_v2_offset,
     6             k_tr1_offset(axis),k_tr2_offset(axis),
     7             k_tr3_offset(axis),k_tr4_offset(axis))
!        write(LuOut,*) "ccsdtq_o4"
              call ccsdtq_o4(d_rr4(axis),d_d1(axis),d_t1,d_t2,d_t3,d_t4,
     &             k_tr4_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset)
!
!        write(LuOut,*) "reconcilefile"
              call reconcilefile(d_rr1(axis),size_tr1(axis))
              call reconcilefile(d_rr2(axis),size_tr2(axis))
              call reconcilefile(d_rr3(axis),size_tr3(axis))
              call reconcilefile(d_rr4(axis),size_tr4(axis))
!        write(LuOut,*) "tce_residual_tr1"
              call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
!        write(LuOut,*) "tce_residual_tr2"
              call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
!        write(LuOut,*) "tce_residual_tr3"
              call tce_residual_tr3(d_rr3(axis),k_tr3_offset(axis),rr3)
!        write(LuOut,*) "tce_residual_tr4"
              call tce_residual_tr4(d_rr4(axis),k_tr4_offset(axis),rr4)
              residual = max(rr1,rr2,rr3,rr4)
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
              if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
              if (residual .lt. thresh) then
                if (nodezero) then
                  write(LuOut,9410)
                endif
                if (ampnorms) then
                  call tce_residual_tr1(d_tr1(axis+dynaxis),
     1                                  k_tr1_offset(axis),rr1)
                  call tce_residual_tr2(d_tr2(axis+dynaxis),
     1                                  k_tr2_offset(axis),rr2)
                  call tce_residual_tr3(d_tr3(axis+dynaxis),
     1                                  k_tr3_offset(axis),rr3)
                  call tce_residual_tr4(d_tr4(axis+dynaxis),
     1                                  k_tr4_offset(axis),rr4)
                  if (nodezero) then
                    write(LuOut,9082) "T(1) singles",rr1
                    write(LuOut,9082) "T(1) doubles",rr2
                    write(LuOut,9082) "T(1) triples",rr3
                    write(LuOut,9082) "T(1) quadruples",r4
                  endif
                endif

                call deletefile(d_rr4(axis))
                call deletefile(d_rr3(axis))
                call deletefile(d_rr2(axis))
                call deletefile(d_rr1(axis))
                call tce_diis_tidy()
                if (save_tr(1)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving T1(1) now...'
                  endif
                  call tr1_restart_save(d_tr1(axis+dynaxis),
     1                 k_tr1_offset(axis),size_tr1(axis),
     2                 axis+dynaxis,handle_tr1(axis),irrep_x)
                endif
                if (save_tr(2)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving T2(1) now...'
                  endif
                  call tr2_restart_save(d_tr2(axis+dynaxis),
     1                 k_tr2_offset(axis),size_tr2(axis),
     2                 axis+dynaxis,handle_tr2(axis),irrep_x)
                endif
                if (save_tr(3)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving T3(1) now...'
                  endif
                  call tr3_restart_save(d_tr3(axis),k_tr3_offset(axis),
     1                 size_tr3(axis),axis,handle_tr3(axis),irrep_x)
                endif
                if (save_tr(4)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving T4(1) now...'
                  endif
                  call tr4_restart_save(d_tr4(axis),k_tr4_offset(axis),
     1                 size_tr4(axis),axis,handle_tr4(axis),irrep_x)
                endif
                goto 2438
              endif
              if (save_tr(1).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving T1(1) now...'
                endif
                call tr1_restart_save(d_tr1(axis+dynaxis),
     1               k_tr1_offset(axis),size_tr1(axis),
     2               axis+dynaxis,handle_tr1(axis),irrep_x)
              endif
              if (save_tr(2).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving T2(1) now...'
                endif
                call tr2_restart_save(d_tr2(axis+dynaxis),
     1               k_tr2_offset(axis),size_tr2(axis),
     2               axis+dynaxis,handle_tr2(axis),irrep_x)
              endif
              if (save_tr(3).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving T3(1) now...'
                endif
                call tr3_restart_save(d_tr3(axis),k_tr3_offset(axis),
     1               size_tr3(axis),axis,handle_tr3(axis),irrep_x)
              endif
              if (save_tr(4).and.(mod(iter,save_interval).eq.0)) then
                if(nodezero) then
                  write(LuOut,*) 'Saving T4(1) now...'
                endif
                call tr4_restart_save(d_tr4(axis),k_tr4_offset(axis),
     1               size_tr4(axis),axis,handle_tr4(axis),irrep_x)
              endif
              call tce_diis2(.false.,iter,.true.,.true.,.true.,.true.,
     1             d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2             size_tr1(axis),
     3             d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4             size_tr2(axis),
     5             d_rr3(axis),d_tr3(axis+dynaxis),k_tr3_offset(axis),
     6             size_tr3(axis),
     7             d_rr4(axis),d_tr4(axis+dynaxis),k_tr4_offset(axis),
     8             size_tr4(axis),omega)
              call deletefile(d_rr4(axis))
              call deletefile(d_rr3(axis))
              call deletefile(d_rr2(axis))
              call deletefile(d_rr1(axis))
              if (nodezero) call util_flush(LuOut)
            enddo ! iter loop
          call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 2438 continue
!
          enddo ! omegasign loop
          endif ! respaxis(axis)
          enddo ! axis loop
!
!         Evaluate Dipole Polarizability
!
          do axisA = 1, 3
          do axisB = 1, axisA
          alpha(axisA,axisB)=0.0d0
          if (respaxis(axisA).and.respaxis(axisB)) then
            irrep_a=sym_abelian_axis(geom,axisA)
            irrep_b=sym_abelian_axis(geom,axisB)
            irrep_y=irrep_g
!
!            if (nodezero) write(LuOut,*) "axisA = ",axisA
!            if (nodezero) write(LuOut,*) "axisB = ",axisB
!
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,size_a0)
!
              alpha1=0.0d0
              alpha2=0.0d0
!
              if (omega1.ne.0.0d0) then
                dynaxis = 3
              else
                dynaxis = 0
              endif
!
              irrep_c=irrep_a
              irrep_oa=irrep_a
              irrep_ob=irrep_b
              irrep_tra=irrep_a
              irrep_trb=irrep_b
!              if (nodezero) write(LuOut,*) "ccsdtq_lr_alpha"
              call ccsdtq_lr_alpha(d_f1,d_a0,d_d1(axisA),d_d1(axisB),
     1             d_t1,d_t2,d_t3,d_t4,
     2             d_tr1(axisA+dynaxis),d_tr2(axisA+dynaxis),
     3             d_tr3(axisA+dynaxis),d_tr4(axisA+dynaxis),
     4             d_tr1(axisB),d_tr2(axisB),d_tr3(axisB),d_tr4(axisB),
     5             d_v2,d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     6             k_f1_offset,k_a0_offset,k_d1_offset(axisA),
     7             k_d1_offset(axisB),k_t1_offset,k_t2_offset,
     8             k_t3_offset,k_t4_offset,k_tr1_offset(axisA),
     9             k_tr2_offset(axisA),k_tr3_offset(axisA),
     1             k_tr4_offset(axisA),k_tr1_offset(axisB),
     2             k_tr2_offset(axisB),k_tr3_offset(axisB),
     3             k_tr4_offset(axisB),k_v2_offset,k_l1_offset,
     4             k_l2_offset,k_l3_offset,k_l4_offset)
!
              call reconcilefile(d_a0,size_a0)
              call get_block(d_a0,alpha1,1,0)
              call tce_zero(d_a0,size_a0)
!              write(LuOut,*) "alpha1 = ",alpha1
!
              if (omega1.eq.0.0d0) then
                alpha2=alpha1
              else
!                if (nodezero) write(LuOut,*) "ccsdtq_lr_alpha"
                call ccsdtq_lr_alpha(d_f1,d_a0,d_d1(axisA),d_d1(axisB),
     1             d_t1,d_t2,d_t3,d_t4,
     2             d_tr1(axisA),d_tr2(axisA),d_tr3(axisA),d_tr4(axisA),
     3             d_tr1(axisB+dynaxis),d_tr2(axisB+dynaxis),
     4             d_tr3(axisB+dynaxis),d_tr4(axisB+dynaxis),
     5             d_v2,d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     6             k_f1_offset,k_a0_offset,k_d1_offset(axisA),
     7             k_d1_offset(axisB),k_t1_offset,k_t2_offset,
     8             k_t3_offset,k_t4_offset,k_tr1_offset(axisA),
     9             k_tr2_offset(axisA),k_tr3_offset(axisA),
     1             k_tr4_offset(axisA),k_tr1_offset(axisB),
     2             k_tr2_offset(axisB),k_tr3_offset(axisB),
     3             k_tr4_offset(axisB),k_v2_offset,k_l1_offset,
     4             k_l2_offset,k_l3_offset,k_l4_offset)
!
                call reconcilefile(d_a0,size_a0)
                call get_block(d_a0,alpha2,1,0)
                call tce_zero(d_a0,size_a0)
              endif ! axisA.eq.axisB
!              write(LuOut,*) "alpha2 = ",alpha2
!
              alpha(axisA,axisB)=(-0.5d0)*(alpha1+alpha2)
!
            call deletefile(d_a0)
!
!            write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
!            if (nodezero) write(LuOut,9020) cpu, wall
!
          endif ! respaxis(axis)
          enddo ! axisB loop
          enddo ! axisA loop
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
!
          do i = 1, 3
            do j = 1, i
              alphacopy(i,j)=alpha(i,j)
              alphacopy(j,i)=alpha(i,j)
            enddo
          enddo
!
          call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
!
          alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
          alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
          alphaani = dsqrt(0.5d0*alphaani)
!
          if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1        .and.nodezero) write(LuOut,9911)
!
          if (nodezero) write(LuOut,9434) "CCSDTQ Linear Response",
     1      afreq(omegacount),auXnm/afreq(omegacount),
     2      alpha(1,1),alpha(2,1),alpha(3,1),
     3      au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4      alpha(2,1),alpha(2,2),alpha(3,2),
     5      au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6      alpha(3,1),alpha(3,2),alpha(3,3),
     7      au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8      alphaevr(1),alphaevr(2),alphaevr(3),
     9      au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1      alphaiso,au2ang3*alphaiso,
     2      alphaani,au2ang3*alphaani
          if (nodezero) write(LuOut,9020) cpu, wall
          call util_flush(LuOut)
!
          enddo ! omegacount loop
!
        endif ! lineresp
!
#endif
      endif
      if ((model.eq.'qcisd').or.do_qcisd) then
! -------------
! QCISD   right
! -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "QCISD"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call qcisd_e(d_e,d_t2,d_v2,k_e_offset,k_t2_offset,
     1                 k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_t1)
          call qcisd_t1(d_f1,d_r1,d_t1,d_t2,d_v2,
     1                  k_f1_offset,k_t1_offset,
     2                  k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_r1,size_t1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call qcisd_t2(d_f1,d_r2,d_t1,d_t2,d_v2,
     1                  k_f1_offset,k_t2_offset,
     2                  k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_r2,size_t2)
          call tce_residual_t1(d_r1,k_t1_offset,r1)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          residual = max(r1,r2)
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "QCISD",corr
              write(LuOut,9080) "QCISD",ref + corr
            endif
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            call tce_diis_tidy()
            if (ampnorms) then
              call tce_residual_t1(d_t1,k_t1_offset,r1)
              call tce_residual_t2(d_t2,k_t2_offset,r2)
              if (nodezero) then
                write(LuOut,9082) "C singles",r1
                write(LuOut,9082) "C doubles",r2
              endif
            endif
            call tce_print_x1(d_t1,k_t1_offset,printtol,irrep_t)
            call tce_print_x2(d_t2,k_t2_offset,printtol,irrep_t)
            if (save_t(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1 now...'
              endif
              call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                             handle_t1,irrep_t)
            endif
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
!            goto 100
            goto 3990
          endif
          if (save_t(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1 now...'
            endif
            call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                           handle_t1,irrep_t)
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          call tce_diis(.false.,iter,.true.,.true.,.false.,.false.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 3990 continue
! -------------------
! QCISD(T) and QCISD[T]
! -------------------
        if (do_pt.or.((perturbative.eq.'(t)').or.
     1                (perturbative.eq.'[t]'))) then
          cpu = - util_cpusec()
          wall = - util_wallsec()
          call ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                d_v2,k_v2_offset,pt3_1,pt3_2,size_t1)
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9075) "QCISD[T]",pt3_1
            write(LuOut,9070) "QCISD[T]",corr+pt3_1
            write(LuOut,9080) "QCISD[T]",ref+corr+pt3_1
            write(LuOut,9075) "QCISD(T)",pt3_2
            write(LuOut,9070) "QCISD(T)",corr+pt3_2
            write(LuOut,9080) "QCISD(T)",ref+corr+pt3_2
            write(LuOut,9020) cpu, wall
          endif
          if (perturbative.eq.'[t]') then
            corr = corr + pt3_1
          elseif (perturbative.eq.'(t)') then
            corr = corr + pt3_2
          endif
        endif
! -------------------------
! CR-QCISD(T)
! -------------------------
        if (do_cr_pt.or.(perturbative.eq.'cr_(t)')) then
          cpu = - util_cpusec()
          wall = - util_wallsec()
!          write(LuOut,*) 'tce_energy: cr_qcisd_t'
          call cr_qcisd_t(d_t1,k_t1_offset,size_t1,
     1                    d_t2,k_t2_offset,size_t2,
     2                    d_f1,k_f1_offset,d_v2,k_v2_offset,pt3_2)
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,9075) "CR-QCISD(T)",pt3_2
            write(LuOut,9070) "CR-QCISD(T)",corr+pt3_2
            write(LuOut,9080) "CR-QCISD(T)",ref+corr+pt3_2
            write(LuOut,9020) cpu, wall
          endif
          corr = corr + pt3_2
        endif
      endif
      if (model.eq.'cisd') then
! -------------
! CISD    right
! -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CISD"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call cisd_e(d_f1,d_e,d_t1,d_t2,d_v2,
     1                k_f1_offset,k_e_offset,
     2                k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_t1)
          call cisd_c1(d_e,d_f1,d_r1,d_t1,d_t2,d_v2,
     1                 k_e_offset,k_f1_offset,k_t1_offset,
     2                 k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_r1,size_t1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call cisd_c2(d_e,d_f1,d_r2,d_t1,d_t2,d_v2,
     1                 k_e_offset,k_f1_offset,k_t2_offset,
     2                 k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_r2,size_t2)
          call tce_residual_t1(d_r1,k_t1_offset,r1)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          residual = max(r1,r2)
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "CISD",corr
              write(LuOut,9080) "CISD",ref + corr
            endif
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            if (ampnorms) then
              call tce_residual_t1(d_t1,k_t1_offset,r1)
              call tce_residual_t2(d_t2,k_t2_offset,r2)
              if (nodezero) then
                write(LuOut,9082) "C singles",r1
                write(LuOut,9082) "C doubles",r2
              endif
            endif
            call tce_print_x1(d_t1,k_t1_offset,printtol,irrep_t)
            call tce_print_x2(d_t2,k_t2_offset,printtol,irrep_t)
            call tce_diis_tidy()
            if (save_t(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1 now...'
              endif
              call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                             handle_t1,irrep_t)
            endif
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
            goto 100
          endif
          if (save_t(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1 now...'
            endif
            call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                           handle_t1,irrep_t)
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          call tce_diis(.false.,iter,.true.,.true.,.false.,.false.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
      endif
      if (model.eq.'cisdt') then
! -------------
! CISDT   right
! -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CISDT"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call cisdt_e(d_f1,d_e,d_t1,d_t2,d_v2,
     1                 k_f1_offset,k_e_offset,
     2                 k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_t1)
          call cisdt_c1(d_e,d_f1,d_r1,d_t1,d_t2,d_t3,d_v2,
     1                  k_e_offset,k_f1_offset,k_t1_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset)
          call reconcilefile(d_r1,size_t1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call cisdt_c2(d_e,d_f1,d_r2,d_t1,d_t2,d_t3,d_v2,
     1                  k_e_offset,k_f1_offset,k_t2_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset)
          call reconcilefile(d_r2,size_t2)
          call tce_filename('r3',filename)
          call createfile(filename,d_r3,size_t3)
          call cisdt_c3(d_e,d_f1,d_r3,d_t1,d_t2,d_t3,d_v2,
     1                  k_e_offset,k_f1_offset,k_t3_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset)
          call reconcilefile(d_r3,size_t3)
          call tce_residual_t1(d_r1,k_t1_offset,r1)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          call tce_residual_t3(d_r3,k_t3_offset,r3)
          residual = max(r1,r2,r3)
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "CISDT",corr
              write(LuOut,9080) "CISDT",ref + corr
            endif
            call deletefile(d_r3)
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            if (ampnorms) then
              call tce_residual_t1(d_t1,k_t1_offset,r1)
              call tce_residual_t2(d_t2,k_t2_offset,r2)
              call tce_residual_t3(d_t3,k_t3_offset,r3)
              if (nodezero) then
                write(LuOut,9082) "C singles",r1
                write(LuOut,9082) "C doubles",r2
                write(LuOut,9082) "C triples",r3
              endif
            endif
            call tce_print_x1(d_t1,k_t1_offset,printtol,irrep_t)
            call tce_print_x2(d_t2,k_t2_offset,printtol,irrep_t)
            call tce_diis_tidy()
            if (save_t(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1 now...'
              endif
              call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                             handle_t1,irrep_t)
            endif
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
            if (save_t(3)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T3 now...'
              endif
              call x3_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                             handle_t3,irrep_t)
            endif
            goto 100
          endif
          if (save_t(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1 now...'
            endif
            call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                           handle_t1,irrep_t)
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          if (save_t(3).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T3 now...'
            endif
            call x3_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                           handle_t3,irrep_t)
          endif
          call tce_diis(.false.,iter,.true.,.true.,.true.,.false.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  d_r3,d_t3,k_t3_offset,size_t3,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r3)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
      endif
      if (model.eq.'cisdtq') then
! -------------
! CISDTQ  right
! -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CISDTQ"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call cisdtq_e(d_f1,d_e,d_t1,d_t2,d_v2,
     1                  k_f1_offset,k_e_offset,
     2                  k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_t1)
          call cisdtq_c1(d_e,d_f1,d_r1,d_t1,d_t2,d_t3,d_v2,
     1                 k_e_offset,k_f1_offset,k_t1_offset,k_t1_offset,
     2                 k_t2_offset,k_t3_offset,k_v2_offset)
          call reconcilefile(d_r1,size_t1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call cisdtq_c2(d_e,d_f1,d_r2,d_t1,d_t2,d_t3,d_t4,d_v2,
     1                 k_e_offset,k_f1_offset,k_t2_offset,k_t1_offset,
     2                 k_t2_offset,k_t3_offset,k_t4_offset,
     3                 k_v2_offset)
          call reconcilefile(d_r2,size_t2)
          call tce_filename('r3',filename)
          call createfile(filename,d_r3,size_t3)
          call cisdtq_c3(d_e,d_f1,d_r3,d_t1,d_t2,d_t3,d_t4,d_v2,
     1                 k_e_offset,k_f1_offset,k_t3_offset,k_t1_offset,
     2                 k_t2_offset,k_t3_offset,k_t4_offset,
     3                 k_v2_offset)
          call reconcilefile(d_r3,size_t3)
          call tce_filename('r4',filename)
          call createfile(filename,d_r4,size_t4)
          call cisdtq_c4(d_e,d_f1,d_r4,d_t2,d_t3,d_t4,d_v2,
     1                 k_e_offset,k_f1_offset,k_t4_offset,k_t2_offset,
     2                 k_t3_offset,k_t4_offset,k_v2_offset)
          call reconcilefile(d_r4,size_t4)
          call tce_residual_t1(d_r1,k_t1_offset,r1)
          call tce_residual_t2(d_r2,k_t2_offset,r2)
          call tce_residual_t3(d_r3,k_t3_offset,r3)
          call tce_residual_t4(d_r4,k_t4_offset,r4)
          residual = max(r1,r2,r3,r4)
          call get_block(d_e,corr,1,0)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "CISDTQ",corr
              write(LuOut,9080) "CISDTQ",ref + corr
            endif
            call deletefile(d_r4)
            call deletefile(d_r3)
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            call tce_print_x1(d_t1,k_t1_offset,printtol,irrep_t)
            call tce_print_x2(d_t2,k_t2_offset,printtol,irrep_t)
            call tce_diis_tidy()
            if (save_t(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1 now...'
              endif
              call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                             handle_t1,irrep_t)
            endif
            if (save_t(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2 now...'
              endif
              call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                             handle_t2,irrep_t)
            endif
            if (save_t(3)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T3 now...'
              endif
              call x3_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                             handle_t3,irrep_t)
            endif
            if (save_t(4)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T4 now...'
              endif
              call x4_restart_save(d_t4,k_t4_offset,size_t4,0,
     1                             handle_t4,irrep_t)
            endif
            goto 100
          endif
          if (save_t(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1 now...'
            endif
            call x1_restart_save(d_t1,k_t1_offset,size_t1,0,
     1                           handle_t1,irrep_t)
          endif
          if (save_t(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2 now...'
            endif
            call x2_restart_save(d_t2,k_t2_offset,size_t2,0,
     1                           handle_t2,irrep_t)
          endif
          if (save_t(3).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T3 now...'
            endif
            call x3_restart_save(d_t3,k_t3_offset,size_t3,0,
     1                           handle_t3,irrep_t)
          endif
          if (save_t(4).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T4 now...'
            endif
            call x4_restart_save(d_t4,k_t4_offset,size_t4,0,
     1                           handle_t4,irrep_t)
          endif
          call tce_diis(.false.,iter,.true.,.true.,.true.,.true.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  d_r3,d_t3,k_t3_offset,size_t3,
     4                  d_r4,d_t4,k_t4_offset,size_t4)
          call deletefile(d_r4)
          call deletefile(d_r3)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
      elseif ((model.eq.'mbpt2').or.
     1         (model.eq.'mbpt3').or.
     2         (model.eq.'mbpt4').or.
     3         (model.eq.'mbpt4sdq_t').or.
     3         (model.eq.'mbpt4sdq')) then
! -------------
! MBPT(2) right
! -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "MBPT(2)"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call mbpt2_e(d_e,d_u2,d_v2,
     1                 k_e_offset,k_u2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_r1)
          call mbpt1_t1(d_f1,d_r1,d_u1,d_u2,k_f1_offset,k_r1_offset,
     1                  k_u1_offset,k_u2_offset)
          call reconcilefile(d_r1,size_r1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_r2)
          call mbpt1_t2(d_f1,d_r2,d_u1,d_u2,d_v2,k_f1_offset,
     1                  k_r2_offset,k_u1_offset,k_u2_offset,
     2                  k_v2_offset)
          call reconcilefile(d_r2,size_r2)
          call tce_residual_t1(d_r1,k_r1_offset,r1)
          call tce_residual_t2(d_r2,k_r2_offset,r2)
          residual = max(r1,r2)
          call get_block(d_e,mbpt2,1,0)
          corr = mbpt2
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,mbpt2,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "MBPT(2)",mbpt2
              write(LuOut,9080) "MBPT(2)",ref + corr
              call util_flush(LuOut)
            endif
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            call tce_diis_tidy()
            if (model.eq.'mbpt2') goto 100
            goto 120
          endif
          call tce_diis(.false.,iter,.true.,.true.,.false.,.false.,
     1                  d_r1,d_u1,k_u1_offset,size_u1,
     2                  d_r2,d_u2,k_u2_offset,size_u2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
! -------------
! MBPT(3) right
! -------------
  120   continue
        call tce_zero(d_w1,size_w1)
        call tce_zero(d_w2,size_w2)
        if (model.eq.'mbpt4') then
          call tce_zero(d_w3,size_w3)
        endif
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "MBPT(3)"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call mbpt3_e(d_e,d_w2,d_v2,
     1                 k_e_offset,k_w2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_r1)
          call mbpt2_t1(d_f1,d_r1,d_u1,d_u2,d_w1,d_w2,d_v2,
     1                  k_f1_offset,k_r1_offset,k_u1_offset,
     2                  k_u2_offset,k_w1_offset,k_w2_offset,
     3                  k_v2_offset)
          call reconcilefile(d_r1,size_r1)
          if (model.eq.'mbpt3') then
            call tce_filename('r2',filename)
            call createfile(filename,d_r2,size_r2)
            call mbpt2_t2_b(d_f1,d_r2,d_u1,d_u2,d_w1,d_w2,d_v2,
     1                      k_f1_offset,k_r2_offset,k_u1_offset,
     2                      k_u2_offset,k_w1_offset,k_w2_offset,
     3                      k_v2_offset)
            call reconcilefile(d_r2,size_r2)
          elseif (model.eq.'mbpt4') then
            call tce_filename('r2',filename)
            call createfile(filename,d_r2,size_r2)
!->d3p975
            call mbpt2_t2(d_f1,d_r2,d_u1,d_u2,d_w2,d_w3,d_v2,
     1                    k_f1_offset,k_r2_offset,k_u1_offset,
     2                    k_u2_offset,k_w2_offset,
     3                    k_w3_offset,k_v2_offset)
!<-d3p975
            call reconcilefile(d_r2,size_r2)
            call tce_filename('r3',filename)
            call createfile(filename,d_r3,size_r3)
!->d3p975
            call mbpt2_t3(d_f1,d_r3,d_u1,d_u2,d_w2,d_w3,d_v2,
     1                    k_f1_offset,k_r3_offset,k_u1_offset,
     2                    k_u2_offset,k_w2_offset,k_w3_offset,
     3                    k_v2_offset)
!<-d3p975
            call reconcilefile(d_r3,size_r3)
!->d3p975
          elseif ((model.eq.'mbpt4sdq').or.
     1             (model.eq.'mbpt4sdq_t')) then
            call tce_filename('r2',filename)
            call createfile(filename,d_r2,size_r2)
            call mbpt2_t2_sd(d_f1,d_r2,d_u1,d_u2,d_w2,d_v2,
     1                       k_f1_offset,k_r2_offset,k_u1_offset,
     1                       k_u2_offset,k_w2_offset,k_v2_offset)
!
            call reconcilefile(d_r2,size_r2)
!<-d3p975
          endif
          call tce_residual_t1(d_r1,k_r1_offset,r1)
          call tce_residual_t2(d_r2,k_r2_offset,r2)
          if (model.eq.'mbpt4') then
            call tce_residual_t3(d_r3,k_r3_offset,r3)
          endif
          if ((model.eq.'mbpt3').or.(model.eq.'mbpt4sdq').or.
     1                              (model.eq.'mbpt4sdq_t')) then
            residual = max(r1,r2)
          else
            residual = max(r1,r2,r3)
          endif
          call get_block(d_e,mbpt3,1,0)
          corr = mbpt2 + mbpt3
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,mbpt3,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9060)
              write(LuOut,9070) "MBPT(3)",mbpt3
              write(LuOut,9080) "MBPT(3)",ref + corr
              call util_flush(LuOut)
            endif
            if (model.eq.'mbpt4') then
              call deletefile(d_r3)
            endif
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            call tce_diis_tidy()
            if (model.eq.'mbpt3') goto 100
            goto 130
          endif
          if (model.eq.'mbpt4') then
!->d3p975
            call tce_diis(.false.,iter,.true.,.true.,.true.,.false.,
     1                    d_r1,d_w1,k_w1_offset,size_w1,
     2                    d_r2,d_w2,k_w2_offset,size_w2,
     3                    d_r3,d_w3,k_w3_offset,size_w3,
     4                    dummy,dummy,dummy,dummy)
!<-d3p975
          else
            call tce_diis(.false.,iter,.true.,.true.,.false.,.false.,
     1                    d_r1,d_w1,k_w1_offset,size_w1,
     2                    d_r2,d_w2,k_w2_offset,size_w2,
     3                    dummy,dummy,dummy,dummy,
     4                    dummy,dummy,dummy,dummy)
          endif
          if (model.eq.'mbpt4') then
            call deletefile(d_r3)
          endif
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
! -------------
! MBPT(4) right
! -------------
  130   continue
        call put_block(d_e2,mbpt2,1,0)
        call tce_zero(d_t1,size_t1)
        call tce_zero(d_t2,size_t2)
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) then
             if(model.eq.'mbpt4') write(LuOut,9050) "MBPT(4)"
             if(model.eq.'mbpt4sdq') write(LuOut,9050) "MBPT(4)(SDQ)"
             if(model.eq.'mbpt4sdq_t') write(LuOut,9050) "MBPT(4)(SDQ)"
          endif
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call mbpt4_e(d_e,d_t2,d_v2,
     1                 k_e_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_r1)
          if(model.eq.'mbpt4')
     1     call mbpt3_t1(d_e2,d_f1,d_r1,d_u1,d_w1,d_w2,d_w3,d_t1,
     1                  d_t2,d_v2,k_e2_offset,k_f1_offset,k_r1_offset,
     2                  k_u1_offset,k_w1_offset,k_w2_offset,
     3                  k_w3_offset,k_t1_offset,k_t2_offset,k_v2_offset)
          if ((model.eq.'mbpt4sdq').or.(model.eq.'mbpt4sdq_t'))
     1     call mbpt3_t1_sd(d_e2,d_f1,d_r1,d_u1,d_w1,d_w2,d_t1,
     1                      d_t2,d_v2,k_e2_offset,k_f1_offset,
     1                      k_r1_offset,k_u1_offset,
     1                      k_w1_offset,k_w2_offset,
     1                      k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_r1,size_r1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_r2)
          if(model.eq.'mbpt4')
     1     call mbpt3_t2(d_e2,d_f1,d_r2,d_u2,d_w1,d_w2,d_w3,
     1                  d_t1,d_t2,d_v2,k_e2_offset,k_f1_offset,
     2                  k_r2_offset,k_u2_offset,k_w1_offset,
     3                  k_w2_offset,k_w3_offset,
     4                  k_t1_offset,k_t2_offset,k_v2_offset)
          if ((model.eq.'mbpt4sdq').or.(model.eq.'mbpt4sdq_t'))
     1      call mbpt3_t2_sd(d_e2,d_f1,d_r2,d_u2,d_w1,d_w2,
     1                       d_t1,d_t2,d_v2,k_e2_offset,k_f1_offset,
     1                       k_r2_offset,k_u2_offset,k_w1_offset,
     1                       k_w2_offset,k_t1_offset,k_t2_offset,
     1                       k_v2_offset)
          call reconcilefile(d_r2,size_r2)
          call tce_residual_t1(d_r1,k_r1_offset,r1)
          call tce_residual_t2(d_r2,k_r2_offset,r2)
          residual = max(r1,r2)
          call get_block(d_e,mbpt4,1,0)
          corr = mbpt2 + mbpt3 + mbpt4
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9100) iter,residual,mbpt4,cpu,wall
          call util_flush(LuOut)
          if (residual .lt. thresh) then
!
            if (model.eq.'mbpt4sdq_t') then
              cpu=-util_cpusec()
              wall=-util_wallsec()
              call ccsd_t(d_u1,k_u1_offset,d_u2,k_u2_offset,
     1                    d_v2,k_v2_offset,pt3_1,pt3_2,size_t1)
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
            endif
!
            if (nodezero) then
              if(model.eq.'mbpt4') then
                write(LuOut,9060)
                write(LuOut,9070) "MBPT(4)",mbpt4
                write(LuOut,9080) "MBPT(4)",ref + corr
              endif
!
              if (model.eq.'mbpt4sdq') then
                write(LuOut,9060)
                write(LuOut,9070) "MBPT(4)(SDQ)",mbpt4
                write(LuOut,9080) "MBPT(4)(SDQ)",ref + corr
              endif
!
              if (model.eq.'mbpt4sdq_t') then
                write(LuOut,9060)
                write(LuOut,9070) "MBPT(4)(SDQ)",mbpt4
                write(LuOut,9080) "MBPT(4)(SDQ)",ref + corr
                write(LuOut,9070) "MBPT(4)",mbpt4 + pt3_1
                write(LuOut,9080) "MBPT(4)",ref + corr + pt3_1
                write(LuOut,9480) 'triples',cpu, wall
              endif
              call util_flush(LuOut)
            endif
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            call tce_diis_tidy()
            goto 100
!            goto 2010
          endif
          call tce_diis(.false.,iter,.true.,.true.,.false.,.false.,
     1                  d_r1,d_t1,k_t1_offset,size_t1,
     2                  d_r2,d_t2,k_t2_offset,size_t2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
          call deletefile(d_r2)
          call deletefile(d_r1)
          call deletefile(d_e)
          if (nodezero) call util_flush(LuOut)
        enddo
        call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
      endif ! model
  100 continue
!
      if (nodezero) call util_flush(LuOut)
!
!     Delete files for energy calculation
!
      if (needw4) then
        call deletefile(d_w4)
        if (.not.ma_pop_stack(l_w4_offset))
     1    call errquit("tce_energy: MA problem",33,MA_ERR)
      endif
      if (needw3) then
        call deletefile(d_w3)
        if (.not.ma_pop_stack(l_w3_offset))
     1    call errquit("tce_energy: MA problem",32,MA_ERR)
      endif
      if (needw2) then
        call deletefile(d_w2)
        if (.not.ma_pop_stack(l_w2_offset))
     1    call errquit("tce_energy: MA problem",31,MA_ERR)
      endif
      if (needw1) then
        call deletefile(d_w1)
        if (.not.ma_pop_stack(l_w1_offset))
     1    call errquit("tce_energy: MA problem",30,MA_ERR)
      endif
      if (needu2) then
        call deletefile(d_u2)
        if (.not.ma_pop_stack(l_u2_offset))
     1    call errquit("tce_energy: MA problem",29,MA_ERR)
      endif
      if (needu1) then
        call deletefile(d_u1)
        if (.not.ma_pop_stack(l_u1_offset))
     1    call errquit("tce_energy: MA problem",34,MA_ERR)
      endif
      if (neede2) then
        call deletefile(d_e2)
        if (.not.ma_pop_stack(l_e2_offset))
     1    call errquit("tce_energy: MA problem",35,MA_ERR)
      endif
      if (do_mbpt) then
        if (needt4) then
          if (.not.ma_pop_stack(l_r4_offset))
     1      call errquit("tce_energy: MA problem",27,MA_ERR)
        endif
        if (needt3) then
          if (.not.ma_pop_stack(l_r3_offset))
     1      call errquit("tce_energy: MA problem",21,MA_ERR)
        endif
        if (needt3a) then
          if (.not.ma_pop_stack(l_r3_offset))
     1      call errquit("tce_energy: MA problem",21,MA_ERR)
        endif
        if (needt2) then
          if (.not.ma_pop_stack(l_r2_offset))
     1      call errquit("tce_energy: MA problem",22,MA_ERR)
        endif
        if (needt1) then
          if (.not.ma_pop_stack(l_r1_offset))
     1      call errquit("tce_energy: MA problem",23,MA_ERR)
        endif
      endif
      if (quadresp) then
        if (.not.ma_pop_stack(l_g0_offset))
     1    call errquit("tce_energy: MA problem",24,MA_ERR)
      endif
      if (leftresp) then
        if (.not.ma_pop_stack(l_b0_offset))
     1    call errquit("tce_energy: MA problem",25,MA_ERR)
      endif
      if (lineresp) then
        if (.not.ma_pop_stack(l_a0_offset))
     1    call errquit("tce_energy: MA problem",26,MA_ERR)
      endif
      if (left) then
        if (multipole.gt.2) then
          if (.not.ma_pop_stack(l_o0_offset))
     1      call errquit("tce_energy: MA problem",300,MA_ERR)
        endif
        if (multipole.gt.1) then
          if (.not.ma_pop_stack(l_q0_offset))
     1      call errquit("tce_energy: MA problem",200,MA_ERR)
        endif
        if (.not.ma_pop_stack(l_d0_offset))
     1    call errquit("tce_energy: MA problem",100,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_e_offset))
     1  call errquit("tce_energy: MA problem",24,MA_ERR)
      if ((nroots.eq.0).or.(.not.excitedfragment)) then
        if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,ref+corr))
     1  call errquit('tce_energy: RTDB problem',0,MA_ERR)
! link to QM/MM
        if (.not.rtdb_put(rtdb,'tce:enegr',mt_dbl,1,ref+corr))
     1  call errquit('tce_energy: RTDB problem',0,MA_ERR)
!
        tce_energy=.true.
      else
!
!     ============================
!     Excited-state CI, CC, & MBPT
!     ============================
!
        irrep_g = 0

!#ifdef EACCSD
      if (.not.rtdb_get(rtdb,'tce:eaccsd',mt_log,1,eaccsd))
     +  eaccsd = .false.
      if(eaccsd) then
       ea_unused_spin=1
       ea_unused_sym=0
       if(nodezero) write(LuOut,*)
       if(nodezero) write(LuOut,*) "EACCSD calculation"
       if(nodezero) call util_flush(LuOut)
      endif
!#endif
!#ifdef IPCCSD
      if (.not.rtdb_get(rtdb,'tce:ipccsd',mt_log,1,ipccsd))
     +  ipccsd = .false.
      if(ipccsd) then
       ip_unused_spin=1
       ip_unused_sym=0
       if(nodezero) write(LuOut,*)
       if(nodezero) write(LuOut,*) "IPCCSD calculation"
       if(nodezero) call util_flush(LuOut)
      endif
!#endif
!
!
           if (ipccsd.and.(eom_solver.eq.2)) eom_solver=1
           if (eaccsd.and.(eom_solver.eq.2)) eom_solver=1
           if (ipccsd.and.(ccsd_var.eq.'ic')) ccsd_var='xx'
           if (eaccsd.and.(ccsd_var.eq.'ic')) ccsd_var='xx'
!
           if((ipccsd.or.eaccsd).and.(.not.restricted))
     1      call errquit("EA-IP-CCSD requires RHF reference",24,MA_ERR) 
!
!
!
!       Restricted closed shell: irrep_g is always zero.
!
        if (ipol.eq.2) then
          do i = 1,2
            do j = 1,nocc(i)
              irrep_g = ieor(irrep_g,int_mb(k_irs(i)+j-1))
            enddo
          enddo
        endif
        call sym_irrepname(geom,irrep_g+1,irrepname)
        if (nodezero.and.util_print('eom',print_default)) then
          write(LuOut,*)
          write(LuOut,9250) irrepname
        endif
        do irrep = 0, nirreps-1
          irrep_x = irrep
          irrep_y = irrep
          call sym_irrepname(geom,ieor(irrep_x,irrep_g)+1,irrepname)
          if ((.not.symmetry).or.(targetsym.eq.irrepname)) then
          call tce_eom_init
          if (nodezero.and.util_print('eom',print_default)) then
            write(LuOut,*)
            write(LuOut,9200) irrepname
          endif
!
          if(nodezero) then
           write(LuOut,9251) hbard
          end if
!
!cx          if(ccsd_var.eq.'ic') eom_solver=2
!
! new EOMCC solver
       if(.not.((model.ne.'ccsd').and.(eom_solver.eq.2))) then
! ccsd_act/eomccsd_act
       if(.not.(model.eq.'ccsd_act')) then
       if (.not.ma_push_get(mt_dbl,hbard*hbard,'hbar',
     1  l_hbar,k_hbar)) call errquit('tce_eom_xdiagon: MA problem',0,
     2  MA_ERR)
       call tce_hbarinit(dbl_mb(k_hbar),hbard)
        end if
        end if
!
!
          if (.not.ma_push_get(mt_dbl,maxtrials,'omegax',
     1      l_omegax,k_omegax))
     2      call errquit('tce_energy: MA problem',1000,MA_ERR)
          if (left) then
            if (.not.ma_push_get(mt_dbl,maxtrials,'omegay',
     1        l_omegay,k_omegay))
     2        call errquit('tce_energy: MA problem',1001,MA_ERR)
            call tce_e_offset(l_d0_offset,k_d0_offset,size_d0)
            call tce_e_offset(l_x0_offset,k_x0_offset,size_x0)
            call tce_e_offset(l_y0_offset,k_y0_offset,size_y0)
          endif
!
          if (needt1) then
           if((.not. eaccsd) .and. (.not. ipccsd)) then
            call tce_x1_offset(l_x1_offset,k_x1_offset,size_x1)
#ifdef IPCCSD
           else if( ipccsd) then
            call tce_ipx1_offset(l_x1_offset,k_x1_offset,size_x1)
#endif
#ifdef EACCSD
           else if (eaccsd) then
            call tce_eax1_offset(l_x1_offset,k_x1_offset,size_x1)
#endif
           endif
            call tce_filename('rx1',filename)
            call createfile(filename,d_rx1,size_x1)
            call gatoeaf(d_rx1)
            if (left) then
              call tce_y1_offset(l_y1_offset,k_y1_offset,size_y1)
              call tce_filename('ry1',filename)
              call createfile(filename,d_ry1,size_y1)
              call gatoeaf(d_ry1)
            endif
            if (nodezero.and.util_print('x1',print_default)) then
              write(LuOut,*)
              write(LuOut,9000) 'x1',size_x1
              call util_flush(LuOut)
            endif
          endif
          if (needt2) then
           if(model.ne.'ccsd_act') then
            if((.not. eaccsd) .and. (.not. ipccsd)) then
             call tce_x2_offset(l_x2_offset,k_x2_offset,size_x2)
#ifdef IPCCSD
            else if(ipccsd) then
             call tce_ipx2_offset(l_x2_offset,k_x2_offset,size_x2)
#endif
#ifdef EACCSD
            else if(eaccsd) then
             call tce_eax2_offset(l_x2_offset,k_x2_offset,size_x2)
#endif
            endif
           else 
! ccsd_act/eomccsd_act
            call tce_x2a_offset(l_x2_offset,k_x2_offset,size_x2)
           end if 
            call tce_filename('rx2',filename)
            call createfile(filename,d_rx2,size_x2)
            call gatoeaf(d_rx2)
            if (left) then
              call tce_y2_offset(l_y2_offset,k_y2_offset,size_y2)
              call tce_filename('ry2',filename)
              call createfile(filename,d_ry2,size_y2)
              call gatoeaf(d_ry2)
            endif
            if (nodezero.and.util_print('x2',print_default)) then
              write(LuOut,*)
              write(LuOut,9000) 'x2',size_x2
              call util_flush(LuOut)
            endif
          endif
          if (needt3a) then
            call tce_x3a_offset(l_x3_offset,k_x3_offset,size_x3)
            call tce_filename('rx3',filename)
            call createfile(filename,d_rx3,size_x3)
            call gatoeaf(d_rx3)
!ccx            if (left) then
!ccx              call tce_y3_offset(l_y3_offset,k_y3_offset,size_y3)
!ccx              call tce_filename('ry3',filename)
!ccx              call createfile(filename,d_ry3,size_y3)
!ccx              call gatoeaf(d_ry3)
!ccx            endif
            if (nodezero.and.util_print('x3',print_default)) then
              write(LuOut,*)
              write(LuOut,9000) 'x3',size_x3
              call util_flush(LuOut)
            endif
          endif
          if (needt3) then
            call tce_x3_offset(l_x3_offset,k_x3_offset,size_x3)
            call tce_filename('rx3',filename)
            call createfile(filename,d_rx3,size_x3)
            call gatoeaf(d_rx3)
            if (left) then
              call tce_y3_offset(l_y3_offset,k_y3_offset,size_y3)
              call tce_filename('ry3',filename)
              call createfile(filename,d_ry3,size_y3)
              call gatoeaf(d_ry3)
            endif
            if (nodezero.and.util_print('x3',print_default)) then
              write(LuOut,*)
              write(LuOut,9000) 'x3',size_x3
              call util_flush(LuOut)
            endif
          endif
          if (needt4) then
            call tce_x4_offset(l_x4_offset,k_x4_offset,size_x4)
            call tce_filename('rx4',filename)
            call createfile(filename,d_rx4,size_x4)
            call gatoeaf(d_rx4)
            if (left) then
              call tce_y4_offset(l_y4_offset,k_y4_offset,size_y4)
              call tce_filename('ry4',filename)
              call createfile(filename,d_ry4,size_y4)
              call gatoeaf(d_ry4)
            endif
            if (nodezero.and.util_print('x4',print_default)) then
              write(LuOut,*)
              write(LuOut,9000) 'x4',size_x4
              call util_flush(LuOut)
            endif
          endif
!
!         ------------------------------
!         Generate initial trial vectors
!         ------------------------------
!
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (model.eq.'ccsd') then


      if(eom_solver.eq.1) then !  ------- EOM-SOLVER ----- (solver 1)
!
!
      if(nodezero) then
       write(LuOut,9510) eom_solver
       call util_flush(LuOut)
      end if
!
      if(ccsd_var.ne.'cc2') then
      ister=0
! rec_s_1
       CALL OFFSET_eomccsd_x1_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs1',-1,MA_ERR)
! rec_s_2
       CALL OFFSET_eomccsd_x1_2_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_2_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs2',-1,MA_ERR)
! rec_s_3
       CALL OFFSET_eomccsd_x1_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs3',-1,MA_ERR)
! rec_s_4
       CALL OFFSET_eomccsd_x1_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs4',-1,MA_ERR)
! rec_s_5
       CALL OFFSET_eomccsd_x1_7_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_7_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_s_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs5',-1,MA_ERR)
! rec_d_1
       CALL OFFSET_eomccsd_x2_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd1',-1,MA_ERR)
! rec_d_2
       CALL OFFSET_eomccsd_x2_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_3_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd2',-1,MA_ERR)
! rec_d_3
       CALL OFFSET_eomccsd_x2_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd3',-1,MA_ERR)
! rec_d_4
       CALL OFFSET_eomccsd_x2_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd4',-1,MA_ERR)
! rec_d_5
       CALL OFFSET_eomccsd_x2_6_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_6_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd5',-1,MA_ERR)
! rec_d_6
       CALL OFFSET_eomccsd_x2_8_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_6,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd6',-1,MA_ERR)
! rec_d_7
       CALL OFFSET_eomccsd_x2_8_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_3_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_7,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd7',-1,MA_ERR)
! rec_d_8
       CALL OFFSET_eomccsd_x2_8_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_4_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_8,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd8',-1,MA_ERR)
! rec_d_9
      CALL OFFSET_eomccsd_x2_10_1_1(l_aux_off,k_aux_off,size_aux)
      CALL TCE_FILENAME('eomccsd_x2_10_1_1_i2',filename)
      CALL CREATEFILE(filename,rec_d_9,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd9',-1,MA_ERR)
      end if 
!
      if(ccsd_var.eq.'lr-ccsd') then
            call tce_x1_offset(l_j1_offset,k_j1_offset,size_j1)
            call tce_filename('j1',filename)
            call createfile(filename,d_j1,size_j1)
!c            call gatoeaf(d_j1)
            call tce_zero(d_j1,size_j1)
!
            call tce_x2_offset(l_j2_offset,k_j2_offset,size_j2)
            call tce_filename('j2',filename)
            call createfile(filename,d_j2,size_j2)
!c            call gatoeaf(d_j2)
            call tce_zero(d_j2,size_j2)
!
            call tce_loc_j1_one(d_j1,k_j1_offset)
            call tce_loc_j2_one(d_j2,k_j2_offset)
!
            call tce_loc_j1(d_j1,k_j1_offset,d_t1,k_t1_offset)
            call tce_loc_j2(d_j2,k_j2_offset,d_t2,k_t2_offset)
!
            call tce_loc_j1_inv(d_j1,k_j1_offset)
            call tce_loc_j2_inv(d_j2,k_j2_offset)
      end if
!
           if((.not. eaccsd) .and. (.not. ipccsd)) then
            call tce_eom_xguess(rtdb,needt1,needt2,needt3,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
#ifdef IPCCSD
           else if(ipccsd) then 
!             if(nodezero) write(LuOut,*) "I am here 2"
!             if(nodezero) call util_flush(LuOut)
            call tce_eom_ipxguess(rtdb,needt1,needt2,needt3,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
#endif
#ifdef EACCSD
           else if(eaccsd) then
!kbn              write(LuOut,*) "I am here 1.2"
!kbn              call util_flush(LuOut)
            call tce_eom_eaxguess(rtdb,needt1,needt2,needt3,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
!kbn              write(LuOut,*) "I am here 1.3"
!kbn              call util_flush(LuOut)
#endif
           endif


            if (nxtrials.eq.0) goto 200
!cxsg
          do ivec=1,nroots_reduced
            call tce_filenameindexed(ivec,'xc1',filename)
            call createfile(filename,xc1(ivec),size_x1)
            call gatoeaf(xc1(ivec))
            call gatodra(xc1(ivec))
            xc1_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc2',filename)
            call createfile(filename,xc2(ivec),size_x2)
            call gatoeaf(xc2(ivec))
            call gatodra(xc2(ivec))
            xc2_exist(ivec) = .true.
          enddo
!cxsg
          if(ccsd_var.eq.'lr-ccsd') then
            modelname = "LR-EOMCCSD right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
          else if(ccsd_var.eq.'cc2') then
            modelname = "CC2 right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
          else
            modelname = "EOM-CCSD right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
          end if
!
!           call hbar_init
            do iter=1,maxiter
              if (nodezero.and.util_print('eom',print_default))
     1          write(LuOut,9210) iter,nxtrials
              do ivec = 1,nxtrials
                if (.not.xp1_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp1',filename)
                  call createfile(filename,xp1(ivec),size_x1)
                  xp1_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
#if defined(EACCSD) || defined(IPCCSD)
        if(eaccsd .or. ipccsd) then
!kbn         if((nodezero) .and. (iter .eq. 1)) 
!kbn     +    write(LuOut,*) "No local T1 and X1 for EACCSD and IPCCSD"
!kbn         if((nodezero) .and. (iter .eq. 1) ) call util_flush(LuOut) 
        else
#endif
! - T1/X1 LOCALIZATION ----------
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,size_x1,'x1_local',
     1      l_x1_local,k_x1_local))
     1      call errquit('x1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
        call ma_zero(dbl_mb(k_x1_local),size_x1)
!    copy d_t1 ==> l_t1_local
!    copy x1(ivec) ==> l_x1_local
!c        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
!c        call ga_get(x1(ivec),1,size_x1,1,1,dbl_mb(k_x1_local),1)
      call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
      call get_block(x1(ivec),dbl_mb(k_x1_local),size_x1,0)
!cx       end if
! -------------------------------
#if defined(EACCSD) || defined(IPCCSD)
        endif ! if(eaccsd .or. ipccsd)
#endif
         if(ccsd_var.eq.'cc2') then
          call cc2_x1(d_f1,xp1(ivec),k_t1_local,d_t2,d_v2,
     &      k_x1_local,x2(ivec),
     &      k_f1_offset,k_x1_offset,k_t1_offset,k_t2_offset,
     &      k_v2_offset,k_x1_offset,k_x2_offset)
         else
          if(ccsd_var.eq.'ic') then
           call eomicsd_x1(d_f1,xp1(ivec),k_t1_local,d_t2,d_v2,
     1      k_x1_local,x2(ivec),k_f1_offset,k_x1_offset,
     2      k_t1_offset,k_t2_offset,k_v2_offset,
     3      k_x1_offset,k_x2_offset,
     4      ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5      rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6      rec_d_6,rec_d_7,rec_d_8,rec_d_9)
          else
#if defined(EACCSD) || defined(IPCCSD)
           if(eaccsd) then
#endif
#ifdef EACCSD
!cx             if(nodezero) write(LuOut,*) "I am here eaccsd_x1"
!cx             if(nodezero) call util_flush(LuOut)
!kbn      SUBROUTINE eaccsd_x1(d_f1,d_i0,d_t1,d_t2,d_v2,d_x1,d_x2,k_f1_offse
!kbn     &t,k_i0_offset,k_t1_offset,k_t2_offset,k_v2_offset,k_x1_offset,k_x2
!kbn     &_offset)
             call eaccsd_x1(d_f1,xp1(ivec),d_t1,d_t2,d_v2,x1(ivec),
     +                      x2(ivec),
     +                      k_f1_offset,k_x1_offset,k_t1_offset,
     +                      k_t2_offset,k_v2_offset,k_x1_offset,
     +                      k_x2_offset)
!cx             if(nodezero) write(LuOut,*) "I am here eaccsd_x1 1"
!cx             if(nodezero) call util_flush(LuOut)
#endif
#if defined(EACCSD) || defined(IPCCSD)
           elseif(ipccsd) then
#endif
!kbn             if(nodezero) write(LuOut,*) "I am here ipccsd_x1"
!kbn             if(nodezero) call util_flush(LuOut)
!kbn      SUBROUTINE ipccsd_x1(d_f1,d_i0,d_t1,d_t2,d_v2,d_x1,d_x2,k_f1_offse
!kbn     &t,k_i0_offset,k_t1_offset,k_t2_offset,k_v2_offset,k_x1_offset,k_x2
!kbn     &_offset)
#ifdef IPCCSD
!cx             if(nodezero) write(LuOut,*) "I am here ipccsd_x1"
!cx             if(nodezero) call util_flush(LuOut)
             call ipccsd_x1(d_f1,xp1(ivec),d_t1,d_t2,d_v2,x1(ivec),
     +                      x2(ivec),
     +                      k_f1_offset,k_x1_offset,k_t1_offset,
     +                      k_t2_offset,k_v2_offset,k_x1_offset,
     +                      k_x2_offset)
!cx             if(nodezero) write(LuOut,*) "I am here ipccsd_x1 1"
!cx             if(nodezero) call util_flush(LuOut)
#endif
#if defined(EACCSD) || defined(IPCCSD)
           else
#endif
!             if(nodezero) write(LuOut,*) "I am here eomccsd x1"
!             if(nodezero) call util_flush(LuOut)
                  call eomccsd_x1_n(d_f1,xp1(ivec),k_t1_local,d_t2,d_v2,
     1              k_x1_local,x2(ivec),k_f1_offset,k_x1_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,
     3              k_x1_offset,k_x2_offset,
     4 ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9)
#if defined(EACCSD) || defined(IPCCSD)
           endif
#endif
          end if
         end if
                  call reconcilefile(xp1(ivec),size_x1)
           if(ccsd_var.eq.'lr-ccsd') then
              call tce_j1_x1_scaling(xp1(ivec),k_x1_offset,
     1             d_j1,k_j1_offset)
           end if
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp1(ivec))
                endif
                if (.not.xp2_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp2',filename)
                  call createfile(filename,xp2(ivec),size_x2)
                  xp2_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
         if(ccsd_var.eq.'cc2') then
           call cc2_x2(d_f1,xp2(ivec),k_t1_local,d_v2,
     &       k_x1_local,x2(ivec),
     &       k_f1_offset,k_x2_offset,k_t1_offset,k_v2_offset,
     &       k_x1_offset,k_x2_offset,size_x1,size_x2)
         else
          if(ccsd_var.eq.'ic') then
           irrep_c=irrep_x
           call tce_filename('c2',filename)
           call createfile(filename,d_c2,size_x2)
           call tce_zero(d_c2,size_x2)
           call c2excit2_loc(d_c2,k_t1_local,k_x1_local,
     &       x2(ivec),k_x2_offset,k_t1_offset,k_x1_offset,
     &       k_x2_offset)
           call eomicsd_x2(d_f1,xp2(ivec),k_t1_local,d_t2,d_v2,
     1       k_x1_local,x2(ivec),k_f1_offset,k_x2_offset,
     2       k_t1_offset,k_t2_offset,k_v2_offset,
     3       k_x1_offset,k_x2_offset,size_x1,size_x2,
     4       ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9,d_c2)
           call deletefile(d_c2)
           else
#if defined(EACCSD) || defined(IPCCSD)
            if(eaccsd) then
#endif
#ifdef EACCSD
!cx             if(nodezero) write(LuOut,*) "I am here eaccsd_x2"
!cx             if(nodezero) call util_flush(LuOut)
!kbn      SUBROUTINE eaccsd_x2(d_f1,d_i0,d_t1,d_t2,d_v2,d_x1,d_x2,k_f1_offse
!kbn     &t,k_i0_offset,k_t1_offset,k_t2_offset,k_v2_offset,k_x1_offset,k_x2
!kbn     &_offset)
          call eaccsd_x2(d_f1,xp2(ivec),d_t1,d_t2,d_v2,x1(ivec),
     +                      x2(ivec),
     +                      k_f1_offset,k_x2_offset,k_t1_offset,
     +                      k_t2_offset,k_v2_offset,k_x1_offset,
     +                      k_x2_offset)
!cx             if(nodezero) write(LuOut,*) "I am here eaccsd_x2 2"
!cx             if(nodezero) call util_flush(LuOut)
#endif
#ifdef EACCSD
            elseif(ipccsd) then
!kbn             if(nodezero) write(LuOut,*) "I am here ipccsd_x2"
!kbn             if(nodezero) call util_flush(LuOut)
!kbn      SUBROUTINE ipccsd_x2(d_f1,d_i0,d_t1,d_t2,d_v2,d_x1,d_x2,k_f1_offse
!kbn     &t,k_i0_offset,k_t1_offset,k_t2_offset,k_v2_offset,k_x1_offset,k_x2
!kbn     &_offset)
!cx             if(nodezero) write(LuOut,*) "I am here ipccsd_x2"
!cx             if(nodezero) call util_flush(LuOut)
                  call ipccsd_x2(d_f1,xp2(ivec),d_t1,d_t2,d_v2,x1(ivec),
     +                           x2(ivec),
     +                           k_f1_offset,k_x2_offset,k_t1_offset,
     +                           k_t2_offset,k_v2_offset,k_x1_offset,
     +                           k_x2_offset)
!cx             if(nodezero) write(LuOut,*) "I am here ipccsd_x2 2"
!cx             if(nodezero) call util_flush(LuOut)
#endif
#if defined(EACCSD) || defined(IPCCSD)
            else
#endif
                  call eomccsd_x2_n(d_f1,xp2(ivec),k_t1_local,d_t2,d_v2,
     1              k_x1_local,x2(ivec),k_f1_offset,k_x2_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,
     3              k_x1_offset,k_x2_offset,size_x1,size_x2,
     4     ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5             rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6             rec_d_6,rec_d_7,rec_d_8,rec_d_9)
#if defined(EACCSD) || defined(IPCCSD)
            endif
#endif
           end if
           ister=1
!
         end if
                  call reconcilefile(xp2(ivec),size_x2)
#if defined(EACCSD) || defined(IPCCSD)
        if(eaccsd .or. ipccsd) then
!kbn         if(nodezero) 
!kbn     +    write(LuOut,*) "No local T1 and X1 for EACCSD and IPCCSD"
!kbn         if(nodezero) call util_flush(LuOut) 
        else
#endif
! - T1/X1 LOCALIZATION --
         if(.not.MA_POP_STACK(l_x1_local))
     &      call errquit('l_x1_local',4,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
! -----------------------
#if defined(EACCSD) || defined(IPCCSD)
        endif
#endif

           if(ccsd_var.eq.'lr-ccsd') then
                  call tce_j2_x2_scaling(xp2(ivec),k_x2_offset,
     1              d_j2,k_j2_offset)
           end if
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp2(ivec))
                endif
              enddo
              if (.not.ma_push_get(mt_dbl,nxtrials,'residual',
     1          l_residual,k_residual))
     2          call errquit('tce_energy: MA problem',101,MA_ERR)
#if defined(EACCSD) || defined(IPCCSD)
             call tce_eom_xdiagon(needt1,needt2,needt3,needt4,
     1          size_x1,size_x2,size_x3,size_x4,
     2          k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
     3          d_rx1,d_rx2,d_rx3,d_rx4,
     4          dbl_mb(k_omegax),dbl_mb(k_residual),k_hbar,iter,
     +          eaccsd,ipccsd)
#else
              call tce_eom_xdiagon(needt1,needt2,needt3,needt4,
     1          size_x1,size_x2,size_x3,size_x4,
     2          k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
     3          d_rx1,d_rx2,d_rx3,d_rx4,
     4          dbl_mb(k_omegax),dbl_mb(k_residual),k_hbar,iter,
     +          eaccsd,ipccsd)
#endif
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
              converged = .true.
              do ivec = 1,nroots_reduced
                if (nodezero.and.(ivec.ne.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegax+ivec-1),
     3            dbl_mb(k_omegax+ivec-1)*au2ev
                if (nodezero.and.(ivec.eq.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegax+ivec-1),
     3            dbl_mb(k_omegax+ivec-1)*au2ev,cpu,wall
                if (nodezero) call util_flush(LuOut)
                if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1            converged = .false.
              enddo
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (.not.ma_pop_stack(l_residual))
     1          call errquit("tce_energy: MA problem",102,MA_ERR)
              if (converged) then
                call tce_eom_xtidy
!               call hbar_tidy
                if (save_t(1)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving X1 now...'
                  endif
                  do ivec = 1,min(nroots_reduced,8)
                    call x1_restart_save(x1(ivec),k_x1_offset,size_x1,
     1                             ivec,handle_x1,irrep_x)
                  enddo
                endif
                if (save_t(2)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving X2 now...'
                  endif
                  do ivec = 1,min(nroots_reduced,8)
                    call x2_restart_save(x2(ivec),k_x2_offset,size_x2,
     1                             ivec,handle_x2,irrep_x)
                  enddo
                endif
                if (nodezero) then
                  write(LuOut,9240)
                  call util_flush(LuOut)
                endif
! ---
cx            if(nodezero) then
cx             write(LuOut,*)'largest EOMCCSD amplitudes: R1 and R2'
cx             call util_flush(LuOut)
cx            end if
             do jvec=1,nroots_reduced
c right EOMCCSD  printout
             if (nodezero) then
               write(LuOut,*)'------------------'
               write(LuOut,9460)
     1           jvec,dbl_mb(k_omegax+jvec-1),
     2           dbl_mb(k_omegax+jvec-1)*au2ev
             endif
c
            if(nodezero) then
             write(LuOut,*)'largest EOMCCSD amplitudes: R1 and R2'
             call util_flush(LuOut)
            end if
c
              if((.not. eaccsd) .and. (.not. ipccsd)) then
               call tce_print_x1(xc1(jvec),k_x1_offset,0.10d0,irrep_x)
               call tce_print_x2(xc2(jvec),k_x2_offset,0.10d0,irrep_x)
!
!QuantumComputing
          if(quantumc) then
          if(nodezero) then
           write(6,*)'Excited state specification:'
           write(6,*)'Rexp(T)|Ref>~(R_1+R_2)|Ref>'
           call util_flush(6)
          endif 
          if(nodezero) write(6,*)'------- a(p)+/a(q) strings ----------'
           xnorm1=0.0d0
           xnorm2=0.0d0
           xnorm=0.0d0
c
           call tce_string_x1_norm(xc1(jvec),k_x1_offset,0.10d0,irrep_x,
     &                             xnorm1)
           call tce_string_x2_norm(xc2(jvec),k_x2_offset,0.10d0,irrep_x,
     &                             xnorm2)
c normalization only on nodezero
           if(nodezero)  then
            xnorm=(1.0d0)/(dsqrt(xnorm1+xnorm2))
           endif
c
          call tce_string_x1(xc1(jvec),k_x1_offset,0.10d0,irrep_x,xnorm)
          call tce_string_x2(xc2(jvec),k_x2_offset,0.10d0,irrep_x,xnorm)
          if(nodezero) write(6,*)'-------------------------------------'
          call util_flush(6)
          endif !quantumc
!QuantumComputing
!
#ifdef IPCCSD
              else if(ipccsd) then
!              write(LuOut,*) "I am here 3.2"
!              call util_flush(LuOut)
               call tce_print_ipx1(xc1(jvec),k_x1_offset,0.10d0,irrep_x)
!              write(LuOut,*) "I am here 3.3"
!              call util_flush(LuOut)
               call tce_print_ipx2(xc2(jvec),k_x2_offset,0.10d0,irrep_x)
!              write(LuOut,*) "I am here 3.4"
!              call util_flush(LuOut)
#endif
#ifdef EACCSD
              else if(eaccsd) then
               call tce_print_eax1(xc1(jvec),k_x1_offset,0.10d0,irrep_x)
               call tce_print_eax2(xc2(jvec),k_x2_offset,0.10d0,irrep_x)
#endif
              endif
              call util_flush(LuOut)
             enddo
! ---
                goto 200
              endif
            enddo
!
           if(ccsd_var.eq.'lr-ccsd')then
            call deletefile(d_j2)
            if (.not.ma_pop_stack(l_j2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
            call deletefile(d_j1)
            if (.not.ma_pop_stack(l_j1_offset))
     1        call errquit("tce_energy: MA problem",37,MA_ERR)
           end if
!
         if(ccsd_var.ne.'cc2') then ! ----
          call deletefile(rec_d_9)
          call deletefile(rec_d_8)
          call deletefile(rec_d_7)
          call deletefile(rec_d_6)
          call deletefile(rec_d_5)
          call deletefile(rec_d_4)
          call deletefile(rec_d_3)
          call deletefile(rec_d_2)
          call deletefile(rec_d_1)
          call deletefile(rec_s_5)
          call deletefile(rec_s_4)
          call deletefile(rec_s_3)
          call deletefile(rec_s_2)
          call deletefile(rec_s_1)
         end if !---------
!
         call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
!
!
!
      else !  ------- EOM-SOLVER --------------- (solver 2)
!
! 
      if(nodezero) then
       write(LuOut,9510) eom_solver
       call util_flush(LuOut)
      end if
! === jaguar ===
      if(read_xa)  then
        nroots_reduced=nroots
        if(nodezero) then
         write(LuOut,9511) nroots_reduced
         call util_flush(LuOut)
        end if
        if (.not.rtdb_get(rtdb,'tce:omg',mt_dbl,nroots,eom_excit)) then
          call errquit('tce_energy: reading omg',nroots,CALC_ERR) 
        end if 
        if (.not.rtdb_get(rtdb,'tce:r0',mt_dbl,nroots,eom_r0)) then
          if(read_in3) then
           call errquit('tce: read_in3 active - r0 reading failed',
     &      nroots,CALC_ERR)
          else
           if(nodezero) then
           write(LuOut,9512)
           call util_flush(LuOut)
           end if 
          end if
        end if
        do ivec=1,nroots
        dbl_mb(k_omegax+ivec-1)=eom_excit(ivec)
        if(nodezero) then
         write(LuOut,*)'read excit. energy: ',dbl_mb(k_omegax+ivec-1)
         call util_flush(LuOut)
        end if
        enddo
!
       if(nroots_reduced.ge.1) then
        unitn=81
        ivec=1
        call tce_filenameindexed(ivec,'xc1',filename)
        call createfile(filename,xc1(ivec),size_x1)
        xc1_exist(ivec) = .true.
        call util_file_name0('x1am1',.false.,.true.,filename,fldgts)
        call read_tensor(filename,xc1(ivec),size_x1,unitn)
        unitn=82
        call tce_filenameindexed(ivec,'xc2',filename)
        call createfile(filename,xc2(ivec),size_x2)
        xc2_exist(ivec) = .true.
        call util_file_name0('x2am1',.false.,.true.,filename,fldgts)
        call read_tensor(filename,xc2(ivec),size_x2,unitn)
                if(nodezero)
     &                  write(LuOut,*)'EOMCCSD root nr. ',ivec
                call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
!
        if(nroots_reduced.ge.2) then
         unitn=83
         ivec=2
         call tce_filenameindexed(ivec,'xc1',filename)
         call createfile(filename,xc1(ivec),size_x1)
         xc1_exist(ivec) = .true.
         call util_file_name0('x1am2',.false.,.true.,filename,fldgts)
         call read_tensor(filename,xc1(ivec),size_x1,unitn)
         unitn=84
         call tce_filenameindexed(ivec,'xc2',filename)
         call createfile(filename,xc2(ivec),size_x2)
         xc2_exist(ivec) = .true.
         call util_file_name0('x2am2',.false.,.true.,filename,fldgts)
         call read_tensor(filename,xc2(ivec),size_x2,unitn)
                if(nodezero)
     &                  write(LuOut,*)'EOMCCSD root nr. ',ivec
                call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
!
         if(nroots_reduced.ge.3) then
          unitn=85
          ivec=3
          call tce_filenameindexed(ivec,'xc1',filename)
          call createfile(filename,xc1(ivec),size_x1)
          xc1_exist(ivec) = .true.
          call util_file_name0('x1am3',.false.,.true.,filename,fldgts)
          call read_tensor(filename,xc1(ivec),size_x1,unitn)
          unitn=86
          call tce_filenameindexed(ivec,'xc2',filename)
          call createfile(filename,xc2(ivec),size_x2)
          xc2_exist(ivec) = .true.
          call util_file_name0('x2am3',.false.,.true.,filename,fldgts)
          call read_tensor(filename,xc2(ivec),size_x2,unitn)
                if(nodezero)
     &                  write(LuOut,*)'EOMCCSD root nr. ',ivec
                call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
!
          if(nroots_reduced.ge.4) then
           unitn=87
           ivec=4
           call tce_filenameindexed(ivec,'xc1',filename)
           call createfile(filename,xc1(ivec),size_x1)
           xc1_exist(ivec) = .true.
           call util_file_name0('x1am4',.false.,.true.,filename,fldgts)
           call read_tensor(filename,xc1(ivec),size_x1,unitn)
           unitn=88
           call tce_filenameindexed(ivec,'xc2',filename)
           call createfile(filename,xc2(ivec),size_x2)
           xc2_exist(ivec) = .true.
           call util_file_name0('x2am4',.false.,.true.,filename,fldgts)
           call read_tensor(filename,xc2(ivec),size_x2,unitn)
                if(nodezero)
     &                  write(LuOut,*)'EOMCCSD root nr. ',ivec
                call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
          end if
         end if
        end if
       end if
!
        go to 200
      end if
! ==============
!
       modelname = "CIS-like initial guess for EOM-CCSD"
        if (nodezero) write(LuOut,9220)
     1  modelname(1:inp_strlen(modelname))
! == jaguar ===
      if(.not.eom_start) then
! =============
        call cis(LuOut,size_x1,k_x1_offset,d_f1,k_f1_offset,d_v2,
     1           k_v2_offset,k_omegax,l_omegax,d_rx1)
      else 
        nroots_reduced=nroots
        nxtrials=nroots
      end if
!
!
! == jaguar ===
        if(.not.eom_start) then
! =============
            if(nodezero) then
             write(LuOut,*)'INITIAL STARTS -----'
             call util_flush(LuOut)
            end if
            do iroot=1,nroots_reduced
             call tce_print_x1(x1_ini(iroot),k_x1_offset,0.10d0,irrep_x)
            enddo
            if(nodezero) then
             write(LuOut,*)'                    '
             write(LuOut,*)'INITIAL STARTS -----'
             call util_flush(LuOut)
            end if
        end if
!
!
            if (nxtrials.eq.0) goto 200
            call tce_eom_init
!
      if(ccsd_var.ne.'cc2') then
      ister=0
! rec_s_1
       CALL OFFSET_eomccsd_x1_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs1',-1,MA_ERR)
! rec_s_2
       CALL OFFSET_eomccsd_x1_2_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_2_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs2',-1,MA_ERR)
! rec_s_3
       CALL OFFSET_eomccsd_x1_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs3',-1,MA_ERR)
! rec_s_4
       CALL OFFSET_eomccsd_x1_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs4',-1,MA_ERR)
! rec_s_5
       CALL OFFSET_eomccsd_x1_7_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_7_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_s_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs5',-1,MA_ERR)
! rec_d_1
       CALL OFFSET_eomccsd_x2_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd1',-1,MA_ERR)
! rec_d_2
       CALL OFFSET_eomccsd_x2_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_3_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd2',-1,MA_ERR)
! rec_d_3
       CALL OFFSET_eomccsd_x2_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd3',-1,MA_ERR)
! rec_d_4
       CALL OFFSET_eomccsd_x2_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd4',-1,MA_ERR)
! rec_d_5
       CALL OFFSET_eomccsd_x2_6_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_6_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd5',-1,MA_ERR)
! rec_d_6
       CALL OFFSET_eomccsd_x2_8_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_6,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd6',-1,MA_ERR)
! rec_d_7
       CALL OFFSET_eomccsd_x2_8_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_3_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_7,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd7',-1,MA_ERR)
! rec_d_8
       CALL OFFSET_eomccsd_x2_8_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_4_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_8,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd8',-1,MA_ERR)
! rec_d_9
      CALL OFFSET_eomccsd_x2_10_1_1(l_aux_off,k_aux_off,size_aux)
      CALL TCE_FILENAME('eomccsd_x2_10_1_1_i2',filename)
      CALL CREATEFILE(filename,rec_d_9,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd9',-1,MA_ERR)
      end if 
!
      if(ccsd_var.eq.'lr-ccsd') then
            call tce_x1_offset(l_j1_offset,k_j1_offset,size_j1)
            call tce_filename('j1',filename)
            call createfile(filename,d_j1,size_j1)
!c            call gatoeaf(d_j1)
            call tce_zero(d_j1,size_j1)
!
            call tce_x2_offset(l_j2_offset,k_j2_offset,size_j2)
            call tce_filename('j2',filename)
            call createfile(filename,d_j2,size_j2)
!c            call gatoeaf(d_j2)
            call tce_zero(d_j2,size_j2)
!
            call tce_loc_j1_one(d_j1,k_j1_offset)
            call tce_loc_j2_one(d_j2,k_j2_offset)
!
            call tce_loc_j1(d_j1,k_j1_offset,d_t1,k_t1_offset)
            call tce_loc_j2(d_j2,k_j2_offset,d_t2,k_t2_offset)
!
            call tce_loc_j1_inv(d_j1,k_j1_offset)
            call tce_loc_j2_inv(d_j2,k_j2_offset)
      end if
!
!
            xdiis=diis
            do imicro=1,xdiis*nroots_reduced
              call tce_filenameindexed(imicro,'xp1',filename)
              call createfile(filename,xp1(imicro),size_x1)
              call gatodra(xp1(imicro))
              xp1_exist(imicro)=.true.
              call tce_filenameindexed(imicro,'xp2',filename)
              call createfile(filename,xp2(imicro),size_x2)
              call gatodra(xp2(imicro))
              xp2_exist(imicro)=.true.
              call tce_filenameindexed(imicro,'x1',filename)
              call createfile(filename,x1(imicro),size_x1)
              call tce_zero(x1(imicro),size_x1)
              call gatodra(x1(imicro))
              x1_exist(imicro)=.true.
              call tce_filenameindexed(imicro,'x2',filename)
              call createfile(filename,x2(imicro),size_x2)
              call tce_zero(x2(imicro),size_x2)
              call gatodra(x2(imicro))
              x2_exist(imicro)=.true.
            enddo
!
          do ivec=1,nroots_reduced
            call tce_filenameindexed(ivec,'xc1',filename)
            call createfile(filename,xc1(ivec),size_x1)
            call gatoeaf(xc1(ivec))
            call gatodra(xc1(ivec))
            xc1_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc2',filename)
            call createfile(filename,xc2(ivec),size_x2)
            call gatoeaf(xc2(ivec))
            call gatodra(xc2(ivec))
            xc2_exist(ivec) = .true.
          enddo
! === jaguar ===
        if(eom_start) then
        if(nodezero) then
         write(LuOut,*)'IN.EOMSD READ IN nroot_reduced=',nroots_reduced
         call util_flush(LuOut)
        end if
        nroots_reduced=nroots
        if(nroots_reduced.ge.1) then
        unitn=81
        ivec=1
        call util_file_name0('x1am1',.false.,.true.,filename,fldgts)
        call read_tensor(filename,xc1(ivec),size_x1,unitn)
        unitn=82
        call util_file_name0('x2am1',.false.,.true.,filename,fldgts)
        call read_tensor(filename,xc2(ivec),size_x2,unitn)
                if(nodezero)
     &                  write(LuOut,*)'INI EOMCCSD root nr. ',ivec
                call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
        if(nroots_reduced.ge.2) then
         unitn=83
         ivec=2
         call util_file_name0('x1am2',.false.,.true.,filename,fldgts)
         call read_tensor(filename,xc1(ivec),size_x1,unitn)
         unitn=84
         call util_file_name0('x2am2',.false.,.true.,filename,fldgts)
         call read_tensor(filename,xc2(ivec),size_x2,unitn)
                if(nodezero)
     &                  write(LuOut,*)'INI EOMCCSD root nr. ',ivec
                call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
         if(nroots_reduced.ge.3) then
          unitn=85
          ivec=3
          call util_file_name0('x1am3',.false.,.true.,filename,fldgts)
          call read_tensor(filename,xc1(ivec),size_x1,unitn)
          unitn=86
          call util_file_name0('x2am3',.false.,.true.,filename,fldgts)
          call read_tensor(filename,xc2(ivec),size_x2,unitn)
                if(nodezero)
     &                  write(LuOut,*)'INI EOMCCSD root nr. ',ivec
                call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
          if(nroots_reduced.ge.4) then
           unitn=87
           ivec=4
           call util_file_name0('x1am4',.false.,.true.,filename,fldgts)
           call read_tensor(filename,xc1(ivec),size_x1,unitn)
           unitn=88
           call util_file_name0('x2am4',.false.,.true.,filename,fldgts)
           call read_tensor(filename,xc2(ivec),size_x2,unitn)
                if(nodezero)
     &                  write(LuOut,*)'INI EOMCCSD root nr. ',ivec
                call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
          end if
         end if
        end if
       end if
       end if ! if eom_start
! ====================
!
          if(ccsd_var.eq.'lr-ccsd') then
            modelname = "LR-EOMCCSD right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
          else if(ccsd_var.eq.'cc2') then
            modelname = "CC2 right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
          else
            modelname = "EOM-CCSD right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
          end if
!
! open k_h_eff  and k_c_eff
!
        if (.not.MA_PUSH_GET(mt_dbl,
     1      xdiis*xdiis*nroots_reduced*nroots_reduced,'h_eff',
     1      l_h_eff,k_h_eff)) call errquit('h_eff',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,xdiis*nroots_reduced*nroots_reduced,
     1      'c_eff',
     1      l_c_eff,k_c_eff)) call errquit('c_eff',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,nroots_reduced,'residual',
     1          l_residual,k_residual)) call errquit('l_res',1,MA_ERR)
!
            iter=1
!
!cx             if(iter.eq.1) then
!cx              do iroot=1,nroots_reduced
!cx               call copyfile(x1_ini(iroot),xc1(iroot),size_x1)
!cx               call tce_zero(xc2(iroot),size_x2)
!cx              enddo
!cx             endif 
!
       if(.not.eom_start) then
             if(iter.eq.1) then
              do iroot=1,nroots_reduced
               call copyfile(x1_ini(iroot),xc1(iroot),size_x1)
               call tce_zero(xc2(iroot),size_x2)
              enddo
             endif
       end if
!
!
  111        continue  !iter continuation
!
             imicro = 1
!
!
            do iroot=1,nroots_reduced
             call copyfile(xc1(iroot),x1(iroot),size_x1)
             call copyfile(xc2(iroot),x2(iroot),size_x2)
            enddo
!
            call norm_x_vector_mr(k_x1_offset,k_x2_offset,
     1                            size_x1,size_x2,nroots_reduced)
!
            do iroot=1,nroots_reduced
             call tce_zero(xp1(iroot),size_x1)
             call tce_zero(xp2(iroot),size_x2)
            enddo
!
             imicro = 1
!
             do i=1,xdiis*nroots_reduced*nroots_reduced
              dbl_mb(k_c_eff+i-1)=0.0d0
             enddo
             do i=1,nroots_reduced
             dbl_mb(k_c_eff+(i-1)*xdiis*nroots_reduced+i-1)=1.0d0
             enddo
             do i=1,xdiis*xdiis*nroots_reduced*nroots_reduced
               dbl_mb(k_h_eff+i-1)=0.0d0
             enddo
!
!    declare h_eff c_eff - MA 
!
      cpu=-util_cpusec()
      wall=-util_wallsec()
!
      do iroot=1,nroots_reduced ! x-x-x-x-x-x-x-x-x
!
! - T1/X1 LOCALIZATION ----------
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,size_x1,'x1_local',
     1      l_x1_local,k_x1_local))
     1      call errquit('x1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
        call ma_zero(dbl_mb(k_x1_local),size_x1)
!c        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
!c        call ga_get(x1(iroot),1,size_x1,1,1,dbl_mb(k_x1_local),1)
        call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
        call get_block(x1(iroot),dbl_mb(k_x1_local),size_x1,0)
! -------------------------------
         if(ccsd_var.eq.'cc2') then
          call cc2_x1(d_f1,xp1(iroot),k_t1_local,d_t2,d_v2,
     &      k_x1_local,x2(iroot),
     &      k_f1_offset,k_x1_offset,k_t1_offset,k_t2_offset,
     &      k_v2_offset,k_x1_offset,k_x2_offset)
         else
          if(ccsd_var.eq.'ic') then
          call eomicsd_x1(d_f1,xp1(iroot),k_t1_local,d_t2,d_v2,
     1      k_x1_local,x2(iroot),k_f1_offset,k_x1_offset,
     2      k_t1_offset,k_t2_offset,k_v2_offset,
     3      k_x1_offset,k_x2_offset,
     4      ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5      rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6      rec_d_6,rec_d_7,rec_d_8,rec_d_9)
          else
          call eomccsd_x1_n(d_f1,xp1(iroot),k_t1_local,d_t2,d_v2,
     1      k_x1_local,x2(iroot),k_f1_offset,k_x1_offset,
     2      k_t1_offset,k_t2_offset,k_v2_offset,
     3      k_x1_offset,k_x2_offset,
     4      ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5      rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6      rec_d_6,rec_d_7,rec_d_8,rec_d_9)
          end if
         end if
         if(ccsd_var.eq.'lr-ccsd') then
            call tce_j1_x1_scaling(xp1(iroot),k_x1_offset,
     1           d_j1,k_j1_offset)
         end if
         if(ccsd_var.eq.'cc2') then
           call cc2_x2(d_f1,xp2(iroot),k_t1_local,d_v2,
     &       k_x1_local,x2(iroot),
     &       k_f1_offset,k_x2_offset,k_t1_offset,k_v2_offset,
     &       k_x1_offset,k_x2_offset,size_x1,size_x2)
         else
          if(ccsd_var.eq.'ic') then
           irrep_c=irrep_x
           call tce_filename('c2',filename)
           call createfile(filename,d_c2,size_x2)
           call tce_zero(d_c2,size_x2)
           call c2excit2_loc(d_c2,k_t1_local,k_x1_local,
     &       x2(iroot),k_x2_offset,k_t1_offset,k_x1_offset,
     &       k_x2_offset)
           call eomicsd_x2(d_f1,xp2(iroot),k_t1_local,d_t2,d_v2,
     1       k_x1_local,x2(iroot),k_f1_offset,k_x2_offset,
     2       k_t1_offset,k_t2_offset,k_v2_offset,
     3       k_x1_offset,k_x2_offset,size_x1,size_x2,
     4       ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9,d_c2)
           call deletefile(d_c2)
          else
           call eomccsd_x2_n(d_f1,xp2(iroot),k_t1_local,d_t2,d_v2,
     1       k_x1_local,x2(iroot),k_f1_offset,k_x2_offset,
     2       k_t1_offset,k_t2_offset,k_v2_offset,
     3       k_x1_offset,k_x2_offset,size_x1,size_x2,
     4       ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9)
          end if
         end if
         if(ccsd_var.eq.'lr-ccsd') then
           call tce_j1_x1_scaling(xp1(iroot),k_x1_offset,
     1          d_j1,k_j1_offset)
         end if
!
         ister=1
! - T1/X1 LOCALIZATION --
         if(.not.MA_POP_STACK(l_x1_local))
     &      call errquit('l_x1_local',4,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
! -----------------------------
         enddo ! x-x-x-x-x-x-x-x-x
!   calculate  H-bar delta x
!
           call form_heff_mr(k_x1_offset,k_x2_offset,size_x1,size_x2,
     &        dbl_mb(k_c_eff),dbl_mb(k_h_eff),dbl_mb(k_omegax),
     &        imicro,xdiis,iter)
!
!
          do imicro=2,xdiis !microloop(imicro)
!
          do iroot=1,nroots_reduced
           ipos=(imicro-1)*nroots_reduced+iroot
           call tce_zero(x1(ipos),size_x1)
           call tce_zero(x2(ipos),size_x2)
           call tce_zero(xp1(ipos),size_x1)
           call tce_zero(xp2(ipos),size_x2)
          enddo
!
          do iroot=1,nroots_reduced
           dbl_mb(k_residual+iroot-1)=0.0d0
          enddo
!
          do iroot=1,nroots_reduced
           call form_correction_mr(k_x1_offset,k_x2_offset,
     1                  size_x1,size_x2,
     1                  dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     1                  dbl_mb(k_omegax+iroot-1),
     1                  imicro,xdiis,dbl_mb(k_residual+iroot-1),iroot)
          enddo
!
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
      if(nodezero) then
                  write(LuOut,9211) iter
             do iroot=1,nroots_reduced
                if(iroot.eq.nroots_reduced) then
                  write(LuOut,9230) dbl_mb(k_residual+iroot-1),
     2            dbl_mb(k_omegax+iroot-1),
     3            dbl_mb(k_omegax+iroot-1)*au2ev,cpu,wall
                else
                  write(LuOut,9230) dbl_mb(k_residual+iroot-1),
     2            dbl_mb(k_omegax+iroot-1),
     3            dbl_mb(k_omegax+iroot-1)*au2ev
                end if
             enddo
      end if
      call util_flush(LuOut)
      cpu=-util_cpusec()
      wall=-util_wallsec()
!
           converged=.true.
!cx           do iroot=1,nroots_reduced
!cx             if(dabs(dbl_mb(k_residual+iroot-1)).gt.thresh) 
!cx     &             converged=.false. 
!cx           enddo 

         if(nroots_reduced.gt.1) then 
           do iroot=1,nroots_reduced-1
             if(dabs(dbl_mb(k_residual+iroot-1)).gt.thresh_eom)
     &             converged=.false.
           enddo
           if(dabs(dbl_mb(k_residual+nroots_reduced-1)).gt.thresh_l)
     7             converged=.false.
         else
           if(dabs(dbl_mb(k_residual+nroots_reduced-1)).gt.thresh_eom)
     7             converged=.false.
         end if
!
            if(converged) then
             do iroot=1,nroots_reduced
              call form_best_guess_mr(k_x1_offset,k_x2_offset,
     &                size_x1,size_x2,
     &                dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     &                imicro-1,
     &                xdiis,iroot)
             enddo
!
               do iroot=1,nroots_reduced
ccx                if(nodezero)
ccx     &                  write(LuOut,*)'EOMCCSD root nr. ',iroot
c right EOMCCSD  printout
             if (nodezero) then
               write(LuOut,*)'------------------'
               write(LuOut,9460)
     1           iroot,dbl_mb(k_omegax+iroot-1),
     2           dbl_mb(k_omegax+iroot-1)*au2ev
             endif
c
               if(nodezero) then
                write(LuOut,*)'largest EOMCCSD amplitudes: R1 and R2'
                call util_flush(LuOut)
               end if
c
                call tce_print_x1(xc1(iroot),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(iroot),k_x2_offset,0.10d0,irrep_x)
               enddo
              go to 112
            end if
! all together
            call deltax_ort_norm_mr(k_x1_offset,k_x2_offset,
     1                              size_x1,size_x2,
     1                              imicro,xdiis)
!
!
!
         do iroot=1,nroots_reduced ! --- iroot
!
          ipos=(imicro-1)*nroots_reduced+iroot
!
! - T1/X1 LOCALIZATION ----------
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,size_x1,'x1_local',
     1      l_x1_local,k_x1_local))
     1      call errquit('x1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
        call ma_zero(dbl_mb(k_x1_local),size_x1)
!c        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
!c        call ga_get(x1(ipos),1,size_x1,1,1,dbl_mb(k_x1_local),1)
        call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
        call get_block(x1(ipos),dbl_mb(k_x1_local),size_x1,0)
! -------------------------------
         if(ccsd_var.eq.'cc2') then
          call cc2_x1(d_f1,xp1(ipos),k_t1_local,d_t2,d_v2,
     &      k_x1_local,x2(ipos),
     &      k_f1_offset,k_x1_offset,k_t1_offset,k_t2_offset,
     &      k_v2_offset,k_x1_offset,k_x2_offset)
         else
          if(ccsd_var.eq.'ic') then
          call eomicsd_x1(d_f1,xp1(ipos),k_t1_local,d_t2,d_v2,
     1      k_x1_local,x2(ipos),k_f1_offset,k_x1_offset,
     2      k_t1_offset,k_t2_offset,k_v2_offset,
     3      k_x1_offset,k_x2_offset,
     4      ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5      rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6      rec_d_6,rec_d_7,rec_d_8,rec_d_9)
          else
          call eomccsd_x1_n(d_f1,xp1(ipos),k_t1_local,d_t2,d_v2,
     1      k_x1_local,x2(ipos),k_f1_offset,k_x1_offset,
     2      k_t1_offset,k_t2_offset,k_v2_offset,
     3      k_x1_offset,k_x2_offset,
     4      ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5      rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6      rec_d_6,rec_d_7,rec_d_8,rec_d_9)
           end if
         end if
         if(ccsd_var.eq.'lr-ccsd') then
            call tce_j1_x1_scaling(xp1(ipos),k_x1_offset,
     1           d_j1,k_j1_offset)
         end if
         if(ccsd_var.eq.'cc2') then
           call cc2_x2(d_f1,xp2(ipos),k_t1_local,d_v2,
     &       k_x1_local,x2(ipos),
     &       k_f1_offset,k_x2_offset,k_t1_offset,k_v2_offset,
     &       k_x1_offset,k_x2_offset,size_x1,size_x2)
         else
          if(ccsd_var.eq.'ic') then
           irrep_c=irrep_x
           call tce_filename('c2',filename)
           call createfile(filename,d_c2,size_x2)
           call tce_zero(d_c2,size_x2)
           call c2excit2_loc(d_c2,k_t1_local,k_x1_local,
     &       x2(ipos),k_x2_offset,k_t1_offset,k_x1_offset,
     &       k_x2_offset)
           call eomicsd_x2(d_f1,xp2(ipos),k_t1_local,d_t2,d_v2,
     1       k_x1_local,x2(ipos),k_f1_offset,k_x2_offset,
     2       k_t1_offset,k_t2_offset,k_v2_offset,
     3       k_x1_offset,k_x2_offset,size_x1,size_x2,
     4       ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9,d_c2)
           call deletefile(d_c2)
          else
           call eomccsd_x2_n(d_f1,xp2(ipos),k_t1_local,d_t2,d_v2,
     1       k_x1_local,x2(ipos),k_f1_offset,k_x2_offset,
     2       k_t1_offset,k_t2_offset,k_v2_offset,
     3       k_x1_offset,k_x2_offset,size_x1,size_x2,
     4       ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9)
          end if
         end if
         if(ccsd_var.eq.'lr-ccsd') then
           call tce_j1_x1_scaling(xp1(ipos),k_x1_offset,
     1          d_j1,k_j1_offset)
         end if
!
         ister=1
! - T1/X1 LOCALIZATION --
         if(.not.MA_POP_STACK(l_x1_local))
     &      call errquit('l_x1_local',4,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
! -----------------------
!
        enddo ! --- iroot 
!
!             calculate  H-bar delta x
!
           call form_heff_mr(k_x1_offset,k_x2_offset,size_x1,size_x2,
     &        dbl_mb(k_c_eff),dbl_mb(k_h_eff),dbl_mb(k_omegax),
     &        imicro,xdiis,iter)
!
           iter=iter+1                              
             if(iter.gt.maxiter) then
              if(nodezero) then
                 write(LuOut,9233) iroot
              end if
             do iroot=1,nroots_reduced
              call form_best_guess_mr(k_x1_offset,k_x2_offset,
     &                size_x1,size_x2,
     &                dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     &                imicro,
     &                xdiis,iroot)
               call tce_print_x1(xc1(iroot),k_x1_offset,0.10d0,irrep_x)
               call tce_print_x2(xc2(iroot),k_x2_offset,0.10d0,irrep_x)
             enddo
              go to 112
             end if
!
          enddo   !microloop(imicro)
!
             do iroot=1,nroots_reduced
              call form_best_guess_mr(k_x1_offset,k_x2_offset,
     &                size_x1,size_x2,
     &                dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     &                xdiis,
     &                xdiis,iroot)
             enddo
!
               go to 111
!
  112         continue !from converged
! === jaguar ===
      if(write_xa)  then
       if(nroots_reduced.ge.1) then
        unitn=81
        call util_file_name0('x1am1',.false.,.true.,filename,fldgts)
        call write_tensor(filename,xc1(1),size_x1,unitn)
        unitn=82
        call util_file_name0('x2am1',.false.,.true.,filename,fldgts)
        call write_tensor(filename,xc2(1),size_x2,unitn)
        if(nroots_reduced.ge.2) then
         unitn=83
         call util_file_name0('x1am2',.false.,.true.,filename,fldgts)
         call write_tensor(filename,xc1(2),size_x1,unitn)
         unitn=84
         call util_file_name0('x2am2',.false.,.true.,filename,fldgts)
         call write_tensor(filename,xc2(2),size_x2,unitn)
         if(nroots_reduced.ge.3) then
          unitn=85
          call util_file_name0('x1am3',.false.,.true.,filename,fldgts)
          call write_tensor(filename,xc1(3),size_x1,unitn)
          unitn=86
          call util_file_name0('x2am3',.false.,.true.,filename,fldgts)
          call write_tensor(filename,xc2(3),size_x2,unitn)
          if(nroots_reduced.ge.4) then
           unitn=87
           call util_file_name0('x1am4',.false.,.true.,filename,fldgts)
           call write_tensor(filename,xc1(4),size_x1,unitn)
           unitn=88
           call util_file_name0('x2am4',.false.,.true.,filename,fldgts)
           call write_tensor(filename,xc2(4),size_x2,unitn)
          end if
         end if
        end if
       end if
      end if
!
         if(.not.MA_POP_STACK(l_residual))
     &      call errquit('l_residual',4,MA_ERR)
         if(.not.MA_POP_STACK(l_c_eff))
     &      call errquit('l_c_eff',4,MA_ERR)
         if(.not.MA_POP_STACK(l_h_eff))
     &      call errquit('l_h_eff',4,MA_ERR)
!
            do imicro=1,xdiis
              call deletefile(xp1(imicro))
              xp1_exist(imicro)=.false.
              call deletefile(xp2(imicro))
              xp2_exist(imicro)=.false.
              call deletefile(x1(imicro))
              x1_exist(imicro)=.false.
              call deletefile(x2(imicro))
              x2_exist(imicro)=.false.
            enddo
!
!
!
! xc1 xc2 are deleted by tce_eom_tidy (later on)
!
!
           if(ccsd_var.eq.'lr-ccsd')then
            call deletefile(d_j2)
            if (.not.ma_pop_stack(l_j2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
            call deletefile(d_j1)
            if (.not.ma_pop_stack(l_j1_offset))
     1        call errquit("tce_energy: MA problem",37,MA_ERR)
           end if
!
!
         if(ccsd_var.ne.'cc2') then ! ----
          call deletefile(rec_d_9)
          call deletefile(rec_d_8)
          call deletefile(rec_d_7)
          call deletefile(rec_d_6)
          call deletefile(rec_d_5)
          call deletefile(rec_d_4)
          call deletefile(rec_d_3)
          call deletefile(rec_d_2)
          call deletefile(rec_d_1)
          call deletefile(rec_s_5)
          call deletefile(rec_s_4)
          call deletefile(rec_s_3)
          call deletefile(rec_s_2)
          call deletefile(rec_s_1)
         end if !---------
!
         do ivec=1,maxtrials
           if(x1_ini_exist(ivec))  then
            call deletefile(x1_ini(ivec))
            x1_ini_exist(ivec) = .false.
           end if
         enddo
!
!
!
         end if ! ------- EOM-SOLVER ---------------
!
!
!
! -------------
! EOM-CCSD    left
! -------------
  200       continue
!
         if(eom_solver.eq.1) then
           if(ccsd_var.eq.'lr-ccsd')then
            call deletefile(d_j2)
            if (.not.ma_pop_stack(l_j2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
            call deletefile(d_j1)
            if (.not.ma_pop_stack(l_j1_offset))
     1        call errquit("tce_energy: MA problem",37,MA_ERR)
           end if
           if(ccsd_var.ne.'cc2') then ! ----
            call deletefile(rec_d_9)
            call deletefile(rec_d_8)
            call deletefile(rec_d_7)
            call deletefile(rec_d_6)
            call deletefile(rec_d_5)
            call deletefile(rec_d_4)
            call deletefile(rec_d_3)
            call deletefile(rec_d_2)
            call deletefile(rec_d_1)
            call deletefile(rec_s_5)
            call deletefile(rec_s_4)
            call deletefile(rec_s_3)
            call deletefile(rec_s_2)
            call deletefile(rec_s_1)
           end if !---------
         end if
!
            if (left) then
            call tce_eom_yguess(needt1,needt2,needt3,needt4,
     1        size_y1,size_y2,size_y3,size_y4,
     2        k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset)
            if (nytrials.eq.0) goto 400
            modelname = "EOM-CCSD left-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
!           call hbar_init
            do iter=1,maxiter
              if (nodezero.and.util_print('eom',print_default))
     1          write(LuOut,9210) iter,nytrials
              do ivec = 1,nytrials
                if (.not.yp1_exist(ivec)) then
                  call tce_filenameindexed(ivec,'yp1',filename)
                  call createfile(filename,yp1(ivec),size_y1)
                  yp1_exist(ivec) = .true.
                  call dratoga(y1(ivec))
                  call dratoga(y2(ivec))
                  call eomccsd_y1(d_f1,yp1(ivec),d_t1,d_t2,d_v2,
     1              y1(ivec),y2(ivec),k_f1_offset,k_y1_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,
     3              k_y1_offset,k_y2_offset)
                  call reconcilefile(yp1(ivec),size_y1)
                  call gatodra(y2(ivec))
                  call gatodra(y1(ivec))
                  call gatodra(yp1(ivec))
                endif
                if (.not.yp2_exist(ivec)) then
                  call tce_filenameindexed(ivec,'yp2',filename)
                  call createfile(filename,yp2(ivec),size_y2)
                  yp2_exist(ivec) = .true.
                  call dratoga(y1(ivec))
                  call dratoga(y2(ivec))
                  call eomccsd_y2(d_f1,yp2(ivec),d_t1,d_t2,d_v2,
     1              y1(ivec),y2(ivec),k_f1_offset,k_y2_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,
     3              k_y1_offset,k_y2_offset)
                  call reconcilefile(yp2(ivec),size_y2)
                  call gatodra(y2(ivec))
                  call gatodra(y1(ivec))
                  call gatodra(yp2(ivec))
                endif
              enddo
              if (.not.ma_push_get(mt_dbl,nytrials,'residual',
     1          l_residual,k_residual))
     2          call errquit('tce_energy: MA problem',101,MA_ERR)
              call tce_eom_ydiagon(needt1,needt2,needt3,needt4,
     1          size_y1,size_y2,size_y3,size_y4,
     2          k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset,
     3          d_ry1,d_ry2,d_ry3,d_ry4,
     4          dbl_mb(k_omegay),dbl_mb(k_residual))
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
              converged = .true.
              do ivec = 1,nroots_reduced
                if (nodezero.and.(ivec.ne.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegay+ivec-1),
     3            dbl_mb(k_omegay+ivec-1)*au2ev
                if (nodezero.and.(ivec.eq.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegay+ivec-1),
     3            dbl_mb(k_omegay+ivec-1)*au2ev,cpu,wall
                if (nodezero) call util_flush(LuOut)
                if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1            converged = .false.
              enddo
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (.not.ma_pop_stack(l_residual))
     1          call errquit("tce_energy: MA problem",102,MA_ERR)
              if (converged) then
                if (nodezero) write(LuOut,9240)
                call tce_eom_ytidy
!               call hbar_tidy
                do ivec = 1,nroots_reduced
                  if (nodezero) write(LuOut,9460)
     1              ivec,dbl_mb(k_omegax+ivec-1),
     2              dbl_mb(k_omegax+ivec-1)*au2ev
                  do jvec = 1,nroots_reduced
                  if (dabs(dbl_mb(k_omegax+ivec-1)-
     1                dbl_mb(k_omegay+jvec-1)).lt.thresh) then
                    call tce_filename('d0',filename)
                    call createfile(filename,d_d0,size_d0)
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
                    call dratoga(yc1(jvec))
                    call dratoga(yc2(jvec))
                    call eomccsd_denominator(d_d0,xc1(ivec),xc2(ivec),
     1                yc1(jvec),yc2(jvec),k_d0_offset,k_x1_offset,
     2                k_x2_offset,k_y1_offset,k_y2_offset)
                    call reconcilefile(d_d0,size_d0)
                    call gatodra(yc2(jvec))
                    call gatodra(yc1(jvec))
                    call gatodra(xc2(ivec))
                    call gatodra(xc1(ivec))
                    call get_block(d_d0,denominator,1,0)
                    call deletefile(d_d0)
                    if (dabs(denominator).gt.thresh) then
                    do axis=1,3
                      irrep_d = sym_abelian_axis(geom,axis)
                      irrep_x = irrep
                      irrep_y = irrep
                      call tce_filename('d0',filename)
                      call createfile(filename,d_d0,size_d0)
                      call tce_filename('x0',filename)
                      call createfile(filename,d_x0,size_x0)
                      call tce_filename('y0',filename)
                      call createfile(filename,d_y0,size_y0)
                      call dratoga(xc1(ivec))
                      call dratoga(xc2(ivec))
                      call dratoga(yc1(jvec))
                      call dratoga(yc2(jvec))
                      call eomccsd_density1(d_d1(axis),d_d0,d_t1,d_t2,
     1                  d_x0,xc1(ivec),xc2(ivec),d_y0,yc1(jvec),
     2                  yc2(jvec),k_d1_offset(axis),k_d0_offset,
     3                  k_t1_offset,k_t2_offset,k_x0_offset,k_x1_offset,
     4                  k_x2_offset,k_y0_offset,k_y1_offset,k_y2_offset)
                      call reconcilefile(d_d0,size_d0)
                      call gatodra(yc2(jvec))
                      call gatodra(yc1(jvec))
                      call gatodra(xc2(ivec))
                      call gatodra(xc1(ivec))
                      call get_block(d_d0,numerator,1,0)
                      call deletefile(d_d0)
                      dipole_exc(axis) = dipole_scf(axis)
     1                                 + numerator/denominator
                      transition(axis) = 1.0d0
                      irrep_x = irrep   !! Caution !!
                      irrep_y = 0       !! Caution !!
                      call tce_filename('d0',filename)
                      call createfile(filename,d_d0,size_d0)
                      call put_block(d_x0,0.0d0,1,0)
                      call put_block(d_y0,1.0d0,1,0)
                      call dratoga(xc1(ivec))
                      call dratoga(xc2(ivec))
                      call eomccsd_density1(d_d1(axis),d_d0,d_t1,d_t2,
     1                  d_x0,xc1(ivec),xc2(ivec),d_y0,d_lambda1,
     2                  d_lambda2,k_d1_offset(axis),k_d0_offset,
     3                  k_t1_offset,k_t2_offset,k_x0_offset,k_x1_offset,
     4                  k_x2_offset,k_y0_offset,k_l1_offset,k_l2_offset)
                      call reconcilefile(d_d0,size_d0)
                      call gatodra(xc2(ivec))
                      call gatodra(xc1(ivec))
                      call get_block(d_d0,numerator,1,0)
                      call deletefile(d_d0)
                      transition(axis) = transition(axis) * numerator
                      irrep_x = 0       !! Caution !!
                      irrep_y = irrep   !! Caution !!
                      call tce_filename('d0',filename)
                      call createfile(filename,d_d0,size_d0)
                      call put_block(d_x0,1.0d0,1,0)
                      call put_block(d_y0,0.0d0,1,0)
                      call tce_t1_offset(l_r1_offset,k_r1_offset,
     1                  size_r1)
                      call tce_filename('r1',filename)
                      call createfile(filename,d_r1,size_r1)
                      call tce_t2_offset(l_r2_offset,k_r2_offset,
     1                  size_r2)
                      call tce_filename('r2',filename)
                      call createfile(filename,d_r2,size_r2)
                      call dratoga(yc1(jvec))
                      call dratoga(yc2(jvec))
                      call eomccsd_density1(d_d1(axis),d_d0,d_t1,d_t2,
     1                  d_x0,d_r1,d_r2,d_y0,yc1(jvec),
     2                  yc2(jvec),k_d1_offset(axis),k_d0_offset,
     3                  k_t1_offset,k_t2_offset,k_x0_offset,k_r1_offset,
     4                  k_r2_offset,k_y0_offset,k_y1_offset,k_y2_offset)
                      call reconcilefile(d_d0,size_d0)
                      call gatodra(yc2(jvec))
                      call gatodra(yc1(jvec))
                      call get_block(d_d0,numerator,1,0)
                      call deletefile(d_r2)
                      if (.not.ma_pop_stack(l_r2_offset))
     1                  call errquit("tce_energy: MA problem",100,
     2                  MA_ERR)
                      call deletefile(d_r1)
                      if (.not.ma_pop_stack(l_r1_offset))
     1                  call errquit("tce_energy: MA problem",101,
     2                  MA_ERR)
                      call deletefile(d_d0)
                      transition(axis) = transition(axis)
     1                                 * numerator / denominator
                      call deletefile(d_y0)
                      call deletefile(d_x0)
                    enddo
                    if (nodezero) write(LuOut,9470) "EOM-CCSD",
     1                      dsqrt(dabs(transition(1))),
     2                      dsqrt(dabs(transition(2))),
     3                      dsqrt(dabs(transition(3))),
     4                      (transition(1)+transition(2)+transition(3))
     5                      *dbl_mb(k_omegax+ivec-1)*2.0d0/3.0d0
                    dipole_exc_t = dsqrt( dipole_exc(1)*dipole_exc(1) +
     1                                    dipole_exc(2)*dipole_exc(2) +
     2                                    dipole_exc(3)*dipole_exc(3) )
                    if (nodezero) write(LuOut,9430) "EOM-CCSD",
     1                      dipole_exc(1),dipole_exc(1)*debye,
     2                      dipole_exc(2),dipole_exc(2)*debye,
     3                      dipole_exc(3),dipole_exc(3)*debye,
     4                      dipole_exc_t ,dipole_exc_t*debye
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
!                    irrep_x = irrep   !! Caution !!
!                    call tce_print_x1(xc1(ivec),k_x1_offset,printtol)
!                    call tce_print_x2(xc2(ivec),k_x2_offset,printtol)
                    call tce_print_x1(xc1(ivec),k_x1_offset,printtol,
     1                                irrep)
                    call tce_print_x2(xc2(ivec),k_x2_offset,printtol,
     1                                irrep)
                    call gatodra(xc2(ivec))
                    call gatodra(xc1(ivec))
                    if (nodezero) call util_flush(LuOut)
                    endif
                  endif
                  enddo
!
!->d3p975
!
!   eom-ccsd 1-particle reduced density matrix/ces
!
                  if (idens) then
                     if (.not.intorb)
     1                  call eomccsd_1prdm(rtdb,dbl_mb(k_omegax+ivec-1),
     1                               d_f1,d_v2,
     2                               k_f1_offset,k_v2_offset,
     2                               d_t1,d_t2,xc1(ivec),xc2(ivec),
     1                               yc1(ivec),yc2(ivec),
     2                               k_t1_offset,k_t2_offset,
     3                               k_x1_offset,k_x2_offset,
     4                               k_y1_offset,k_y2_offset)
                  endif
! 
!<-d3p975
!
                enddo
                modelname = "EOM-CCSD"
                goto 400
              endif
            enddo
            call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
            endif


! -------------------------
! CR-EOMCCSD(T)
! -------------------------
        if ((perturbative.eq.'creom_(t)').or.
     &      (perturbative.eq.'emb1').or.(perturbative.eq.'emb2')) then
        irrep_c = irrep_x  !important
        go to 8922
!cccx      if(irrep_x.eq.irrep_g) then !ground-state calculations
!      if(irrep_x.eq.irrep_t) then !ground-state calculations
!            call tce_x2_offset(l_c2_offset,k_c2_offset,size_c2)
!            call tce_filename('c2',filename)
!            call createfile(filename,d_c2,size_c2)
!cc            call gatoeaf(d_c2)
!            call tce_zero(d_c2,size_c2)
!       call tce_e_offset(l_e_offset,k_e_offset,size_e)
!       call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
!       call tce_e_offset(l_ex2_offset,k_ex2_offset,size_ex2)
!          call tce_filename('e',filename)
!          call createfile(filename,d_e,size_e)
!          call tce_zero(d_e,size_e)
!          call tce_filename('ex1',filename)
!          call createfile(filename,d_ex1,size_ex1)
!          call tce_zero(d_ex1,size_ex1)
!          call tce_filename('ex2',filename)
!          call createfile(filename,d_ex2,size_ex2)
!          call tce_zero(d_ex2,size_ex2)
!      pt3_6=0.0d0
!      call cr_eomccsd_t_gr(d_t1,k_t1_offset,d_t2,k_t2_offset,
!     2                     d_f1,k_f1_offset,d_v2,k_v2_offset,
!     3                     d_e,k_e_offset,
!     4                     d_ex1,k_ex1_offset,size_ex1,
!     5                     d_ex2,k_ex2_offset,size_ex2,
!     7                     d_c2,k_c2_offset,size_c2,
!     8                     pt3_6)
!       energy_grcr=ref+corr+pt3_6
!       if (nodezero) then
!       write(LuOut,9080) "ground-state CR-CCSD(T)IA",ref+corr+pt3_6
!       write(LuOut,*)
!       call util_flush(LuOut)
!       end if
!c
!          call deletefile(d_ex2)
!          call deletefile(d_ex1)
!          call deletefile(d_e)
!          if (.not.ma_pop_stack(l_ex2_offset))
!     1     call errquit("tce_energy: MA problem",24,MA_ERR)
!          if (.not.ma_pop_stack(l_ex1_offset))
!     1     call errquit("tce_energy: MA problem",24,MA_ERR)
!          if (.not.ma_pop_stack(l_e_offset))
!     1     call errquit("tce_energy: MA problem",24,MA_ERR)
!c
!            call deletefile(d_c2)
!            if (.not.ma_pop_stack(l_c2_offset))
!     1        call errquit("tce_energy: MA problem",36,MA_ERR)
!c
!      end if
! UWAGA
! ACHTUNG
! BHUMAHE
! ATTENTION
 8922  continue
!
        DO 3443,ivec = 1,nroots_reduced
          cpu = - util_cpusec()
          wall = - util_wallsec()
! d_e,d_ex1,d_ex2 global arrays
!
            call tce_x1_offset(l_c1_offset,k_c1_offset,size_c1)
            call tce_filename('c1',filename)
            call createfile(filename,d_c1,size_c1)
!c            call gatoeaf(d_c1)
            call tce_zero(d_c1,size_c1)
!
            call tce_x2_offset(l_c2_offset,k_c2_offset,size_c2)
            call tce_filename('c2',filename)
            call createfile(filename,d_c2,size_c2)
!c            call gatoeaf(d_c2)
            call tce_zero(d_c2,size_c2)
!
      pt3_5 = 0.0d0
      pt3_6 = 0.0d0
      excit = dbl_mb(k_omegax+ivec-1)
      call tce_e_offset(l_e_offset,k_e_offset,size_e)
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call tce_e_offset(l_ex2_offset,k_ex2_offset,size_ex2)
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call tce_zero(d_e,size_e)
          call tce_filename('ex1',filename)
          call createfile(filename,d_ex1,size_ex1)
          call tce_zero(d_ex1,size_ex1)
          call tce_filename('ex2',filename)
          call createfile(filename,d_ex2,size_ex2)
          call tce_zero(d_ex2,size_ex2)
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
      if(perturbative.eq.'creom_(t)') then
      if(.not.slice) then
      call cr_eomccsd_t(xc1(ivec),k_x1_offset,xc2(ivec),k_x2_offset,
     1                        d_t1,k_t1_offset,d_t2,k_t2_offset,
     2                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                        d_e,k_e_offset,
     4                        d_ex1,k_ex1_offset,size_ex1,
     5                        d_ex2,k_ex2_offset,size_ex2,
     6                        d_c1,k_c1_offset,size_c1,
     7                        d_c2,k_c2_offset,size_c2,
     8                        excit,pt3_6,pt3_5,size_t1,size_x1,
     9                        r0xx)
      else
      if(.not.d4d5) then
      call cr_eomccsd_t_6dts(xc1(ivec),k_x1_offset,xc2(ivec),
     1                        k_x2_offset,
     1                        d_t1,k_t1_offset,d_t2,k_t2_offset,
     2                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                        d_e,k_e_offset,
     4                        d_ex1,k_ex1_offset,size_ex1,
     5                        d_ex2,k_ex2_offset,size_ex2,
     6                        d_c1,k_c1_offset,size_c1,
     7                        d_c2,k_c2_offset,size_c2,
     8        excit,pt3_6,pt3_5,size_t1,size_x1,r0xx,xmem)
      else
      call cr_eomccsd_t_6dts_d4d5_parallel(xc1(ivec),
     1                        k_x1_offset,xc2(ivec),
     1                        k_x2_offset,
     1                        d_t1,k_t1_offset,d_t2,k_t2_offset,
     2                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                        d_e,k_e_offset,
     4                        d_ex1,k_ex1_offset,size_ex1,
     5                        d_ex2,k_ex2_offset,size_ex2,
     6                        d_c1,k_c1_offset,size_c1,
     7                        d_c2,k_c2_offset,size_c2,
     8        excit,pt3_6,pt3_5,size_t1,size_x1,r0xx,xmem)
      end if
      end if
                    call gatodra(xc2(ivec))
                    call gatodra(xc1(ivec))
      elseif(perturbative.eq.'emb1') then
      call emb1(xc1(ivec),k_x1_offset,xc2(ivec),k_x2_offset,
     1                        d_t1,k_t1_offset,d_t2,k_t2_offset,
     2                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                        d_e,k_e_offset,
     4                        d_ex1,k_ex1_offset,size_ex1,
     5                        d_ex2,k_ex2_offset,size_ex2,
     6                        d_c1,k_c1_offset,size_c1,
     7                        d_c2,k_c2_offset,size_c2,
     8                        excit,pt3_6,pt3_5,size_t1,size_x1)
      else       !emb2
      call emb2(xc1(ivec),k_x1_offset,xc2(ivec),k_x2_offset,
     1                        d_t1,k_t1_offset,d_t2,k_t2_offset,
     2                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                        d_e,k_e_offset,
     4                        d_ex1,k_ex1_offset,size_ex1,
     5                        d_ex2,k_ex2_offset,size_ex2,
     6                        d_c1,k_c1_offset,size_c1,
     7                        d_c2,k_c2_offset,size_c2,
     8                        excit,pt3_6,pt3_5,size_t1,size_x1)
      end if
!
          call deletefile(d_ex2)
          call deletefile(d_ex1)
          call deletefile(d_e)
          if (.not.ma_pop_stack(l_ex2_offset))
     1     call errquit("tce_energy: MA problem",24,MA_ERR)
          if (.not.ma_pop_stack(l_ex1_offset))
     1     call errquit("tce_energy: MA problem",24,MA_ERR)
          if (.not.ma_pop_stack(l_e_offset))
     1     call errquit("tce_energy: MA problem",24,MA_ERR)
!
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
       if(perturbative.eq.'creom_(t)') then
          if (nodezero) then
           write(LuOut,9080)"CR-EOMCCSD(T) ",ref+corr+excit+pt3_6
! UWAGA
! ACHTUNG
! BHUMAHE
! ATTENTION
!           write(LuOut,9085)"CR-EOMCCSD(T)IA - CR-CCSD(T)IA",
!     &                      (ref+corr+excit+pt3_6-energy_grcr)*au2ev
           write(LuOut,9085)"CR-EOMCCSD(T) ",
     &                      (excit+pt3_6)*au2ev
           write(LuOut,*)
           call util_flush(LuOut)
          endif
         end if
       if(perturbative.eq.'emb1') then
          if (nodezero) then
           write(LuOut,9080)"r-CR-EOMCCSD(T)-I",ref+corr+excit+pt3_6
! UWAGA
! ACHTUNG
! BHUMAHE
! ATTENTION
           write(LuOut,9085)"delta-r-CR-EOMCCSD(T)-I        ",
     &                      (excit+pt3_6)*au2ev
           write(LuOut,*)
           call util_flush(LuOut)
          endif
         end if
       if(perturbative.eq.'emb2') then
          if (nodezero) then
           write(LuOut,9080)"r-CR-EOMCCSD(T)-II",ref+corr+excit+pt3_6
! UWAGA
! ACHTUNG
! BHUMAHE
! ATTENTION
           write(LuOut,9085)"delta-r-CR-EOMCCSD(T)-II        ",
     &                      (excit+pt3_6)*au2ev
           write(LuOut,*)
           call util_flush(LuOut)
          endif
         end if

! redefining the vector with excitation energies
          dbl_mb(k_omegax+ivec-1)=dbl_mb(k_omegax+ivec-1)+pt3_6
!
            call deletefile(d_c2)
            if (.not.ma_pop_stack(l_c2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
            call deletefile(d_c1)
            if (.not.ma_pop_stack(l_c1_offset))
     1        call errquit("tce_energy: MA problem",37,MA_ERR)
!
 3443   CONTINUE
        end if !cr-eom
!
!
!
!
!
        else if (model.eq.'ccsd_act') then
! == jaguar ====
      if(read_xa)  then
        nroots_reduced=nroots
        if(nodezero) then
         write(LuOut,*)'nroots_reduced: ',nroots_reduced
         call util_flush(LuOut)
        end if
        if (.not.rtdb_get(rtdb,'tce:omg',mt_dbl,nroots,eom_excit)) then
          call errquit('tce_energy: reading omg',nroots,CALC_ERR) 
        end if 
        if (.not.rtdb_get(rtdb,'tce:r0',mt_dbl,nroots,eom_r0)) then
          if(read_in3) then
           call errquit('tce: read_in3 active - r0 reading failed',
     &      nroots,CALC_ERR)
          else
           if(nodezero) then
           write(LuOut,*)'no R0-s restarted !!!!!'
           call util_flush(LuOut)
           end if 
          end if
        end if
        do ivec=1,nroots
        dbl_mb(k_omegax+ivec-1)=eom_excit(ivec)
        if(nodezero) then
         write(LuOut,*)'read excitation energy: ',
     1                  dbl_mb(k_omegax+ivec-1)
         call util_flush(LuOut)
        end if
        enddo
!
       if(nroots_reduced.ge.1) then
        unitn=81
        ivec=1
        call tce_filenameindexed(ivec,'xc1',filename)
        call createfile(filename,xc1(ivec),size_x1)
        xc1_exist(ivec) = .true.
        call util_file_name0('x1am1',.false.,.true.,filename,fldgts)
        call read_tensor(filename,xc1(ivec),size_x1,unitn)
        unitn=82
        call tce_filenameindexed(ivec,'xc2',filename)
        call createfile(filename,xc2(ivec),size_x2)
        xc2_exist(ivec) = .true.
        call util_file_name0('x2am1',.false.,.true.,filename,fldgts)
        call read_tensor(filename,xc2(ivec),size_x2,unitn)
          if(nodezero)
     &            write(LuOut,*)'EOMCCSD root nr. ',ivec
          call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
          call tce_print_x2a(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
!
        if(nroots_reduced.ge.2) then
         unitn=83
         ivec=2
         call tce_filenameindexed(ivec,'xc1',filename)
         call createfile(filename,xc1(ivec),size_x1)
         xc1_exist(ivec) = .true.
         call util_file_name0('x1am2',.false.,.true.,filename,fldgts)
         call read_tensor(filename,xc1(ivec),size_x1,unitn)
         unitn=84
         call tce_filenameindexed(ivec,'xc2',filename)
         call createfile(filename,xc2(ivec),size_x2)
         xc2_exist(ivec) = .true.
         call util_file_name0('x2am2',.false.,.true.,filename,fldgts)
         call read_tensor(filename,xc2(ivec),size_x2,unitn)
           if(nodezero)
     &         write(LuOut,*)'EOMCCSD root nr. ',ivec
           call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
           call tce_print_x2a(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
!
         if(nroots_reduced.ge.3) then
          unitn=85
          ivec=3
          call tce_filenameindexed(ivec,'xc1',filename)
          call createfile(filename,xc1(ivec),size_x1)
          xc1_exist(ivec) = .true.
          call util_file_name0('x1am3',.false.,.true.,filename,fldgts)
          call read_tensor(filename,xc1(ivec),size_x1,unitn)
          unitn=86
          call tce_filenameindexed(ivec,'xc2',filename)
          call createfile(filename,xc2(ivec),size_x2)
          xc2_exist(ivec) = .true.
          call util_file_name0('x2am3',.false.,.true.,filename,fldgts)
          call read_tensor(filename,xc2(ivec),size_x2,unitn)
            if(nodezero)
     &              write(LuOut,*)'EOMCCSD root nr. ',ivec
            call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
            call tce_print_x2a(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
!
          if(nroots_reduced.ge.4) then
           unitn=87
           ivec=4
           call tce_filenameindexed(ivec,'xc1',filename)
           call createfile(filename,xc1(ivec),size_x1)
           xc1_exist(ivec) = .true.
           call util_file_name0('x1am4',.false.,.true.,filename,fldgts)
           call read_tensor(filename,xc1(ivec),size_x1,unitn)
           unitn=88
           call tce_filenameindexed(ivec,'xc2',filename)
           call createfile(filename,xc2(ivec),size_x2)
           xc2_exist(ivec) = .true.
           call util_file_name0('x2am4',.false.,.true.,filename,fldgts)
           call read_tensor(filename,xc2(ivec),size_x2,unitn)
             if(nodezero)
     &              write(LuOut,*)'EOMCCSD root nr. ',ivec
             call tce_print_x1(xc1(ivec),k_x1_offset,0.10d0,irrep_x)
             call tce_print_x2a(xc2(ivec),k_x2_offset,0.10d0,irrep_x)
          end if
         end if
        end if
       end if
!
        go to 2001
      end if
! ==============

! -------------------------------------
!    ccsd_act/eomccsd_act
! -------------------------------------
!
       modelname = "CIS ini. guess for act-EOM-CCSD"
        if (nodezero) write(LuOut,9220)
     1  modelname(1:inp_strlen(modelname))
!
        call cis(LuOut,size_x1,k_x1_offset,d_f1,k_f1_offset,d_v2,
     1           k_v2_offset,k_omegax,l_omegax,d_rx1)
!
!
!
           if(nodezero) then
            write(LuOut,*)'INITIAL STARTS -----'
            call util_flush(LuOut)
           end if
           do iroot=1,nroots_reduced
            call tce_print_x1(x1_ini(iroot),k_x1_offset,0.10d0,irrep_x)
           enddo
           if(nodezero) then
            write(LuOut,*)'                    '
            write(LuOut,*)'INITIAL STARTS -----'
            call util_flush(LuOut)
           end if
!
!
!
            if (nxtrials.eq.0) goto 2001
            call tce_eom_init
!
      ister=0
! rec_s_1
       CALL OFFSET_eomccsd_x1a_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs1',-1,MA_ERR)
! rec_s_2
       CALL OFFSET_eomccsd_x1a_2_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_2_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs2',-1,MA_ERR)
! rec_s_3
       CALL OFFSET_eomccsd_x1a_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs3',-1,MA_ERR)
! rec_s_4
       CALL OFFSET_eomccsd_x1a_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs4',-1,MA_ERR)
! rec_s_5
       CALL OFFSET_eomccsd_x1a_7_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_7_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_s_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs5',-1,MA_ERR)
! rec_d_1
       CALL OFFSET_eomccsd_x2a_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd1',-1,MA_ERR)
! rec_d_2
       CALL OFFSET_eomccsd_x2a_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_3_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd2',-1,MA_ERR)
! rec_d_3
       CALL OFFSET_eomccsd_x2a_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd3',-1,MA_ERR)
! rec_d_4
       CALL OFFSET_eomccsd_x2a_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd4',-1,MA_ERR)
! rec_d_5
       CALL OFFSET_eomccsd_x2a_6_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_6_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd5',-1,MA_ERR)
! rec_d_6
       CALL OFFSET_eomccsd_x2a_8_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_6,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd6',-1,MA_ERR)
! rec_d_7
       CALL OFFSET_eomccsd_x2a_8_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_3_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_7,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd7',-1,MA_ERR)
! rec_d_8
       CALL OFFSET_eomccsd_x2a_8_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_4_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_8,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd8',-1,MA_ERR)
! rec_d_9
      CALL OFFSET_eomccsd_x2a_10_1_1(l_aux_off,k_aux_off,size_aux)
      CALL TCE_FILENAME('eomccsd_x2_10_1_1_i2',filename)
      CALL CREATEFILE(filename,rec_d_9,size_aux)
      IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd9',-1,MA_ERR)
!
            xdiis=diis
            do imicro=1,xdiis*nroots_reduced
              call tce_filenameindexed(imicro,'xp1',filename)
              call createfile(filename,xp1(imicro),size_x1)
              call gatodra(xp1(imicro))
              xp1_exist(imicro)=.true.
              call tce_filenameindexed(imicro,'xp2',filename)
              call createfile(filename,xp2(imicro),size_x2)
              call gatodra(xp2(imicro))
              xp2_exist(imicro)=.true.
              call tce_filenameindexed(imicro,'x1',filename)
              call createfile(filename,x1(imicro),size_x1)
              call tce_zero(x1(imicro),size_x1)
              call gatodra(x1(imicro))
              x1_exist(imicro)=.true.
              call tce_filenameindexed(imicro,'x2',filename)
              call createfile(filename,x2(imicro),size_x2)
              call tce_zero(x2(imicro),size_x2)
              call gatodra(x2(imicro))
              x2_exist(imicro)=.true.
            enddo
!
          do ivec=1,nroots_reduced
            call tce_filenameindexed(ivec,'xc1',filename)
            call createfile(filename,xc1(ivec),size_x1)
            call gatoeaf(xc1(ivec))
            call gatodra(xc1(ivec))
            xc1_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc2',filename)
            call createfile(filename,xc2(ivec),size_x2)
            call gatoeaf(xc2(ivec))
            call gatodra(xc2(ivec))
            xc2_exist(ivec) = .true.
          enddo
!
            modelname = "EOM-CCSD-act"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
!
! open k_h_eff  and k_c_eff
!
        if (.not.MA_PUSH_GET(mt_dbl,
     1      xdiis*xdiis*nroots_reduced*nroots_reduced,'h_eff',
     1      l_h_eff,k_h_eff)) call errquit('h_eff',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,xdiis*nroots_reduced*nroots_reduced,
     1      'c_eff',
     1      l_c_eff,k_c_eff)) call errquit('c_eff',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,nroots_reduced,'residual',
     1          l_residual,k_residual)) call errquit('l_res',1,MA_ERR)
!
!mmm          do iroot=1,nroots_reduced  ! ---- main loop over the roots
!
            iter=1
!
!
             if(iter.eq.1) then
              do iroot=1,nroots_reduced
               call copyfile(x1_ini(iroot),xc1(iroot),size_x1)
               call tce_zero(xc2(iroot),size_x2)
              enddo
             endif 
!
!
 4111         continue  !iter continuation
!
             imicro = 1
!
            do iroot=1,nroots_reduced
             call copyfile(xc1(iroot),x1(iroot),size_x1)
             call copyfile(xc2(iroot),x2(iroot),size_x2)
            enddo
!
            call norm_x_vector_act_mr(k_x1_offset,k_x2_offset,
     1                            size_x1,size_x2,nroots_reduced)
!
            do iroot=1,nroots_reduced
             call tce_zero(xp1(iroot),size_x1)
             call tce_zero(xp2(iroot),size_x2)
            enddo
!
             imicro = 1
!
             do i=1,xdiis*nroots_reduced*nroots_reduced
              dbl_mb(k_c_eff+i-1)=0.0d0
             enddo
             do i=1,nroots_reduced
             dbl_mb(k_c_eff+(i-1)*xdiis*nroots_reduced+i-1)=1.0d0
             enddo
             do i=1,xdiis*xdiis*nroots_reduced*nroots_reduced
               dbl_mb(k_h_eff+i-1)=0.0d0
             enddo
!
      cpu=-util_cpusec()
      wall=-util_wallsec()
!
      do iroot=1,nroots_reduced ! x-x-x-x-x-x-x-x-x
!
! - T1/X1 LOCALIZATION ----------
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,size_x1,'x1_local',
     1      l_x1_local,k_x1_local))
     1      call errquit('x1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
        call ma_zero(dbl_mb(k_x1_local),size_x1)
        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
        call ga_get(x1(iroot),1,size_x1,1,1,dbl_mb(k_x1_local),1)
! -------------------------------
          call eomccsd_x1a(d_f1,xp1(iroot),k_t1_local,d_t2,d_v2,
     1      k_x1_local,x2(iroot),k_f1_offset,k_x1_offset,
     2      k_t1_offset,k_t2_offset,k_v2_offset,
     3      k_x1_offset,k_x2_offset,
     4      ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5      rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6      rec_d_6,rec_d_7,rec_d_8,rec_d_9)
          call eomccsd_x2a(d_f1,xp2(iroot),k_t1_local,d_t2,d_v2,
     1       k_x1_local,x2(iroot),k_f1_offset,k_x2_offset,
     2       k_t1_offset,k_t2_offset,k_v2_offset,
     3       k_x1_offset,k_x2_offset,size_x1,size_x2,
     4       ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9)
!
         ister=1
! - T1/X1 LOCALIZATION --
         if(.not.MA_POP_STACK(l_x1_local))
     &      call errquit('l_x1_local',4,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
! -----------------------------
         enddo ! x-x-x-x-x-x-x-x-x
!   calculate  H-bar delta x
!
         call form_heff_act_mr(k_x1_offset,k_x2_offset,size_x1,size_x2,
     &        dbl_mb(k_c_eff),dbl_mb(k_h_eff),dbl_mb(k_omegax),
     &        imicro,xdiis,iter)
!
          do imicro=2,xdiis !microloop(imicro)
!
          do iroot=1,nroots_reduced
           ipos=(imicro-1)*nroots_reduced+iroot
           call tce_zero(x1(ipos),size_x1)
           call tce_zero(x2(ipos),size_x2)
           call tce_zero(xp1(ipos),size_x1)
           call tce_zero(xp2(ipos),size_x2)
          enddo
!
          do iroot=1,nroots_reduced
           dbl_mb(k_residual+iroot-1)=0.0d0
          enddo
!
          do iroot=1,nroots_reduced
           call form_correction_act_mr(k_x1_offset,k_x2_offset,
     1                  size_x1,size_x2,
     1                  dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     1                  dbl_mb(k_omegax+iroot-1),
     1                  imicro,xdiis,dbl_mb(k_residual+iroot-1),iroot)
          enddo
!
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
      if(nodezero) then
                  write(LuOut,9211) iter
             do iroot=1,nroots_reduced
                  write(LuOut,9230) dbl_mb(k_residual+iroot-1),
     2            dbl_mb(k_omegax+iroot-1),
     3            dbl_mb(k_omegax+iroot-1)*au2ev,cpu,wall
             enddo
      end if
      call util_flush(LuOut)
      cpu=-util_cpusec()
      wall=-util_wallsec()
!
           converged=.true.
           do iroot=1,nroots_reduced-1
             if(dabs(dbl_mb(k_residual+iroot-1)).gt.thresh_eom)
     &             converged=.false.
           enddo
           if(dabs(dbl_mb(k_residual+nroots_reduced-1)).gt.thresh_l)
     7             converged=.false.
!
            if(converged) then
             do iroot=1,nroots_reduced
              call form_best_guess_act_mr(k_x1_offset,k_x2_offset,
     &                size_x1,size_x2,
     &                dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     &                imicro-1,
     &                xdiis,iroot)
             enddo
!
             do iroot=1,nroots_reduced
              if(nodezero)
     &                write(LuOut,*)'EOMCCSD root nr. ',iroot
              call tce_print_x1(xc1(iroot),k_x1_offset,0.10d0,irrep_x)
              call tce_print_x2a(xc2(iroot),k_x2_offset,0.10d0,irrep_x)
             enddo
             go to 4112
            end if
! all together
            call deltax_ort_norm_act_mr(k_x1_offset,k_x2_offset,
     1                              size_x1,size_x2,
     1                              imicro,xdiis)
!
         do iroot=1,nroots_reduced ! --- iroot
!
          ipos=(imicro-1)*nroots_reduced+iroot
!
! - T1/X1 LOCALIZATION ----------
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,size_x1,'x1_local',
     1      l_x1_local,k_x1_local))
     1      call errquit('x1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
        call ma_zero(dbl_mb(k_x1_local),size_x1)
        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
        call ga_get(x1(ipos),1,size_x1,1,1,dbl_mb(k_x1_local),1)
! -------------------------------
          call eomccsd_x1a(d_f1,xp1(ipos),k_t1_local,d_t2,d_v2,
     1      k_x1_local,x2(ipos),k_f1_offset,k_x1_offset,
     2      k_t1_offset,k_t2_offset,k_v2_offset,
     3      k_x1_offset,k_x2_offset,
     4      ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5      rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6      rec_d_6,rec_d_7,rec_d_8,rec_d_9)
          call eomccsd_x2a(d_f1,xp2(ipos),k_t1_local,d_t2,d_v2,
     1       k_x1_local,x2(ipos),k_f1_offset,k_x2_offset,
     2       k_t1_offset,k_t2_offset,k_v2_offset,
     3       k_x1_offset,k_x2_offset,size_x1,size_x2,
     4       ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9)
!
         ister=1
! - T1/X1 LOCALIZATION --
         if(.not.MA_POP_STACK(l_x1_local))
     &      call errquit('l_x1_local',4,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
! -----------------------
!
        enddo ! --- iroot 
!
!             calculate  H-bar delta x
!
         call form_heff_act_mr(k_x1_offset,k_x2_offset,size_x1,size_x2,
     &        dbl_mb(k_c_eff),dbl_mb(k_h_eff),dbl_mb(k_omegax),
     &        imicro,xdiis,iter)
!
           iter=iter+1                              
             if(iter.gt.maxiter) then
              if(nodezero) then
                 write(LuOut,9233) iroot
              end if
!
             do iroot=1,nroots_reduced
              call form_best_guess_act_mr(k_x1_offset,k_x2_offset,
     &                size_x1,size_x2,
     &                dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     &                imicro,
     &                xdiis,iroot)
              call tce_print_x1(xc1(iroot),k_x1_offset,0.10d0,irrep_x)
              call tce_print_x2a(xc2(iroot),k_x2_offset,0.10d0,irrep_x)
             enddo
!
              go to 4112
             end if
!
          enddo   !microloop(imicro)
!
             do iroot=1,nroots_reduced
              call form_best_guess_act_mr(k_x1_offset,k_x2_offset,
     &                size_x1,size_x2,
     &                dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     &                xdiis,
     &                xdiis,iroot)
             enddo
!
               go to 4111
!
 4112          continue !from converged
!mmm            enddo ! ---- main loop over the roots (iroot)
!
! === jaguar ===
      if(write_xa)  then
       if(nroots_reduced.ge.1) then
        unitn=81
        call util_file_name0('x1am1',.false.,.true.,filename,fldgts)
        call write_tensor(filename,xc1(1),size_x1,unitn)
        unitn=82
        call util_file_name0('x2am1',.false.,.true.,filename,fldgts)
        call write_tensor(filename,xc2(1),size_x2,unitn)
        if(nroots_reduced.ge.2) then
         unitn=83
         call util_file_name0('x1am2',.false.,.true.,filename,fldgts)
         call write_tensor(filename,xc1(2),size_x1,unitn)
         unitn=84
         call util_file_name0('x2am2',.false.,.true.,filename,fldgts)
         call write_tensor(filename,xc2(2),size_x2,unitn)
         if(nroots_reduced.ge.3) then
          unitn=85
          call util_file_name0('x1am3',.false.,.true.,filename,fldgts)
          call write_tensor(filename,xc1(3),size_x1,unitn)
          unitn=86
          call util_file_name0('x2am3',.false.,.true.,filename,fldgts)
          call write_tensor(filename,xc2(3),size_x2,unitn)
          if(nroots_reduced.ge.4) then
           unitn=87
           call util_file_name0('x1am4',.false.,.true.,filename,fldgts)
           call write_tensor(filename,xc1(4),size_x1,unitn)
           unitn=88
           call util_file_name0('x2am4',.false.,.true.,filename,fldgts)
           call write_tensor(filename,xc2(4),size_x2,unitn)
          end if
         end if
        end if
       end if
      end if
! ==============
!
         if(.not.MA_POP_STACK(l_residual))
     &      call errquit('l_residual',4,MA_ERR)
         if(.not.MA_POP_STACK(l_c_eff))
     &      call errquit('l_c_eff',4,MA_ERR)
         if(.not.MA_POP_STACK(l_h_eff))
     &      call errquit('l_h_eff',4,MA_ERR)
!
            do imicro=1,xdiis
              call deletefile(xp1(imicro))
              xp1_exist(imicro)=.false.
              call deletefile(xp2(imicro))
              xp2_exist(imicro)=.false.
              call deletefile(x1(imicro))
              x1_exist(imicro)=.false.
              call deletefile(x2(imicro))
              x2_exist(imicro)=.false.
            enddo
!
!
!
! xc1 xc2 are deleted by tce_eom_tidy (later on)
!
!
!
          call deletefile(rec_d_9)
          call deletefile(rec_d_8)
          call deletefile(rec_d_7)
          call deletefile(rec_d_6)
          call deletefile(rec_d_5)
          call deletefile(rec_d_4)
          call deletefile(rec_d_3)
          call deletefile(rec_d_2)
          call deletefile(rec_d_1)
          call deletefile(rec_s_5)
          call deletefile(rec_s_4)
          call deletefile(rec_s_3)
          call deletefile(rec_s_2)
          call deletefile(rec_s_1)
!
         do ivec=1,maxtrials
           if(x1_ini_exist(ivec))  then
            call deletefile(x1_ini(ivec))
            x1_ini_exist(ivec) = .false.
           end if
         enddo
!
!
! --------------------------------------
!
 2001   continue
! -------------------------
! CR-EOMCCSD(T)-ACT
! -------------------------
        if ((perturbative.eq.'creom(t)a').or.
     &      (perturbative.eq.'creom(t)asim')) THEN
        irrep_c = irrep_x  !important
        DO 3449,ivec = 1,nroots_reduced
          cpu = - util_cpusec()
          wall = - util_wallsec()
! d_e,d_ex1,d_ex2 global arrays
!
            call tce_x1_offset(l_c1_offset,k_c1_offset,size_c1)
            call tce_filename('c1',filename)
            call createfile(filename,d_c1,size_c1)
!c            call gatoeaf(d_c1)
            call tce_zero(d_c1,size_c1)
!
            call tce_x2a_offset(l_c2_offset,k_c2_offset,size_c2)
            call tce_filename('c2',filename)
            call createfile(filename,d_c2,size_c2)
!c            call gatoeaf(d_c2)
            call tce_zero(d_c2,size_c2)
!
      pt3_5 = 0.0d0
      pt3_6 = 0.0d0
      excit = dbl_mb(k_omegax+ivec-1)
      call tce_e_offset(l_e_offset,k_e_offset,size_e)
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call tce_e_offset(l_ex2_offset,k_ex2_offset,size_ex2)
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call tce_zero(d_e,size_e)
          call tce_filename('ex1',filename)
          call createfile(filename,d_ex1,size_ex1)
          call tce_zero(d_ex1,size_ex1)
          call tce_filename('ex2',filename)
          call createfile(filename,d_ex2,size_ex2)
          call tce_zero(d_ex2,size_ex2)
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
      if(perturbative.eq.'creom(t)a') then
      call cr_eomccsd_t_act(
     1              xc1(ivec),k_x1_offset,xc2(ivec),k_x2_offset,
     1                        d_t1,k_t1_offset,d_t2,k_t2_offset,
     2                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                        d_e,k_e_offset,
     4                        d_ex1,k_ex1_offset,size_ex1,
     5                        d_ex2,k_ex2_offset,size_ex2,
     6                        d_c1,k_c1_offset,size_c1,
     7                        d_c2,k_c2_offset,size_c2,
     8                        excit,pt3_6,pt3_5,size_t1,size_x1)
                    call gatodra(xc2(ivec))
                    call gatodra(xc1(ivec))
      end if
      if(perturbative.eq.'creom(t)asim') then
      call cr_eomccsd_tasim(
     1              xc1(ivec),k_x1_offset,xc2(ivec),k_x2_offset,
     1                        d_t1,k_t1_offset,d_t2,k_t2_offset,
     2                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                        d_e,k_e_offset,
     4                        d_ex1,k_ex1_offset,size_ex1,
     5                        d_ex2,k_ex2_offset,size_ex2,
     6                        d_c1,k_c1_offset,size_c1,
     7                        d_c2,k_c2_offset,size_c2,
     8                        excit,pt3_6,pt3_5,size_t1,size_x1)
                    call gatodra(xc2(ivec))
                    call gatodra(xc1(ivec))
      end if
!
          call deletefile(d_ex2)
          call deletefile(d_ex1)
          call deletefile(d_e)
          if (.not.ma_pop_stack(l_ex2_offset))
     1     call errquit("tce_energy: MA problem",24,MA_ERR)
          if (.not.ma_pop_stack(l_ex1_offset))
     1     call errquit("tce_energy: MA problem",24,MA_ERR)
          if (.not.ma_pop_stack(l_e_offset))
     1     call errquit("tce_energy: MA problem",24,MA_ERR)
!
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
         if(perturbative.eq.'creom(t)a') then
          if (nodezero) then
          write(LuOut,9080)"CR-EOMCCSD(T)-ACT(II) ",ref+corr+excit+pt3_6
          write(LuOut,9085)"CR-EOMCCSD(T)-ACT(II) ",
     &                      (excit+pt3_6)*au2ev
           write(LuOut,*)
           call util_flush(LuOut)
          endif
         end if
         if(perturbative.eq.'creom(t)asim') then
          if (nodezero) then
      write(LuOut,9080)"CR-EOMCCSD(T)-ACT(II)-SIMP",ref+corr+excit+pt3_6
      write(LuOut,9085)"CR-EOMCCSD(T)-ACT(II)-SIMP",
     &                      (excit+pt3_6)*au2ev
           write(LuOut,*)
           call util_flush(LuOut)
          endif
         end if
! redefining the vector with excitation energies
          dbl_mb(k_omegax+ivec-1)=dbl_mb(k_omegax+ivec-1)+pt3_6
!
            call deletefile(d_c2)
            if (.not.ma_pop_stack(l_c2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
            call deletefile(d_c1)
            if (.not.ma_pop_stack(l_c1_offset))
     1        call errquit("tce_energy: MA problem",37,MA_ERR)
!
 3449   CONTINUE
        end if !triples corr. to EOMCCSd
!
!
!
!
!
          elseif (model.eq.'ccsdt') then
#if !defined(CCSDT)
      CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
! -------------
! EOM-CCSDT   right
! -------------
            call tce_eom_xguess(rtdb,needt1,needt2,needt3,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
!edo            if (nxtrials.eq.0) goto 210
            if (nxtrials.eq.0) then
               ccsdtleft=.true.
               goto 400
            endif
!cxsg
          do ivec=1,nroots_reduced
            call tce_filenameindexed(ivec,'xc1',filename)
            call createfile(filename,xc1(ivec),size_x1)
            call gatoeaf(xc1(ivec))
            call gatodra(xc1(ivec))
            xc1_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc2',filename)
            call createfile(filename,xc2(ivec),size_x2)
            call gatoeaf(xc2(ivec))
            call gatodra(xc2(ivec))
            xc2_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc3',filename)
            call createfile(filename,xc3(ivec),size_x3)
            call gatoeaf(xc3(ivec))
            call gatodra(xc3(ivec))
            xc3_exist(ivec) = .true.
          enddo
!cxsg
            modelname = "EOM-CCSDT right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
!           call hbar_init
            do iter=1,maxiter
              if (nodezero.and.util_print('eom',print_default))
     1          write(LuOut,9210) iter,nxtrials
              do ivec = 1,nxtrials
                if (.not.xp1_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp1',filename)
                  call createfile(filename,xp1(ivec),size_x1)
                  xp1_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call eomccsdt_x1(d_f1,xp1(ivec),d_t1,d_t2,d_v2,
     1              x1(ivec),x2(ivec),x3(ivec),k_f1_offset,k_x1_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,k_x1_offset,
     3              k_x2_offset,k_x3_offset)
                  call reconcilefile(xp1(ivec),size_x1)
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp1(ivec))
                endif
                if (.not.xp2_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp2',filename)
                  call createfile(filename,xp2(ivec),size_x2)
                  xp2_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call eomccsdt_x2(d_f1,xp2(ivec),d_t1,d_t2,d_t3,d_v2,
     1              x1(ivec),x2(ivec),x3(ivec),k_f1_offset,k_x2_offset,
     2              k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     3              k_x1_offset,k_x2_offset,k_x3_offset)
                  call reconcilefile(xp2(ivec),size_x2)
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp2(ivec))
                endif
                if (.not.xp3_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp3',filename)
                  call createfile(filename,xp3(ivec),size_x3)
                  xp3_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call eomccsdt_x3(d_f1,xp3(ivec),d_t1,d_t2,d_t3,d_v2,
     1              x1(ivec),x2(ivec),x3(ivec),k_f1_offset,k_x3_offset,
     2              k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     3              k_x1_offset,k_x2_offset,k_x3_offset)
                  call reconcilefile(xp3(ivec),size_x3)
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp3(ivec))
                endif
              enddo
              if (.not.ma_push_get(mt_dbl,nxtrials,'residual',
     1          l_residual,k_residual))
     2          call errquit('tce_energy: MA problem',101,MA_ERR)
              call tce_eom_xdiagon(needt1,needt2,needt3,needt4,
     1          size_x1,size_x2,size_x3,size_x4,
     2          k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
     3          d_rx1,d_rx2,d_rx3,d_rx4,
     4          dbl_mb(k_omegax),dbl_mb(k_residual),k_hbar,iter,
     +          eaccsd,ipccsd)
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
              converged = .true.
              do ivec = 1,nroots_reduced
                if (nodezero.and.(ivec.ne.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegax+ivec-1),
     3            dbl_mb(k_omegax+ivec-1)*au2ev
                if (nodezero.and.(ivec.eq.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegax+ivec-1),
     3            dbl_mb(k_omegax+ivec-1)*au2ev,cpu,wall
                if (nodezero) call util_flush(LuOut)
                if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1            converged = .false.
              enddo
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (.not.ma_pop_stack(l_residual))
     1          call errquit("tce_energy: MA problem",102,MA_ERR)
              if (converged) then
                call tce_eom_xtidy
!               call hbar_tidy
                if (nodezero) then
                  write(LuOut,9240)
                  call util_flush(LuOut)
                endif
! ---
            if(nodezero) then
             write(LuOut,*)'largest EOMCCSDT amplitudes: R1 and R2'
             call util_flush(LuOut)
            end if
             do jvec=1,nroots_reduced
!              call tce_print_x1(xc1(jvec),k_x1_offset,0.10d0)
!              call tce_print_x2(xc2(jvec),k_x2_offset,0.10d0)
              call tce_print_x1(xc1(jvec),k_x1_offset,printtol,irrep_x)
              call tce_print_x2(xc2(jvec),k_x2_offset,printtol,irrep_x)
              call util_flush(LuOut)
             enddo
! ---
                ccsdtleft=.true.
                goto 400
!edo                goto 210
              endif
            enddo
            call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
#endif
!
!
!
!
          elseif (model.eq.'ccsdta') then
#if !defined(CCSDT)
      CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
! -----------------------
! EOM-CCSDt  right
! -----------------------
            call tce_eom_xguess(rtdb,needt1,needt2,needt3a,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
!ccx            if (nxtrials.eq.0) goto 210
!cxsg
          do ivec=1,nroots_reduced
            call tce_filenameindexed(ivec,'xc1',filename)
            call createfile(filename,xc1(ivec),size_x1)
            call gatoeaf(xc1(ivec))
            call gatodra(xc1(ivec))
            xc1_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc2',filename)
            call createfile(filename,xc2(ivec),size_x2)
            call gatoeaf(xc2(ivec))
            call gatodra(xc2(ivec))
            xc2_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc3',filename)
            call createfile(filename,xc3(ivec),size_x3)
            call gatoeaf(xc3(ivec))
            call gatodra(xc3(ivec))
            xc3_exist(ivec) = .true.
          enddo
!cxsg
            modelname = "EOM-CCSDt right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
!           call hbar_init
            do iter=1,maxiter
              if (nodezero.and.util_print('eom',print_default))
     1          write(LuOut,9210) iter,nxtrials
              do ivec = 1,nxtrials
                if (.not.xp1_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp1',filename)
                  call createfile(filename,xp1(ivec),size_x1)
                  xp1_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call eomccsdt_x1a(d_f1,xp1(ivec),d_t1,d_t2,d_v2,
     1              x1(ivec),x2(ivec),x3(ivec),k_f1_offset,k_x1_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,k_x1_offset,
     3              k_x2_offset,k_x3_offset)
                  call reconcilefile(xp1(ivec),size_x1)
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp1(ivec))
                endif
                if (.not.xp2_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp2',filename)
                  call createfile(filename,xp2(ivec),size_x2)
                  xp2_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call eomccsdt_x2a(d_f1,xp2(ivec),d_t1,d_t2,d_t3,d_v2,
     1              x1(ivec),x2(ivec),x3(ivec),k_f1_offset,k_x2_offset,
     2              k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     3              k_x1_offset,k_x2_offset,k_x3_offset)
                  call reconcilefile(xp2(ivec),size_x2)
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp2(ivec))
                endif
                if (.not.xp3_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp3',filename)
                  call createfile(filename,xp3(ivec),size_x3)
                  xp3_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call eomccsdt_x3a(d_f1,xp3(ivec),d_t1,d_t2,d_t3,d_v2,
     1              x1(ivec),x2(ivec),x3(ivec),k_f1_offset,k_x3_offset,
     2              k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     3              k_x1_offset,k_x2_offset,k_x3_offset)
                  call reconcilefile(xp3(ivec),size_x3)
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp3(ivec))
                endif
              enddo
              if (.not.ma_push_get(mt_dbl,nxtrials,'residual',
     1          l_residual,k_residual))
     2          call errquit('tce_energy: MA problem',101,MA_ERR)
              call tce_eom_xdiagon(needt1,needt2,needt3a,needt4,
     1          size_x1,size_x2,size_x3,size_x4,
     2          k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
     3          d_rx1,d_rx2,d_rx3,d_rx4,
     4          dbl_mb(k_omegax),dbl_mb(k_residual),k_hbar,iter,
     +          eaccsd,ipccsd)
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
              converged = .true.
              do ivec = 1,nroots_reduced
                if (nodezero.and.(ivec.ne.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegax+ivec-1),
     3            dbl_mb(k_omegax+ivec-1)*au2ev
                if (nodezero.and.(ivec.eq.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegax+ivec-1),
     3            dbl_mb(k_omegax+ivec-1)*au2ev,cpu,wall
                if (nodezero) call util_flush(LuOut)
                if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1            converged = .false.
              enddo
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (.not.ma_pop_stack(l_residual))
     1          call errquit("tce_energy: MA problem",102,MA_ERR)
              if (converged) then
                call tce_eom_xtidy
!               call hbar_tidy
                if (nodezero) then
                  write(LuOut,9240)
                  call util_flush(LuOut)
                endif
! ---
            if(nodezero) then
             write(LuOut,*)'largest EOMCCSDt amplitudes: R1 and R2 '
             call util_flush(LuOut)
            end if
             do jvec=1,nroots_reduced
!              call tce_print_x1(xc1(jvec),k_x1_offset,0.10d0)
!              call tce_print_x2(xc2(jvec),k_x2_offset,0.10d0)
              call tce_print_x1(xc1(jvec),k_x1_offset,printtol,irrep_x)
              call tce_print_x2(xc2(jvec),k_x2_offset,printtol,irrep_x)
              call util_flush(LuOut)
             enddo
! ---
                ccsdtleft=.true.
                goto 400
!edo                goto 210
              endif
            enddo
            call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
#endif
!
!
!
!
          elseif (model.eq.'ccsdtq') then
#if !defined(CCSDTQ)
      CALL ERRQUIT('CCSDTQ not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
! -------------
! EOM-CCSDTQ  right
! -------------
            call tce_eom_xguess(rtdb,needt1,needt2,needt3,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
            if (nxtrials.eq.0) goto 220
!cxsg
          do ivec=1,nroots_reduced
            call tce_filenameindexed(ivec,'xc1',filename)
            call createfile(filename,xc1(ivec),size_x1)
            call gatoeaf(xc1(ivec))
            call gatodra(xc1(ivec))
            xc1_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc2',filename)
            call createfile(filename,xc2(ivec),size_x2)
            call gatoeaf(xc2(ivec))
            call gatodra(xc2(ivec))
            xc2_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc3',filename)
            call createfile(filename,xc3(ivec),size_x3)
            call gatoeaf(xc3(ivec))
            call gatodra(xc3(ivec))
            xc3_exist(ivec) = .true.
            call tce_filenameindexed(ivec,'xc4',filename)
            call createfile(filename,xc4(ivec),size_x4)
            call gatoeaf(xc4(ivec))
            call gatodra(xc4(ivec))
            xc4_exist(ivec) = .true.
          enddo
!cxsg
            modelname = "EOM-CCSDTQ right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
!           call hbar_init
            do iter=1,maxiter
              if (nodezero.and.util_print('eom',print_default))
     1          write(LuOut,9210) iter,nxtrials
              do ivec = 1,nxtrials
                if (.not.xp1_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp1',filename)
                  call createfile(filename,xp1(ivec),size_x1)
                  xp1_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call dratoga(x4(ivec))
                  call eomccsdtq_x1(d_f1,xp1(ivec),d_t1,d_t2,d_v2,
     1              x1(ivec),x2(ivec),x3(ivec),k_f1_offset,k_x1_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,k_x1_offset,
     3              k_x2_offset,k_x3_offset)
                  call reconcilefile(xp1(ivec),size_x1)
                  call gatodra(x4(ivec))
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp1(ivec))
                endif
                if (.not.xp2_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp2',filename)
                  call createfile(filename,xp2(ivec),size_x2)
                  xp2_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call dratoga(x4(ivec))
                  call eomccsdtq_x2(d_f1,xp2(ivec),d_t1,d_t2,d_t3,d_v2,
     1              x1(ivec),x2(ivec),x3(ivec),x4(ivec),k_f1_offset,
     2              k_x2_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     3              k_v2_offset,k_x1_offset,k_x2_offset,k_x3_offset,
     4              k_x4_offset)
                  call reconcilefile(xp2(ivec),size_x2)
                  call gatodra(x4(ivec))
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp2(ivec))
                endif
                if (.not.xp3_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp3',filename)
                  call createfile(filename,xp3(ivec),size_x3)
                  xp3_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call dratoga(x4(ivec))
                  call eomccsdtq_x3(d_f1,xp3(ivec),d_t1,d_t2,d_t3,d_t4,
     1              d_v2,x1(ivec),x2(ivec),x3(ivec),x4(ivec),
     2              k_f1_offset,k_x3_offset,k_t1_offset,k_t2_offset,
     3              k_t3_offset,k_t4_offset,k_v2_offset,k_x1_offset,
     4              k_x2_offset,k_x3_offset,k_x4_offset)
                  call reconcilefile(xp3(ivec),size_x3)
                  call gatodra(x4(ivec))
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp3(ivec))
                endif
                if (.not.xp4_exist(ivec)) then
                  call tce_filenameindexed(ivec,'xp4',filename)
                  call createfile(filename,xp4(ivec),size_x4)
                  xp4_exist(ivec) = .true.
                  call dratoga(x1(ivec))
                  call dratoga(x2(ivec))
                  call dratoga(x3(ivec))
                  call dratoga(x4(ivec))
                  call eomccsdtq_x4(d_f1,xp4(ivec),d_t1,d_t2,d_t3,d_t4,
     1              d_v2,x1(ivec),x2(ivec),x3(ivec),x4(ivec),
     2              k_f1_offset,k_x4_offset,k_t1_offset,k_t2_offset,
     3              k_t3_offset,k_t4_offset,k_v2_offset,k_x1_offset,
     4              k_x2_offset,k_x3_offset,k_x4_offset)
                  call reconcilefile(xp4(ivec),size_x4)
                  call gatodra(x4(ivec))
                  call gatodra(x3(ivec))
                  call gatodra(x2(ivec))
                  call gatodra(x1(ivec))
                  call gatodra(xp4(ivec))
                endif
              enddo
              if (.not.ma_push_get(mt_dbl,nxtrials,'residual',
     1          l_residual,k_residual))
     2          call errquit('tce_energy: MA problem',101,MA_ERR)
              call tce_eom_xdiagon(needt1,needt2,needt3,needt4,
     1          size_x1,size_x2,size_x3,size_x4,
     2          k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
     3          d_rx1,d_rx2,d_rx3,d_rx4,
     4          dbl_mb(k_omegax),dbl_mb(k_residual),k_hbar,iter,
     +          eaccsd,ipccsd)
	      cpu=cpu+util_cpusec()
	      wall=wall+util_wallsec()
	      converged = .true.
	      do ivec = 1,nroots_reduced
		if (nodezero.and.(ivec.ne.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegax+ivec-1),
     3            dbl_mb(k_omegax+ivec-1)*au2ev
		if (nodezero.and.(ivec.eq.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegax+ivec-1),
     3            dbl_mb(k_omegax+ivec-1)*au2ev,cpu,wall
                if (nodezero) call util_flush(LuOut)
		if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1            converged = .false.
	      enddo
	      cpu=-util_cpusec()
	      wall=-util_wallsec()
	      if (.not.ma_pop_stack(l_residual))
     1          call errquit("tce_energy: MA problem",102,MA_ERR)
	      if (converged) then
		call tce_eom_xtidy
!               call hbar_tidy
		if (nodezero) then
		  write(LuOut,9240)
		  call util_flush(LuOut)
		endif
		goto 220
	      endif
	    enddo
	    call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
! -------------
! EOM-CCSDTQ  left
! -------------
  220       continue
	    if (left) then
	    call tce_eom_yguess(needt1,needt2,needt3,needt4,
     1        size_y1,size_y2,size_y3,size_y4,
     2        k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset)
	    if (nytrials.eq.0) goto 400
	    modelname = "EOM-CCSDTQ left-hand side"
	    if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
!           call hbar_init
	    do iter=1,maxiter
	      if (nodezero.and.util_print('eom',print_default))
     1          write(LuOut,9210) iter,nytrials
	      do ivec = 1,nytrials
		if (.not.yp1_exist(ivec)) then
		  call tce_filenameindexed(ivec,'yp1',filename)
		  call createfile(filename,yp1(ivec),size_y1)
		  yp1_exist(ivec) = .true.
		  call dratoga(y1(ivec))
		  call dratoga(y2(ivec))
		  call dratoga(y3(ivec))
		  call dratoga(y4(ivec))
!ccx                  call eomccsdtq_y1(d_f1,yp1(ivec),d_t1,d_t2,d_t3,d_t4,
!ccx     1              d_v2,y1(ivec),y2(ivec),y3(ivec),y4(ivec),
!ccx     2              k_f1_offset,k_y1_offset,k_t1_offset,k_t2_offset,
!ccx     3              k_t3_offset,k_t4_offset,k_v2_offset,k_y1_offset,
!ccx     4              k_y2_offset,k_y3_offset,k_y4_offset)
		  call reconcilefile(yp1(ivec),size_y1)
		  call gatodra(y4(ivec))
		  call gatodra(y3(ivec))
		  call gatodra(y2(ivec))
		  call gatodra(y1(ivec))
		  call gatodra(yp1(ivec))
		endif
		if (.not.yp2_exist(ivec)) then
		  call tce_filenameindexed(ivec,'yp2',filename)
		  call createfile(filename,yp2(ivec),size_y2)
		  yp2_exist(ivec) = .true.
		  call dratoga(y1(ivec))
		  call dratoga(y2(ivec))
		  call dratoga(y3(ivec))
		  call dratoga(y4(ivec))
!ccx                  call eomccsdtq_y2(d_f1,yp2(ivec),d_t1,d_t2,d_t3,d_t4,
!ccx     1              d_v2,y1(ivec),y2(ivec),y3(ivec),y4(ivec),
!ccx     2              k_f1_offset,k_y2_offset,k_t1_offset,k_t2_offset,
!ccx     3              k_t3_offset,k_t4_offset,k_v2_offset,k_y1_offset,
!ccx     4              k_y2_offset,k_y3_offset,k_y4_offset)
		  call reconcilefile(yp2(ivec),size_y2)
		  call gatodra(y4(ivec))
		  call gatodra(y3(ivec))
		  call gatodra(y2(ivec))
		  call gatodra(y1(ivec))
		  call gatodra(yp2(ivec))
		endif
		if (.not.yp3_exist(ivec)) then
		  call tce_filenameindexed(ivec,'yp3',filename)
		  call createfile(filename,yp3(ivec),size_y3)
		  yp3_exist(ivec) = .true.
		  call dratoga(y1(ivec))
		  call dratoga(y2(ivec))
		  call dratoga(y3(ivec))
		  call dratoga(y4(ivec))
!ccx                    call eomccsdtq_y3(d_f1,yp3(ivec),d_t1,d_t2,d_t3,
!ccx     1              d_v2,y1(ivec),y2(ivec),y3(ivec),y4(ivec),
!ccx     2              k_f1_offset,k_y3_offset,k_t1_offset,k_t2_offset,
!ccx     3              k_t3_offset,k_v2_offset,k_y1_offset,
!ccx     4              k_y2_offset,k_y3_offset,k_y4_offset)
		  call reconcilefile(yp3(ivec),size_y3)
		  call gatodra(y4(ivec))
		  call gatodra(y3(ivec))
		  call gatodra(y2(ivec))
		  call gatodra(y1(ivec))
		  call gatodra(yp3(ivec))
		endif
		if (.not.yp4_exist(ivec)) then
		  call tce_filenameindexed(ivec,'yp4',filename)
		  call createfile(filename,yp4(ivec),size_y4)
		  yp4_exist(ivec) = .true.
		  call dratoga(y1(ivec))
		  call dratoga(y2(ivec))
		  call dratoga(y3(ivec))
		  call dratoga(y4(ivec))
!ccx                  call eomccsdtq_y4(d_f1,yp4(ivec),d_t1,d_t2,
!ccx     1              d_v2,y2(ivec),y3(ivec),y4(ivec),
!ccx     2              k_f1_offset,k_y4_offset,k_t1_offset,k_t2_offset,
!ccx     3              k_v2_offset,k_y2_offset,k_y3_offset,k_y4_offset)
		  call reconcilefile(yp4(ivec),size_y4)
		  call gatodra(y4(ivec))
		  call gatodra(y3(ivec))
		  call gatodra(y2(ivec))
		  call gatodra(y1(ivec))
		  call gatodra(yp4(ivec))
		endif
	      enddo
	      if (.not.ma_push_get(mt_dbl,nytrials,'residual',
     1          l_residual,k_residual))
     2          call errquit('tce_energy: MA problem',101,MA_ERR)
	      call tce_eom_ydiagon(needt1,needt2,needt3,needt4,
     1          size_y1,size_y2,size_y3,size_y4,
     2          k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset,
     3          d_ry1,d_ry2,d_ry3,d_ry4,
     4          dbl_mb(k_omegay),dbl_mb(k_residual))
	      cpu=cpu+util_cpusec()
	      wall=wall+util_wallsec()
	      converged = .true.
	      do ivec = 1,nroots_reduced
		if (nodezero.and.(ivec.ne.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegay+ivec-1),
     3            dbl_mb(k_omegay+ivec-1)*au2ev
		if (nodezero.and.(ivec.eq.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegay+ivec-1),
     3            dbl_mb(k_omegay+ivec-1)*au2ev,cpu,wall
                if (nodezero) call util_flush(LuOut)
		if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1            converged = .false.
	      enddo
	      cpu=-util_cpusec()
	      wall=-util_wallsec()
	      if (.not.ma_pop_stack(l_residual))
     1          call errquit("tce_energy: MA problem",102,MA_ERR)
	      if (converged) then
		if (nodezero) write(LuOut,9240)
		call tce_eom_ytidy
!               call hbar_tidy
		do ivec = 1,nroots_reduced
		  if (nodezero) write(LuOut,9460)
     1              ivec,dbl_mb(k_omegax+ivec-1),
     2              dbl_mb(k_omegax+ivec-1)*au2ev
		  do jvec = 1,nroots_reduced
		  if (dabs(dbl_mb(k_omegax+ivec-1)-
     1                dbl_mb(k_omegay+jvec-1)).lt.thresh) then
		    call tce_filename('d0',filename)
		    call createfile(filename,d_d0,size_d0)
		    call dratoga(xc1(ivec))
		    call dratoga(xc2(ivec))
		    call dratoga(xc3(ivec))
		    call dratoga(xc4(ivec))
		    call dratoga(yc1(jvec))
		    call dratoga(yc2(jvec))
		    call dratoga(yc3(jvec))
		    call dratoga(yc4(jvec))
!ccx                    call eomccsdtq_denominator(d_d0,xc1(ivec),xc2(ivec),
!ccx     1                xc3(ivec),xc4(ivec),yc1(jvec),yc2(jvec),yc3(jvec),
!ccx     2                yc4(jvec),k_d0_offset,k_x1_offset,k_x2_offset,
!ccx     3                k_x3_offset,k_x4_offset,k_y1_offset,k_y2_offset,
!ccx     4                k_y3_offset,k_y4_offset)
		    call reconcilefile(d_d0,size_d0)
		    call gatodra(yc4(jvec))
		    call gatodra(yc3(jvec))
		    call gatodra(yc2(jvec))
		    call gatodra(yc1(jvec))
		    call gatodra(xc4(ivec))
		    call gatodra(xc3(ivec))
		    call gatodra(xc2(ivec))
		    call gatodra(xc1(ivec))
		    call get_block(d_d0,denominator,1,0)
                    call deletefile(d_d0)
		    if (dabs(denominator).gt.thresh) then
		    do axis=1,3
		      irrep_d=sym_abelian_axis(geom,axis)
                      irrep_x = irrep
                      irrep_y = irrep
                      call tce_filename('d0',filename)
                      call createfile(filename,d_d0,size_d0)
		      call tce_filename('x0',filename)
		      call createfile(filename,d_x0,size_x0)
		      call tce_filename('y0',filename)
		      call createfile(filename,d_y0,size_y0)
		      call dratoga(xc1(ivec))
		      call dratoga(xc2(ivec))
		      call dratoga(xc3(ivec))
		      call dratoga(xc4(ivec))
		      call dratoga(yc1(jvec))
		      call dratoga(yc2(jvec))
		      call dratoga(yc3(jvec))
		      call dratoga(yc4(jvec))
!ccx                      call eomccsdtq_density1(d_d1(axis),d_d0,d_t1,d_t2,
!ccx     1                  d_t3,d_t4,d_x0,xc1(ivec),xc2(ivec),xc3(ivec),
!ccx     2                  xc4(ivec),d_y0,yc1(jvec),yc2(jvec),yc3(jvec),
!ccx     3                  yc4(jvec),k_d1_offset(axis),k_d0_offset,
!ccx     4                  k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
!ccx     5                  k_x0_offset,k_x1_offset,k_x2_offset,k_x3_offset,
!ccx     6                  k_x4_offset,k_y0_offset,k_y1_offset,k_y2_offset,
!ccx     7                  k_y3_offset,k_y4_offset)
		      call reconcilefile(d_d0,size_d0)
		      call gatodra(yc4(jvec))
		      call gatodra(yc3(jvec))
		      call gatodra(yc2(jvec))
		      call gatodra(yc1(jvec))
		      call gatodra(xc4(ivec))
		      call gatodra(xc3(ivec))
		      call gatodra(xc2(ivec))
		      call gatodra(xc1(ivec))
		      call get_block(d_d0,numerator,1,0)
                      call deletefile(d_d0)
		      dipole_exc(axis) = dipole_scf(axis)
     1                                 + numerator/denominator
		      transition(axis) = 1.0d0
                      irrep_x = irrep
                      irrep_y = 0
                      call tce_filename('d0',filename)
                      call createfile(filename,d_d0,size_d0)
		      call put_block(d_x0,0.0d0,1,0)
		      call put_block(d_y0,1.0d0,1,0)
		      call dratoga(xc1(ivec))
		      call dratoga(xc2(ivec))
		      call dratoga(xc3(ivec))
		      call dratoga(xc4(ivec))
!ccx                      call eomccsdtq_density1(d_d1(axis),d_d0,d_t1,d_t2,
!ccx     1                  d_t3,d_t4,d_x0,xc1(ivec),xc2(ivec),xc3(ivec),
!ccx     2                  xc4(ivec),d_y0,d_lambda1,d_lambda2,d_lambda3,
!ccx     3                  d_lambda4,k_d1_offset(axis),k_d0_offset,
!ccx     4                  k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
!ccx     5                  k_x0_offset,k_x1_offset,k_x2_offset,k_x3_offset,
!ccx     6                  k_x4_offset,k_y0_offset,k_l1_offset,k_l2_offset,
!ccx     7                  k_l3_offset,k_l4_offset)
		      call reconcilefile(d_d0,size_d0)
		      call gatodra(xc4(ivec))
		      call gatodra(xc3(ivec))
		      call gatodra(xc2(ivec))
		      call gatodra(xc1(ivec))
		      call get_block(d_d0,numerator,1,0)
                      call deletefile(d_d0)
		      transition(axis) = transition(axis) * numerator
                      irrep_x = 0
                      irrep_y = irrep
                      call tce_filename('d0',filename)
                      call createfile(filename,d_d0,size_d0)
		      call put_block(d_x0,1.0d0,1,0)
		      call put_block(d_y0,0.0d0,1,0)
                      call tce_t1_offset(l_r1_offset,k_r1_offset,
     1                  size_r1)
                      call tce_filename('r1',filename)
                      call createfile(filename,d_r1,size_r1)
                      call tce_t2_offset(l_r2_offset,k_r2_offset,
     1                  size_r2)
                      call tce_filename('r2',filename)
                      call createfile(filename,d_r2,size_r2)
                      call tce_t3_offset(l_r3_offset,k_r3_offset,
     1                  size_r3)
                      call tce_filename('r3',filename)
                      call createfile(filename,d_r3,size_r3)
                      call tce_t4_offset(l_r4_offset,k_r4_offset,
     1                  size_r4)
                      call tce_filename('r4',filename)
                      call createfile(filename,d_r4,size_r4)
                      call dratoga(yc1(jvec))
                      call dratoga(yc2(jvec))
                      call dratoga(yc3(jvec))
                      call dratoga(yc4(jvec))
!ccx                      call eomccsdtq_density1(d_d1(axis),d_d0,d_t1,d_t2,
!ccx     1                  d_t3,d_t4,d_x0,d_r1,d_r2,d_r3,
!ccx     2                  d_r4,d_y0,yc1(jvec),yc2(jvec),yc3(jvec),
!ccx     3                  yc4(jvec),k_d1_offset(axis),k_d0_offset,
!ccx     4                  k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
!ccx     5                  k_x0_offset,k_r1_offset,k_r2_offset,k_r3_offset,
!ccx     6                  k_r4_offset,k_y0_offset,k_y1_offset,k_y2_offset,
!ccx     7                  k_y3_offset,k_y4_offset)
                      call reconcilefile(d_d0,size_d0)
                      call gatodra(yc4(jvec))
                      call gatodra(yc3(jvec))
                      call gatodra(yc2(jvec))
                      call gatodra(yc1(jvec))
                      call get_block(d_d0,numerator,1,0)
                      call deletefile(d_d0)
                      transition(axis) = transition(axis)
     1                                 * numerator / denominator
                      call deletefile(d_r4)
                      if (.not.ma_pop_stack(l_r4_offset))
     1                  call errquit("tce_energy: MA problem",100,
     2                  MA_ERR)
                      call deletefile(d_r3)
                      if (.not.ma_pop_stack(l_r3_offset))
     1                  call errquit("tce_energy: MA problem",101,
     2                  MA_ERR)
                      call deletefile(d_r2)
                      if (.not.ma_pop_stack(l_r2_offset))
     1                  call errquit("tce_energy: MA problem",102,
     2                  MA_ERR)
                      call deletefile(d_r1)
                      if (.not.ma_pop_stack(l_r1_offset))
     1                  call errquit("tce_energy: MA problem",103,
     2                  MA_ERR)
                      call deletefile(d_y0)
                      call deletefile(d_x0)
                    enddo
                    if (nodezero) write(LuOut,9470) "EOM-CCSDTQ",
     1                      dsqrt(dabs(transition(1))),
     2                      dsqrt(dabs(transition(2))),
     3                      dsqrt(dabs(transition(3))),
     4                      (transition(1)+transition(2)+transition(3))
     5                      *dbl_mb(k_omegax+ivec-1)*2.0d0/3.0d0
                    dipole_exc_t = dsqrt( dipole_exc(1)*dipole_exc(1) +
     1                                    dipole_exc(2)*dipole_exc(2) +
     2                                    dipole_exc(3)*dipole_exc(3) )
                    if (nodezero) write(LuOut,9430) "EOM-CCSDTQ",
     1                      dipole_exc(1),dipole_exc(1)*debye,
     2                      dipole_exc(2),dipole_exc(2)*debye,
     3                      dipole_exc(3),dipole_exc(3)*debye,
     4                      dipole_exc_t ,dipole_exc_t*debye
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
                    irrep_x = irrep   !! Caution !!
!                    call tce_print_x1(xc1(ivec),k_x1_offset,printtol)
!                    call tce_print_x2(xc2(ivec),k_x2_offset,printtol)
              call tce_print_x1(xc1(ivec),k_x1_offset,printtol,irrep_x)
              call tce_print_x2(xc2(ivec),k_x2_offset,printtol,irrep_x)
                    call gatodra(xc2(ivec))
                    call gatodra(xc1(ivec))
                    if (nodezero) call util_flush(LuOut)
                    endif
                  endif
                  enddo
                enddo
                modelname = "EOM-CCSDTQ"
                goto 400
              endif
            enddo
            call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
            endif
#endif
          else
            call errquit("Not yet implemented",0,CAPMIS_ERR)
          endif
!
!         --------------------
!         Write energy to RTDB
!         --------------------
!
  400     continue
! -------------
! EOM-CCSDT   left
! -------------
          if(ccsdtleft) then
#if !defined(CCSDT)
        CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
            if (left) then
            call tce_eom_yguess(needt1,needt2,needt3,needt4,
     1        size_y1,size_y2,size_y3,size_y4,
     2        k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset)
            if (nytrials.eq.0) goto 4400
            modelname = "EOM-CCSDT left-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
!           call hbar_init
            do iter=1,maxiter
              if (nodezero.and.util_print('eom',print_default))
     1          write(LuOut,9210) iter,nytrials
              do ivec = 1,nytrials
                if (.not.yp1_exist(ivec)) then
                  call tce_filenameindexed(ivec,'yp1',filename)
                  call createfile(filename,yp1(ivec),size_y1)
                  yp1_exist(ivec) = .true.
                  call dratoga(y1(ivec))
                  call dratoga(y2(ivec))
                  call dratoga(y3(ivec))
                  call eomccsdt_y1(d_f1,yp1(ivec),d_t1,d_t2,d_t3,
     1              d_v2,y1(ivec),y2(ivec),y3(ivec),k_f1_offset,
     2              k_y1_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     3              k_v2_offset,k_y1_offset,k_y2_offset,k_y3_offset)
                  call reconcilefile(yp1(ivec),size_y1)
                  call gatodra(y3(ivec))
                  call gatodra(y2(ivec))
                  call gatodra(y1(ivec))
                  call gatodra(yp1(ivec))
                endif
                if (.not.yp2_exist(ivec)) then
                  call tce_filenameindexed(ivec,'yp2',filename)
                  call createfile(filename,yp2(ivec),size_y2)
                  yp2_exist(ivec) = .true.
                  call dratoga(y1(ivec))
                  call dratoga(y2(ivec))
                  call dratoga(y3(ivec))
                  call eomccsdt_y2(d_f1,yp2(ivec),d_t1,d_t2,d_t3,
     1              d_v2,y1(ivec),y2(ivec),y3(ivec),k_f1_offset,
     2              k_y2_offset,k_t1_offset,k_t2_offset,k_t3_offset,
     3              k_v2_offset,k_y1_offset,k_y2_offset,k_y3_offset)
                  call reconcilefile(yp2(ivec),size_y2)
                  call gatodra(y3(ivec))
                  call gatodra(y2(ivec))
                  call gatodra(y1(ivec))
                  call gatodra(yp2(ivec))
                endif
                if (.not.yp3_exist(ivec)) then
                  call tce_filenameindexed(ivec,'yp3',filename)
                  call createfile(filename,yp3(ivec),size_y3)
                  yp3_exist(ivec) = .true.
                  call dratoga(y1(ivec))
                  call dratoga(y2(ivec))
                  call dratoga(y3(ivec))
                  call eomccsdt_y3(d_f1,yp3(ivec),d_t1,d_t2,
     1              d_v2,y1(ivec),y2(ivec),y3(ivec),k_f1_offset,
     2              k_y3_offset,k_t1_offset,k_t2_offset,
     3              k_v2_offset,k_y1_offset,k_y2_offset,k_y3_offset)
                  call reconcilefile(yp3(ivec),size_y3)
                  call gatodra(y3(ivec))
                  call gatodra(y2(ivec))
                  call gatodra(y1(ivec))
                  call gatodra(yp3(ivec))
                endif
              enddo
              if (.not.ma_push_get(mt_dbl,nytrials,'residual',
     1          l_residual,k_residual))
     2          call errquit('tce_energy: MA problem',101,MA_ERR)
              call tce_eom_ydiagon(needt1,needt2,needt3,needt4,
     1          size_y1,size_y2,size_y3,size_y4,
     2          k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset,
     3          d_ry1,d_ry2,d_ry3,d_ry4,
     4          dbl_mb(k_omegay),dbl_mb(k_residual))
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
              converged = .true.
              do ivec = 1,nroots_reduced
                if (nodezero.and.(ivec.ne.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegay+ivec-1),
     3            dbl_mb(k_omegay+ivec-1)*au2ev
                if (nodezero.and.(ivec.eq.nroots_reduced))
     1            write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2            dbl_mb(k_omegay+ivec-1),
     3            dbl_mb(k_omegay+ivec-1)*au2ev,cpu,wall
                if (nodezero) call util_flush(LuOut)
                if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1            converged = .false.
              enddo
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (.not.ma_pop_stack(l_residual))
     1          call errquit("tce_energy: MA problem",102,MA_ERR)
              if (converged) then
                if (nodezero) write(LuOut,9240)
                call tce_eom_ytidy
!               call hbar_tidy
                do ivec = 1,nroots_reduced
                  if (nodezero) write(LuOut,9460)
     1              ivec,dbl_mb(k_omegax+ivec-1),
     2              dbl_mb(k_omegax+ivec-1)*au2ev
                  do jvec = 1,nroots_reduced
                  if (dabs(dbl_mb(k_omegax+ivec-1)-
     1                dbl_mb(k_omegay+jvec-1)).lt.thresh) then
                    call tce_filename('d0',filename)
                    call createfile(filename,d_d0,size_d0)
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
                    call dratoga(xc3(ivec))
                    call dratoga(yc1(jvec))
                    call dratoga(yc2(jvec))
                    call dratoga(yc3(jvec))
                    call eomccsdt_denominator(d_d0,xc1(ivec),xc2(ivec),
     1                xc3(ivec),yc1(jvec),yc2(jvec),yc3(ivec),
     2                k_d0_offset,k_x1_offset,k_x2_offset,k_x3_offset,
     3                k_y1_offset,k_y2_offset,k_y3_offset)
                    call reconcilefile(d_d0,size_d0)
                    call gatodra(yc3(jvec))
                    call gatodra(yc2(jvec))
                    call gatodra(yc1(jvec))
                    call gatodra(xc3(ivec))
                    call gatodra(xc2(ivec))
                    call gatodra(xc1(ivec))
                    call get_block(d_d0,denominator,1,0)
                    call deletefile(d_d0)
                    if (dabs(denominator).gt.thresh) then
                    do axis=1,3
                      irrep_d=sym_abelian_axis(geom,axis)
                      irrep_x = irrep
                      irrep_y = irrep
                      call tce_filename('d0',filename)
                      call createfile(filename,d_d0,size_d0)
                      call tce_filename('x0',filename)
                      call createfile(filename,d_x0,size_x0)
                      call tce_filename('y0',filename)
                      call createfile(filename,d_y0,size_y0)
                      call dratoga(xc1(ivec))
                      call dratoga(xc2(ivec))
                      call dratoga(xc3(ivec))
                      call dratoga(yc1(jvec))
                      call dratoga(yc2(jvec))
                      call dratoga(yc3(jvec))
                      call eomccsdt_density1(d_d1(axis),d_d0,d_t1,d_t2,
     1                  d_t3,d_x0,xc1(ivec),xc2(ivec),xc3(ivec),d_y0,
     2                  yc1(jvec),yc2(jvec),yc3(jvec),k_d1_offset(axis),
     3                  k_d0_offset, k_t1_offset,k_t2_offset,
     4                  k_t3_offset,k_x0_offset,k_x1_offset,
     5                  k_x2_offset,k_x3_offset,k_y0_offset,
     6                  k_y1_offset,k_y2_offset,k_y3_offset)
                      call reconcilefile(d_d0,size_d0)
                      call gatodra(yc3(jvec))
                      call gatodra(yc2(jvec))
                      call gatodra(yc1(jvec))
                      call gatodra(xc3(ivec))
                      call gatodra(xc2(ivec))
                      call gatodra(xc1(ivec))
                      call get_block(d_d0,numerator,1,0)
                      call deletefile(d_d0)
                      dipole_exc(axis) = dipole_scf(axis)
     1                                 + numerator/denominator
                      transition(axis) = 1.0d0
                      irrep_x = irrep
                      irrep_y = 0
                      call tce_filename('d0',filename)
                      call createfile(filename,d_d0,size_d0)
                      call put_block(d_x0,0.0d0,1,0)
                      call put_block(d_y0,1.0d0,1,0)
                      call dratoga(xc1(ivec))
                      call dratoga(xc2(ivec))
                      call dratoga(xc3(ivec))
                      call eomccsdt_density1(d_d1(axis),d_d0,d_t1,d_t2,
     1                  d_t3,d_x0,xc1(ivec),xc2(ivec),xc3(ivec),d_y0,
     2                  d_lambda1,d_lambda2,d_lambda3,k_d1_offset(axis),
     3                  k_d0_offset,k_t1_offset,k_t2_offset,
     4                  k_t3_offset,k_x0_offset,k_x1_offset,
     5                  k_x2_offset,k_x3_offset,k_y0_offset,
     6                  k_l1_offset,k_l2_offset,k_l3_offset)
                      call reconcilefile(d_d0,size_d0)
                      call gatodra(xc3(ivec))
                      call gatodra(xc2(ivec))
                      call gatodra(xc1(ivec))
                      call get_block(d_d0,numerator,1,0)
                      call deletefile(d_d0)
                      transition(axis) = transition(axis) * numerator
                      irrep_x = 0
                      irrep_y = irrep
                      call tce_filename('d0',filename)
                      call createfile(filename,d_d0,size_d0)
                      call put_block(d_x0,1.0d0,1,0)
                      call put_block(d_y0,0.0d0,1,0)
                      call tce_t1_offset(l_r1_offset,k_r1_offset,
     1                  size_r1)
                      call tce_filename('r1',filename)
                      call createfile(filename,d_r1,size_r1)
                      call tce_t2_offset(l_r2_offset,k_r2_offset,
     1                  size_r2)
                      call tce_filename('r2',filename)
                      call createfile(filename,d_r2,size_r2)
                      call tce_t3_offset(l_r3_offset,k_r3_offset,
     1                  size_r3)
                      call tce_filename('r3',filename)
                      call createfile(filename,d_r3,size_r3)
                      call dratoga(yc1(jvec))
                      call dratoga(yc2(jvec))
                      call dratoga(yc3(jvec))
                      call eomccsdt_density1(d_d1(axis),d_d0,d_t1,d_t2,
     1                  d_t3,d_x0,d_r1,d_r2,d_r3,d_y0,
     2                  yc1(jvec),yc2(jvec),yc3(jvec),k_d1_offset(axis),
     3                  k_d0_offset,k_t1_offset,k_t2_offset,
     4                  k_t3_offset,k_x0_offset,k_r1_offset,
     5                  k_r2_offset,k_r3_offset,k_y0_offset,
     6                  k_y1_offset,k_y2_offset,k_y3_offset)
                      call reconcilefile(d_d0,size_d0)
                      call gatodra(yc3(jvec))
                      call gatodra(yc2(jvec))
                      call gatodra(yc1(jvec))
                      call get_block(d_d0,numerator,1,0)
                      call deletefile(d_d0)
                      transition(axis) = transition(axis)
     1                                 * numerator / denominator
                      call deletefile(d_r3)
                      if (.not.ma_pop_stack(l_r3_offset))
     1                  call errquit("tce_energy: MA problem",100,
     2                  MA_ERR)
                      call deletefile(d_r2)
                      if (.not.ma_pop_stack(l_r2_offset))
     1                  call errquit("tce_energy: MA problem",101,
     2                  MA_ERR)
                      call deletefile(d_r1)
                      if (.not.ma_pop_stack(l_r1_offset))
     1                  call errquit("tce_energy: MA problem",102,
     2                  MA_ERR)
                      call deletefile(d_y0)
                      call deletefile(d_x0)
                    enddo
                    if (nodezero) write(LuOut,9470) "EOM-CCSDT",
     1                      dsqrt(dabs(transition(1))),
     2                      dsqrt(dabs(transition(2))),
     3                      dsqrt(dabs(transition(3))),
     4                      (transition(1)+transition(2)+transition(3))
     5                      *dbl_mb(k_omegax+ivec-1)*2.0d0/3.0d0
                    dipole_exc_t = dsqrt( dipole_exc(1)*dipole_exc(1) +
     1                                    dipole_exc(2)*dipole_exc(2) +
     2                                    dipole_exc(3)*dipole_exc(3) )
                    if (nodezero) write(LuOut,9430) "EOM-CCSDT",
     1                      dipole_exc(1),dipole_exc(1)*debye,
     2                      dipole_exc(2),dipole_exc(2)*debye,
     3                      dipole_exc(3),dipole_exc(3)*debye,
     4                      dipole_exc_t ,dipole_exc_t*debye
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
                    irrep_x = irrep   !! Caution !!
!                    call tce_print_x1(xc1(ivec),k_x1_offset,printtol)
!                    call tce_print_x2(xc2(ivec),k_x2_offset,printtol)
              call tce_print_x1(xc1(ivec),k_x1_offset,printtol,irrep_x)
              call tce_print_x2(xc2(ivec),k_x2_offset,printtol,irrep_x)
                    call gatodra(xc2(ivec))
                    call gatodra(xc1(ivec))
                    if (nodezero) call util_flush(LuOut)
                    endif
                  endif
                  enddo
                enddo
                modelname = "EOM-CCSDT"
                goto 4400
              endif
            enddo
            call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
            endif
#endif
         endif
!     done with 2nd ccsdtleft
 4400    ccsdtleft=.false.
          call tce_eom_tidy
          if ((targetsym.eq.irrepname)
     1      .and.(target.le.nroots_reduced)) then
            if (nodezero.and.util_print('excited state',print_low)) then
              write(LuOut,*)
              write(LuOut,9300) target
              write(LuOut,9340) targetsym
              write(LuOut,*)
! === jaguar ===
              if(.not.(read_ta.or.read_xa.or.read_in3)) then 
! ==============
               write(LuOut,9310) modelname(1:inp_strlen(modelname)),
     1           ref+corr
               write(LuOut,9320) modelname(1:inp_strlen(modelname)),
     1           dbl_mb(k_omegax+target-1)
               write(LuOut,9330) modelname(1:inp_strlen(modelname)),
     1           ref+corr+dbl_mb(k_omegax+target-1)
              end if !jaguar
            endiF
            if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,
     1        ref+corr+dbl_mb(k_omegax+target-1)))
     2        call errquit('tce_energy: RTDB problem',0,RTDB_ERR)
            tce_energy=.true.
          elseif ((targetsym.eq.'none').and.(irrep_x.eq.0)
     1      .and.(target.le.nroots_reduced)) then
            if (nodezero.and.util_print('excited state',print_low)) then
!cx              write(LuOut,*)
!cx              write(LuOut,9310) modelname(1:inp_strlen(modelname)),
!cx     1          ref+corr
!cx              write(LuOut,9320) modelname(1:inp_strlen(modelname)),
!cx     1          dbl_mb(k_omegax+target-1)
!cx              write(LuOut,9330) modelname(1:inp_strlen(modelname)),
!cx     1          ref+corr+dbl_mb(k_omegax+target-1)
            endif
            if (.not.rtdb_put(rtdb,'tce:omega',mt_dbl,1,
     1        dbl_mb(k_omegax+target-1)*au2ev))
     2        call errquit('tce_energy: RTDB problem',0,RTDB_ERR)
            if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,
     1        ref+corr+dbl_mb(k_omegax+target-1)))
     2        call errquit('tce_energy: RTDB problem',0,RTDB_ERR)
            tce_energy=.true.
          endif
! link to QM/MM (no symmetry OR a-symmetry) -------
          if(irrep.eq.0) then
            if (.not.rtdb_put(rtdb,'tce:eneex',mt_dbl,nroots_reduced,
     1        dbl_mb(k_omegax)))
     2        call errquit('tce_energy: RTDB problem',0,RTDB_ERR)
          end if
! -------------------------------------------------
!
!         ------------------------------------
!         Terminate excited state calculations
!         ------------------------------------
!
          if (needt4) then
            if (left) then
              call deletefile(d_ry4)
              if (.not.ma_pop_stack(l_y4_offset))
     1          call errquit("tce_energy: MA problem",34,MA_ERR)
            endif
            call deletefile(d_rx4)
            if (.not.ma_pop_stack(l_x4_offset))
     1        call errquit("tce_energy: MA problem",34,MA_ERR)
          endif
          if (needt3) then
            if (left) then
              call deletefile(d_ry3)
              if (.not.ma_pop_stack(l_y3_offset))
     1          call errquit("tce_energy: MA problem",35,MA_ERR)
            endif
            call deletefile(d_rx3)
            if (.not.ma_pop_stack(l_x3_offset))
     1        call errquit("tce_energy: MA problem",35,MA_ERR)
          endif
!cc
          if (needt3a) then
!ccx           if (left) then
!ccx              call deletefile(d_ry3)
!ccx              if (.not.ma_pop_stack(l_y3_offset))
!ccx     1          call errquit("tce_energy: MA problem",35,MA_ERR)
!ccx            endif
            call deletefile(d_rx3)
            if (.not.ma_pop_stack(l_x3_offset))
     1        call errquit("tce_energy: MA problem",35,MA_ERR)
          endif
!cc
          if (needt2) then
            if (left) then
              call deletefile(d_ry2)
              if (.not.ma_pop_stack(l_y2_offset))
     1          call errquit("tce_energy: MA problem",36,MA_ERR)
            endif
            call deletefile(d_rx2)
            if (.not.ma_pop_stack(l_x2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
          endif
          if (needt1) then
            if (left) then
              call deletefile(d_ry1)
              if (.not.ma_pop_stack(l_y1_offset))
     1          call errquit("tce_energy: MA problem",37,MA_ERR)
            endif
            call deletefile(d_rx1)
            if (.not.ma_pop_stack(l_x1_offset))
     1        call errquit("tce_energy: MA problem",37,MA_ERR)
          endif
          if (left) then
            if (.not.ma_pop_stack(l_y0_offset))
     1        call errquit("tce_energy: MA problem",38,MA_ERR)
            if (.not.ma_pop_stack(l_x0_offset))
     1        call errquit("tce_energy: MA problem",38,MA_ERR)
            if (.not.ma_pop_stack(l_d0_offset))
     1        call errquit("tce_energy: MA problem",38,MA_ERR)
            if (.not.ma_pop_stack(l_omegay))
     1        call errquit("tce_energy: MA problem",1003,MA_ERR)
          endif
          if (.not.ma_pop_stack(l_omegax))
     1      call errquit("tce_energy: MA problem",1002,MA_ERR)
! new EOMCC solver
       if(.not.((model.ne.'ccsd').and.(eom_solver.eq.2))) then
! ccsd_act/eomccsd_act
       if(.not.(model.eq.'ccsd_act')) then
       if (.not.ma_pop_stack(l_hbar))
     1  call errquit('tce_eom_xdiagon: MA problem',12,MA_ERR)
      end if
      end if
          endif
        enddo
      endif
!
!     ===========================
!     End of all TCE calculations
!     ===========================
!
!->d3p975
      do i = 1, ipol
        if (.not.ga_destroy(g_movecs(i)))
     1    call errquit('tce_energy: GA problem',0,GA_ERR)
      enddo
!<-d3p975
!
      if (io_v2.eq.IOSF) then
        call deletefile_sf(d_v2)
      else
! 
! === jaguar ===
        if(.not.read_in3) then
         call deletefile(d_v2)
        end if
! ==============
      endif
c
      if(intorb.and.(model.eq.'ccsdt')) then 
       call deletefile(d_v2aux)
      endif
c
      if (left) then
        if (multipole.gt.2) then
          do axis=10,1,-1
            call deletefile(d_o1(axis))
          enddo
        endif
        if (multipole.gt.1) then
          do axis=6,1,-1
            call deletefile(d_q1(axis))
          enddo
        endif
        do axis=3,1,-1
          call deletefile(d_d1(axis))
        enddo
      endif
! === jaguar ===
      if(.not.read_in3) then
      call deletefile(d_f1)
      end if
! ==============
      if (needt4) then
        if (leftresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_yr4(axis))
            call deletefile(d_yr4(axis+3))
            if (.not.ma_pop_stack(l_yr4_offset(axis)))
     1        call errquit("tce_energy: MA problem",86,MA_ERR)
          endif
          enddo ! axis
        endif ! leftresp
        if (lineresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_tr4(axis))
            call deletefile(d_tr4(axis+3))
            if (.not.ma_pop_stack(l_tr4_offset(axis)))
     1        call errquit("tce_energy: MA problem",89,MA_ERR)
          endif
          enddo ! axis
        endif ! lineresp
        if (left) then
          call deletefile(d_lambda4)
          if (.not.ma_pop_stack(l_l4_offset))
     1      call errquit("tce_energy: MA problem",26,MA_ERR)
        endif
        call deletefile(d_t4)
        if (.not.ma_pop_stack(l_t4_offset))
     1    call errquit("tce_energy: MA problem",26,MA_ERR)
      endif
      if (needt3) then
        if (leftresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_yr3(axis))
            call deletefile(d_yr3(axis+3))
            if (.not.ma_pop_stack(l_yr3_offset(axis)))
     1        call errquit("tce_energy: MA problem",86,MA_ERR)
          endif
          enddo ! axis
        endif ! leftresp
        if (lineresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_tr3(axis))
            call deletefile(d_tr3(axis+3))
            if (.not.ma_pop_stack(l_tr3_offset(axis)))
     1        call errquit("tce_energy: MA problem",88,MA_ERR)
          endif
          enddo ! axis
        endif ! lineresp
        if (left) then
          call deletefile(d_lambda3)
          if (.not.ma_pop_stack(l_l3_offset))
     1      call errquit("tce_energy: MA problem",25,MA_ERR)
        endif
        call deletefile(d_t3)
        if (.not.ma_pop_stack(l_t3_offset))
     1    call errquit("tce_energy: MA problem",25,MA_ERR)
      endif
      if (needt3a) then
!ccx        if (left) then
!ccx          call deletefile(d_lambda3)
!ccx          if (.not.ma_pop_stack(l_l3_offset))
!ccx     1      call errquit("tce_energy: MA problem",25,MA_ERR)
!ccx        endif
        call deletefile(d_t3)
        if (.not.ma_pop_stack(l_t3_offset))
     1    call errquit("tce_energy: MA problem",25,MA_ERR)
      endif
      if (needt2) then
        if (leftresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_yr2(axis))
            call deletefile(d_yr2(axis+3))
            if (.not.ma_pop_stack(l_yr2_offset(axis)))
     1        call errquit("tce_energy: MA problem",86,MA_ERR)
          endif
          enddo ! axis
        endif ! leftresp
        if (lineresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_tr2(axis))
            call deletefile(d_tr2(axis+3))
            if (.not.ma_pop_stack(l_tr2_offset(axis)))
     1        call errquit("tce_energy: MA problem",87,MA_ERR)
          endif
          enddo ! axis
        endif ! lineresp
        if (left) then
          call deletefile(d_lambda2)
          if (.not.ma_pop_stack(l_l2_offset))
     1      call errquit("tce_energy: MA problem",0,MA_ERR)
        endif
        call deletefile(d_t2)
        if (.not.ma_pop_stack(l_t2_offset))
     1    call errquit("tce_energy: MA problem",0,MA_ERR)
      endif
      if (needt1) then
        if (leftresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_yr1(axis))
            call deletefile(d_yr1(axis+3))
            if (.not.ma_pop_stack(l_yr1_offset(axis)))
     1        call errquit("tce_energy: MA problem",86,MA_ERR)
          endif
          enddo ! axis
        endif ! leftresp
        if (lineresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_tr1(axis))
            call deletefile(d_tr1(axis+3))
            if (.not.ma_pop_stack(l_tr1_offset(axis)))
     1        call errquit("tce_energy: MA problem",86,MA_ERR)
          endif
          enddo ! axis
        endif ! lineresp
        if (left) then
          call deletefile(d_lambda1)
          if (.not.ma_pop_stack(l_l1_offset))
     1      call errquit("tce_energy: MA problem",1,MA_ERR)
        endif
        call deletefile(d_t1)
        if (.not.ma_pop_stack(l_t1_offset))
     1    call errquit("tce_energy: MA problem",1,MA_ERR)
      endif
      if(intorb) then
       if(intorb.and.(model.eq.'ccsdt')) then 
        if(.not.ma_pop_stack(l_v2aux_offset))
     1   call errquit("tce_energy: MA problem",2,MA_ERR)
       endif
! === jaguar ===
       if(.not.read_in3) then ! --- read_in3
! ==============
       if (.not.ma_pop_stack(l_v2spin_offset))
     1  call errquit("tce_energy: MA problem",59,MA_ERR)
!--       if (.not.ma_pop_stack(l_o2ot_offset))
!--     1  call errquit("tce_energy: MA problem",58,MA_ERR)
!--       if (.not.ma_pop_stack(l_o2o_offset))
!--     1  call errquit("tce_energy: MA problem",57,MA_ERR)
       if (.not.ma_pop_stack(l_v2_alpha_offset))
     1  call errquit("tce_energy: MA problem",56,MA_ERR)
       end if ! read_in3
      end if
      if(.not.intorb) then
      if (.not.ma_pop_stack(l_v2_offset))
     1  call errquit("tce_energy: MA problem",2,MA_ERR)
      end if
! === jaguar ===
      if(.not.read_in3) then ! --- read_in3 ---
      if (.not.ma_pop_stack(l_f1_offset))
     1  call errquit("tce_energy: MA problem",3,MA_ERR)
      end if
! ==============
      if (left) then
        if (multipole.gt.2) then
          do axis=10,1,-1
            if (.not.ma_pop_stack(l_o1_offset(axis)))
     1        call errquit("tce_energy: MA problem",3030+axis,MA_ERR)
          enddo
        endif
        if (multipole.gt.1) then
          do axis=6,1,-1
            if (.not.ma_pop_stack(l_q1_offset(axis)))
     1        call errquit("tce_energy: MA problem",3020+axis,MA_ERR)
          enddo
        endif
        do axis=3,1,-1
          if (.not.ma_pop_stack(l_d1_offset(axis)))
     1      call errquit("tce_energy: MA problem",3010+axis,MA_ERR)
        enddo
      endif
      if(intorb) then
! === jaguar ===
       if(.not.read_in3) then ! ---
! ==============
       if (.not.ma_pop_stack(l_offset_alpha))
     1  call errquit("tce_energy: MA problem",55,MA_ERR)
       if (.not.ma_pop_stack(l_range_alpha))
     1  call errquit("tce_energy: MA problem",54,MA_ERR)
       if (.not.ma_pop_stack(l_sym_alpha))
     1  call errquit("tce_energy: MA problem",53,MA_ERR)
       if (.not.ma_pop_stack(l_spin_alpha))
     1  call errquit("tce_energy: MA problem",52,MA_ERR)
       if (.not.ma_pop_stack(l_b2am))
     1  call errquit("tce_energy: MA problem",51,MA_ERR)
       end if ! --- read_in3
      end if
! ---
! check if read_in3 is needed here
! ---
      if(model.eq."ccsd_act") THEN
      if (.not.ma_pop_stack(l_active_o))
     1  call errquit("tce_energy: MA problem",35,MA_ERR)
      end if
      if(activecalc.or.(.not.intorb).or.(model.eq."ccsd_act")) THEN
      if (.not.ma_pop_stack(l_active))
     1  call errquit("tce_energy: MA problem",34,MA_ERR)
      end if
! === jaguar ===
      if(.not.read_in3) then !---
! ==============
      if (.not.ma_pop_stack(l_alpha))
     1  call errquit("tce_energy: MA problem",4,MA_ERR)
      if (.not.ma_pop_stack(l_offset))
     1  call errquit("tce_energy: MA problem",5,MA_ERR)
      if (.not.ma_pop_stack(l_range))
     1  call errquit("tce_energy: MA problem",6,MA_ERR)
      if (.not.ma_pop_stack(l_sym))
     1  call errquit("tce_energy: MA problem",7,MA_ERR)
      if (.not.ma_pop_stack(l_spin))
     1  call errquit("tce_energy: MA problem",8,MA_ERR)
      if(.not.ma_pop_stack(l_mo_index))
     1  call errquit("tce_energy: ma problem",85,ma_err)
      if (.not.ma_pop_stack(l_evl_sorted))
     1  call errquit("tce_energy: MA problem",9,MA_ERR)
      if (.not.ma_pop_stack(l_irs_sorted))
     1  call errquit("tce_energy: MA problem",10,MA_ERR)
      if (.not.ma_pop_stack(l_spin_sorted))
     1  call errquit("tce_energy: MA problem",11,MA_ERR)
      if (.not.ma_pop_stack(l_movecs_sorted))
     1  call errquit("tce_energy: MA problem",12,MA_ERR)
! ====
      end if
! ====
!
!     ===============
!     Destroy a mutex
!     ===============
!
      if(need_mutex) then
         if (.not.ga_destroy_mutexes())
     1        call errquit('tce_energy: GA problem',1,GA_ERR)
      endif
!     
!     =========
!     Terminate
!     =========
!
      if(.not.read_in3) then
       call tce_tidy(rtdb)
       call util_print_pop
      end if
      return
!
!     ======
!     Format
!     ======
!
 9000 format(1x,A,' file size   = ',i16)
 9010 format(1x,A,' file name   = ',A)
 9090 format(1x,A,' file handle = ',i10)
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9480 format(1x,'Cpu & wall time / sec for ',A,2f15.1)
 9110 format(1x,'MBPT(0) energy / hartree             = ',f25.15)
 9030 format(/,1x,'MBPT(2) correlation energy / hartree = ',f25.15)
 9040 format(1x,'MBPT(2) total energy / hartree       = ',f25.15)
 9050 format(/,1x,A,' iterations',/,
     1  1x,'--------------------------------------------------------',/
     2  1x,'Iter          Residuum       Correlation     Cpu    Wall',/
     3  1x,'--------------------------------------------------------')
 9060 format(
     1  1x,'--------------------------------------------------------',/
     2  1x,'Iterations converged')
 9070 format(1x,A,' correlation energy / hartree = ',f25.15)
 9075 format(1x,A,'  correction energy / hartree = ',f25.15)
 9080 format(1x,A,' total energy / hartree       = ',f25.15)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9083 format(1x,'T1 diagnostic        = ',f25.15)
 9085 format(1x,A,' excitation energy (eV)       = ',f12.5)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9121 format(/,1x,A)
 9122 format(1x,A,i4)
 9123 format(1x,'4-index algorithm nr.',i4,1x,'is used')
 9124 format(1x,'imaxsize = ',2x,i6)
 9125 format(1x,'imaxsize ichop = ',2x,2i6)
! 9126 format(1x,'imaxsize = ',2x,i6)
 9250 format(1x,'Ground-state symmetry is ',A4)
 9251 format(1x,'Dim. of EOMCC iter. space ',2x,i6)
 9200 format(1x,'=========================================',/,
     1       1x,'Excited-state calculation ( ',A4,'symmetry)',/,
     2       1x,'=========================================')
 9210 format(/,1x,'Iteration ',i3,' using ',i4,' trial vectors')
 9220 format(/,1x,A,' iterations',/,1x,
     1'--------------------------------------------------------------'
     2,/,1x,
     3'     Residuum       Omega / hartree  Omega / eV    Cpu    Wall'
     4,/,1x,
     5'--------------------------------------------------------------')
 9211 format(/,1x,'Iteration ',i5)
 9230 format(1x,f17.13,f18.13,f11.5,2f8.1)
 9233 format(/,1x,'Iterations for root',1x,i3,1x,'failed to converge')
 9240 format(1x,
     1'--------------------------------------------------------------'
     2,/,1x,'Iterations converged')
 9300 format(' Target root     =',i3)
 9340 format(' Target symmetry =',1x,a4)
 9310 format(1x,A,' ground state energy / hartree  =',f25.15)
 9320 format(1x,A,' excitation energy / hartree    =',f25.15)
 9330 format(1x,A,' excited state energy / hartree =',f25.15)
 9400 format(/,1x,A,' iterations',/,
     1  1x,'---------------------------------------------',/
     2  1x,'Iter          Residuum            Cpu    Wall',/
     3  1x,'---------------------------------------------')
 9410 format(
     1  1x,'---------------------------------------------',/
     2  1x,'Iterations converged')
 9420 format(1x,i4,f25.13,2f8.1)
 9430 format(/,1x,A,' dipole moments / hartree & Debye',/,
     1  1x,'------------------------------------',/
     2  1x,'  X   ',2f15.7,/
     3  1x,'  Y   ',2f15.7,/
     4  1x,'  Z   ',2f15.7,/
     5  1x,'Total ',2f15.7,/
     6  1x,'------------------------------------')
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
 9432 format(/,1x,A,' ',A,' moments / hartree & Buckingham',/,
     1  1x,'--------------------------------',/
     2  1x,'XX ',2f15.7,/
     3  1x,'XY ',2f15.7,/
     4  1x,'XZ ',2f15.7,/
     5  1x,'YY ',2f15.7,/
     6  1x,'YZ ',2f15.7,/
     7  1x,'ZZ ',2f15.7,/
     8  1x,'--------------------------------')
 9433 format(/,1x,A,' ',A,' moments / hartree & Octos',/,
     1  1x,'--------------------------------',/
     2  1x,'XXX ',2f15.7,/
     3  1x,'XXY ',2f15.7,/
     4  1x,'XXZ ',2f15.7,/
     5  1x,'XYY ',2f15.7,/
     6  1x,'XYZ ',2f15.7,/
     7  1x,'XZZ ',2f15.7,/
     8  1x,'YYY ',2f15.7,/
     9  1x,'YYZ ',2f15.7,/
     1  1x,'YZZ ',2f15.7,/
     2  1x,'ZZZ ',2f15.7,/
     3  1x,'--------------------------------')
 9434 format(/,1x,A,' polarizability / au ',/
     1  1x,'Frequency  = ',f15.7,' / au',/
     1  1x,'Wavelength = ',f15.7,' / nm',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     2  1x,'                    atomic units (bohr^3)       '
     2    ,'       |                   angstroms^3           ',/
     2  1x,'                 X              Y              Z',
     2  1x,'      |             X              Y              Z',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     4  1x,'X      ',3f15.7,3x,'|',3f15.7,/
     5  1x,'Y      ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Z      ',3f15.7,3x,'|',3f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     6  1x,'Eigs = ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Isotropic   = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     6  1x,'Anisotropic = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------')
 9435 format(/,1x,A,' C6 coefficients ',/
     1  1x,'--------------------------------',/
     2  1x,'C6(XX)  ',f15.7,/
     3  1x,'C6(YY)  ',f15.7,/
     4  1x,'C6(ZZ)  ',f15.7,/
     5  1x,'C6(XY)  ',f15.7,/
     6  1x,'C6(XZ)  ',f15.7,/
     7  1x,'C6(YZ)  ',f15.7,/
     8  1x,'C6(AVG) ',f15.7,/
     9  1x,'C6(ANI) ',f15.7,/
     1  1x,'--------------------------------')
 9911 format(/,1x,'Warning: you have not solved ',
     1            'the response equations for all axes.  ',
     2            'Please analyze the results carefully as ',
     3            'the average and anisotropic polarizabilities ',
     4            'are surely wrong.',/)
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
 9450 format(1x,'(T) & [T] code has been provided by ',
     1  'A.A.Auer (Waterloo)')
 9460 format(/,1x,'Excited state root',i3,/
     1       1x,'Excitation energy / hartree =',f25.15,/
     2       1x,'                  / eV      =',f25.15)
 9470 format(/,1x,A,' transition moments / hartree',/
     1  1x,'--------------------------------------------',/
     2  1x,'X',f13.7,' Y',f13.7,' Z',f13.7,/
     3  1x,'Oscillator Strength            ',f13.7,/
     4  1x,'--------------------------------------------')
 9501 format(/,1x,A,' hyperpolarizability / au ')
! Use 9431 for frequency
 9502 format(1x,'-----------------------------------------------')
 9503 format(1x,'beta(',a1,',',a1,',',a1,') = ',f18.9,a)
 9504 format(1x,'Static Hyperpolarizability')
 9505 format(1x,'Second Harmonic Generation (SHG)')
 9506 format(1x,'Optical Rectification (OR)')
 9507 format(1x,'Electro-Optic Pockels Effect (EOPE)')
 9508 format(1x,'beta(',f1.5,',',f1.5,',',f1.5,')')
 9510 format(/,1x,'EOMCCSD SOLVER TYPE ',i2)
 9511 format(1x,'NUMBER OF ROOTS ',i2)
 9512 format(1x,'NO R0 RESTARTED')
 9515 format(3x,'TCE RESTART OPTIONS')
 9520 format(3x,'     L-threshold : ',d9.2)
 9521 format(3x,'     EOM-threshold : ',d9.2)
!QuantumComputing
 4881  format(4i5,2x,f18.10)
 4882  format(2i5,2x,f18.10)
 6887 format(1x,' EHF(A)             = ',f18.8)
 6888 format(1x,' EHF(total)         = ',f18.8)
 6889 format(1x,' Shift (HFtot-HFA)  = ',f18.8)
!QuantumComputing
      end
      subroutine tce_hbarinit(hbar,hbard)
      implicit none
      integer hbard
      double precision hbar(hbard,hbard)
!
      integer i,j
!
      do i=1,hbard
         do j=1,hbard
            hbar(j,i)=0d0
        enddo
      enddo
      do i=1,hbard
         hbar(i,i)=1.0d+8
      enddo
      return
      end
