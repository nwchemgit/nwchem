      logical function tce_energy(rtdb,excitedfragment)
c
c $Id: tce_energy.F,v 1.164 2009-03-16 03:49:57 jhammond Exp $
c 
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
c#include "paw_params.fh"
c
c     CI, CC, & MBPT
c
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
c --- relevant if intorb = .true. ----
      integer lax_v2_alpha_offset     ! Offset for 2e integrals if intorb=.true.
      integer kax_v2_alpha_offset     ! Offset for 2e integrals if intorb=.true.
      integer lax_v2_alpha_offset_ga  ! Offset for 2e integrals if intorb=.true.
      integer kax_v2_alpha_offset_ga  ! Offset for 2e integrals if intorb=.true.
      integer lax_v2_alpha_offset_sf  ! Offset for 2e integrals if intorb=.true.
      integer kax_v2_alpha_offset_sf  ! Offset for 2e integrals if intorb=.true.
c      integer l_o2o_offset       ! Offset for offset to offset table
c      integer k_o2o_offset       ! Offset for offset to offset table
c      integer l_o2ot_offset      ! orbitals for offset to offset tables
c      integer k_o2ot_offset      ! orbitals for offset to offset tables
      integer lx_o2o_offset,kx_o2o_offset
      integer lx_o2ot_offset,kx_o2ot_offset
      integer lx_v2spin_offset,kx_v2spin_offset
      integer length
c      integer d_v2orb
c ------------------------------------
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
c FOR LRH-SD
      integer d_j1,d_j2
      integer l_j1_offset,l_j2_offset
      integer k_j1_offset,k_j2_offset
      integer size_j1,size_j2
c FOR CR-EOMCCSD(T) ---------------------
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
c  NEW SOLVER ---------------------------
      integer l_hbar,k_hbar,ihbar
      integer imicro,xdiis,iroot
      integer l_c_eff,k_c_eff
      integer l_h_eff,k_h_eff
      integer ipos,ipos1,ipos2
c - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer l_x1_local,k_x1_local
c ---------------------------------------
c - intermediates reuse --
      integer l_aux_off,k_aux_off,size_aux
      integer ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5
      integer rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,rec_d_6
      integer rec_d_7,rec_d_8,rec_d_9
c ---------------------------------------
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
c --- debug ---
c      integer dup1,dup2,dup3,dup4,dup5,dup6
c -------------
c
c     CC Lambda equation
c
#if !defined(NEWDIPOLE)
      integer g_aod1(3)         ! GA handle for AO dipole matrices
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
      parameter (octo=0.711688d0)
      parameter (buck=1.344911d0)
      parameter (debye=2.54176568d0)
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
      double precision dipole_scf(3)     ! Dipole moment
      double precision dipole_scf_t      ! Dipole moment
      double precision dipole_cor(3)     ! Dipole moment
      double precision dipole_cor_t      ! Dipole moment
      double precision dipole_exc(3)     ! Dipole moment
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
c
c     CC Response
c
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
c
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
c     RESPONSE AMPLITUDES
      integer d_tr1(9)             ! SF handle for tr1 amplitudes
      integer l_tr1_offset(3)      ! Offset for tr1 file
      integer k_tr1_offset(3)      ! Offset for tr1 file
      integer size_tr1(3)          ! File size in doubles
      integer d_tr2(9)             ! SF handle for tr2 amplitudes
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
c     RESPONSE RESIDUALS
      integer d_rr1(3)             ! SF handle for rr1 amplitudes
      integer d_rr2(3)             ! SF handle for rr2 amplitudes
      integer d_rr3(3)             ! SF handle for rr3 amplitudes
      integer d_rr4(3)             ! SF handle for rr4 amplitudes
c     =====================================================
c               CC-T(1) Imaginary Auxilary Amplitudes
c               Use: ir = A*x, rr = A*ir (rr=A*A*x)
c     =====================================================
      integer d_denom1,d_denom2
      integer d_ir1(3)             ! SF handle for ir1 amplitudes
      integer d_ir2(3)             ! SF handle for ir2 amplitudes
      integer d_ir3(3)             ! SF handle for ir3 amplitudes
      integer d_ir4(3)             ! SF handle for ir4 amplitudes
      character*4 ir1filename(3) ! File name stub
      data ir1filename/'ir1x','ir1y','ir1z'/
      character*4 ir2filename(3) ! File name stub
      data ir2filename/'ir2x','ir2y','ir1z'/
      character*4 ir3filename(3) ! File name stub
      data ir3filename/'ir3x','ir3y','ir3z'/
      character*4 ir4filename(3) ! File name stub
      data ir4filename/'ir4x','ir4y','ir4z'/
c     LAMBDA RESPONSE AMPLITUDES
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
c     LAMBDA RESPONSE RESIDUALS
      integer d_yrr1(3)             ! SF handle for yrr1 amplitudes
      integer d_yrr2(3)             ! SF handle for yrr2 amplitudes
      integer d_yrr3(3)             ! SF handle for yrr3 amplitudes
      integer d_yrr4(3)             ! SF handle for yrr4 amplitudes
c     OTHER RESPONSE STUFF
c      double precision rr1        ! Residual tr1
c      double precision rr2        ! Residual tr2
c      double precision rr3        ! Residual tr3
c      double precision rr4        ! Residual tr4
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
c      double precision gamm(3,3,3,3)! Dipole second hyperpolarizability tensor (gamma not safe variable name)
c      double precision gamm1        ! Dipole second hyperpolarizability term 1
c      double precision gamm2        ! Dipole second hyperpolarizability term 2
c      double precision gamm3        ! Dipole second hyperpolarizability term 3
c      double precision gamm4        ! Dipole second hyperpolarizability term 4
      double precision rr1,rr2,rr3,rr4
      double precision t1diagn ! T1 diagnostic
      double precision alphacopy(3,3)   ! Dipole polarizability tensor copy
      double precision alphaiso         ! Isotropic dipole polarizability
      double precision alphaani         ! Anisotropic dipole polarizability
      double precision alphaevr(3)      ! Dipole polarizability tensor eigenvalues (real)
      double precision alphaevi(3)      ! Dipole polarizability tensor eigenvalues (imag)
      double precision aeigvec(3,3)     ! dummy dgeev eigenvectors
c      integer iter             ! Iteration counter
c
c     Excited-state CI, CC, & MBPT
c
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
      logical converged
      double precision pt3_1   ! [T] energy
      double precision pt3_2   ! (T)-[T] energy
c LR-CCSD(T) method
      double precision pt3_3,pt3_4,pt3_5,pt3_6
c --------------------------------------------
c LR-CCSD(TQ) method
      double precision pt4_lr
c --------------------------------------------
c
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
      logical do_ccd, do_lccd, do_ccsd, do_lccsd, do_qcisd
      logical do_pt, do_lam_pt, do_cr_pt, do_lcr_pt, do_2t_pt
c
c     Fragment extension
c
      logical excitedfragment
c
      logical ccsdtleft
      logical tce_skip_scf
      nodezero=(ga_nodeid().eq.0)
      ccsdtleft=.false.
c
c     =========================
c     Ground-state HF/DFT first
c     =========================
c
      tce_energy = .false.
      if (.not.rtdb_get(rtdb,'tce:reference',mt_int,1,reference)) then
        reference = 1
        if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1    call errquit('tce_energy: failed to write reference',0,
     2    RTDB_ERR)
      endif
c
c     TCE restart precludes SCF calculation
c
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
c
c     ===================
c     Print utility start
c     ===================
c
      call util_print_push
      call util_print_rtdb_load(rtdb,'tce')
c
c     ==========
c     Initialize
c     ==========
c
      call tce_init(rtdb)
      if (.not.rtdb_get(rtdb,'tce:densmat',mt_log,1,idens))
     1   call errquit('tce_energy:rtdb_get densmat failed',0, RTDB_ERR)
      if (nodezero) call util_flush(LuOut)
c      if (lineresp) left=.true.
c      if (imagresp) lineresp=.false.
c
c      if (nodezero) then
c        write(LuOut,*) "lineresp = ",lineresp
c        write(LuOut,*) "leftresp = ",leftresp
c        write(LuOut,*) "quadresp = ",quadresp
c        write(LuOut,*) "anumfreq = ",anumfreq
c        do i=1,anumfreq
c          write(LuOut,*) "afreq(",i,") = ",afreq(i)
c        enddo
c        write(LuOut,*) "respaxis(1) = ",respaxis(1)
c        write(LuOut,*) "respaxis(2) = ",respaxis(2)
c        write(LuOut,*) "respaxis(3) = ",respaxis(3)
c        print*,"respaxis ",respaxis(1),respaxis(2),respaxis(3)
c        call util_flush(LuOut)
c      endif
c
c     =================
c     Totally symmetric  RESPONSE WILL CHANGE THESE - BE CAREFUL
c     =================
c
      irrep_x = 0
      irrep_y = 0
      irrep_c = 0
      irrep_g = 0
c
c     ==============
c     Create a mutex
c     ==============
c
      if (.not.ga_create_mutexes(1))
     1  call errquit('tce_energy: GA problem',0,GA_ERR)
c
c     =========================
c     Reorder and tile orbitals
c     =========================
c
      call tce_tile(rtdb)
      call sf_test
c
c     ================================
c     Parallel integral transformation
c     ================================
c
      if (left) then
#if !defined(NEWDIPOLE)
        cpu = - util_cpusec()
        wall = - util_wallsec()
        call tce_nud1(dipole_scf)
        call tce_aod1(rtdb,g_aod1,dipole_scf)
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
          write(LuOut,*)
          write(LuOut,9020) cpu, wall
          call util_flush(LuOut)
        endif
        do axis = 1, 3
          irrep_d=sym_abelian_axis(geom,axis)
          call sym_irrepname(geom,irrep_d+1,irrepname)
          if (nodezero.and.util_print('mod1',print_default)) then
            write(LuOut,*)
            write(LuOut,9440) axisname(axis),irrepname
          endif
          call tce_mod1_offset(l_d1_offset(axis),k_d1_offset(axis),
     1      size_d1(axis))
          call tce_filename(dfilename(axis),filename)
          call createfile(filename,d_d1(axis),size_d1(axis))
          if (nodezero.and.util_print('mod1',print_default)) then
            write(LuOut,*)
            write(LuOut,9000) 'dipole',size_d1(axis)
            write(LuOut,9010) 'dipole',filename(1:20)
          endif
        enddo
        if (.not.read_integrals(3)) then
          call tce_mod1(g_aod1,d_d1)
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
#endif
c
c       new dipole integrals
c
#if defined(NEWDIPOLE)
        if (multipole.ge.1) then
          do i = 1, multipole
            if (i.eq.1) then
              call hnd_mtpole2(rtdb,ao_bas_han,geom,1,dipole_scf)
            elseif (i.eq.2) then
              call hnd_mtpole2(rtdb,ao_bas_han,geom,2,quadrupole_scf)
            elseif (i.eq.3) then
              call hnd_mtpole2(rtdb,ao_bas_han,geom,3,octupole_scf)
            endif
c
c            print*,'multipole = ', multipole
c            print*,'dipole_scf = ',dipole_scf
c
            j = (i+1)*(i+2)/2 ! 3, 6, 10
            do axis = 1, j
              irrep_d=sym_abelian_axis(geom,axis)
              call sym_irrepname(geom,irrep_d+1,irrepname)
              if (nodezero.and.util_print('mod1',print_default)) then
                write(LuOut,*)
                if (i.eq.1) write(LuOut,9440) axisname(axis),irrepname
                if (i.eq.2) write(LuOut,9440) axisname(axis+3),irrepname
                if (i.eq.3) write(LuOut,9440) axisname(axis+9),irrepname
              endif
              if (i.eq.1) then
                call tce_mod1_offset_new(l_d1_offset(axis),
     1            k_d1_offset(axis),size_d1(axis),irrep_d)
                call tce_filename(dfilename(axis),filename)
                call createfile(filename,d_d1(axis),size_d1(axis))
                if (nodezero.and.util_print('mod1',print_default)) then
                  if (i.eq.1) write(LuOut,*)
                  write(LuOut,9000) 'dipole',size_d1(axis)
                  write(LuOut,9010) 'dipole',filename(1:20)
                endif
              elseif (i.eq.2) then
                call tce_mod1_offset_new(l_q1_offset(axis),
     1            k_q1_offset(axis),size_q1(axis),irrep_d)
                call tce_filename(qfilename(axis),filename)
                call createfile(filename,d_q1(axis),size_q1(axis))
                if (nodezero.and.util_print('moq1',print_default)) then
                  if (i.eq.1) write(LuOut,*)
                  write(LuOut,9000) 'quadrupole',size_q1(axis)
                  write(LuOut,9010) 'quadrupole',filename(1:20)
                endif
              elseif (i.eq.3) then
                call tce_mod1_offset_new(l_o1_offset(axis),
     1            k_o1_offset(axis),size_o1(axis),irrep_d)
                call tce_filename(ofilename(axis),filename)
                call createfile(filename,d_o1(axis),size_o1(axis))
                if (nodezero.and.util_print('moo1',print_default)) then
                  if (i.eq.1) write(LuOut,*)
                  write(LuOut,9000) 'octupole',size_o1(axis)
                  write(LuOut,9010) 'octupole',filename(1:20)
                endif
              endif
            enddo
            if (i.eq.1) then
              if (.not.read_integrals(3)) then
                call tce_mod1_new(rtdb,d_d1,1)
                if (save_integrals(3)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving dipole integrals now...'
                  endif
                  do axis = 1, j
                    irrep_d=sym_abelian_axis(geom,axis)
                    call d1_restart_save(d_d1(axis),k_d1_offset(axis),
     &                   size_d1(axis),axis,handle_d1(axis),irrep_d)
                  enddo
                endif
              elseif (read_integrals(3)) then
                if(nodezero) then
                  write(LuOut,*) 'Reading in dipole integrals now...'
                endif
                do axis = 1, j
                  irrep_d=sym_abelian_axis(geom,axis)
                  call d1_restart_read(d_d1(axis),k_d1_offset(axis),
     &                 size_d1(axis),axis,handle_d1(axis),irrep_d)
                enddo
              endif
              do axis = 1, j
                call gatoeaf(d_d1(axis))
              enddo
            elseif (i.eq.2) then
              if (.not.read_integrals(4)) then
                call tce_mod1_new(rtdb,d_q1,2)
                if (save_integrals(4)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving quadrupole integrals now...'
                  endif
                  do axis = 1, j
                    irrep_d=sym_abelian_axis(geom,axis)
                    call d1_restart_save(d_d1(axis),k_d1_offset(axis),
     &                   size_d1(axis),axis,handle_d1(axis),irrep_d)
                  enddo
                endif
              elseif (read_integrals(4)) then
                if(nodezero) then
                  write(LuOut,*)'Reading in quadrupole integrals now...'
                endif
                do axis = 1, j
                  irrep_d=sym_abelian_axis(geom,axis)
                  call d1_restart_read(d_q1(axis),k_q1_offset(axis),
     &                 size_q1(axis),axis,handle_q1(axis),irrep_d)
                enddo
              endif
              do axis = 1, j
                call gatoeaf(d_q1(axis))
              enddo
            elseif (i.eq.3) then
              if (.not.read_integrals(5)) then
                call tce_mod1_new(rtdb,d_o1,3)
                if (save_integrals(5)) then
                  if(nodezero) then
                    write(LuOut,*) 'Saving octupole integrals now...'
                  endif
                  do axis = 1, j
                    irrep_d=sym_abelian_axis(geom,axis)
                    call d1_restart_save(d_o1(axis),k_o1_offset(axis),
     &                   size_o1(axis),axis,handle_o1(axis),irrep_d)
                  enddo
                endif
              elseif (read_integrals(5)) then
                if(nodezero) then
                  write(LuOut,*) 'Reading in octupole integrals now...'
                endif
                do axis = 1, j
                  irrep_d=sym_abelian_axis(geom,axis)
                  call d1_restart_read(d_o1(axis),k_o1_offset(axis),
     &                 size_o1(axis),axis,handle_o1(axis),irrep_d)
                enddo
              endif
              do axis = 1, j
                call gatoeaf(d_o1(axis))
              enddo
            endif !(i.eq.N)
          enddo !(do i)
        endif !(multipole.ge.1)
#endif
      endif
c
c     1-e integrals
c
      cpu = - util_cpusec()
      wall = - util_wallsec()
      if (.not.rtdb_get(rtdb,'tce:recompf',mt_log,1,recompf))
     1  call errquit('tce_energy: failed reading from rtdb',0,
     2  RTDB_ERR)
      if (recompf.and.(.not.read_integrals(1))) then
        call tce_ao1e_fock2e(rtdb,g_ao1e)
      endif
c
c->d3p975
c      do i = 1, ipol
c        if (.not.ga_destroy(g_movecs(i)))
c     1    call errquit('tce_energy: GA problem',0,GA_ERR)
c      enddo
c<-d3p975
c
      call tce_filename('f1',filename)
      call tce_mo1e_offset(l_f1_offset,k_f1_offset,size_1e)
      call createfile(filename,d_f1,size_1e)
      if (.not.read_integrals(1)) then
        if (recompf) then
          call tce_mo1e(g_ao1e,d_f1,k_f1_offset)
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
        write(LuOut,9010) '1-e',filename(1:20)
c       write(LuOut,9090) '1-e',d_f1
        write(LuOut,9020) cpu, wall
        call util_flush(LuOut)
      endif
#ifdef MIRROR_FOCK
      if (nodezero) write(LuOut,*) 'Mirroring Fock matrix'
c
c     allocate space for local copy
c
      if (.not.ma_push_get(mt_dbl,size_1e,'f1temp',l_f1temp,k_f1temp))
     1    call errquit('tce_energy: MA problem',3737,MA_ERR)
c
c     make a mirror GA for Fock matrix
c
      call tce_filename('f1temp',filename)
      call createfile_mirror(filename,d_f1temp,size_1e)
c
c     copy the entire distributed GA to local buffer
c
      call ga_get(d_f1,1,size_1e,1,1,dbl_mb(k_f1temp),1)
      call ga_sync()
c
c     put it back in the mirrored GA
c
      call ga_put(d_f1temp,1,size_1e,1,1,dbl_mb(k_f1temp),1)
      call ga_sync() ! no need to merge_mirrored since copying entire array
c
c     delete the distributed version
c
      if (.not.ga_destroy(d_f1))
     1    call errquit('tce_energy: GA problem',3737,GA_ERR)
c
c     delete local temporary space
c
      if (.not.ma_pop_stack(l_f1temp))
     1  call errquit('tce_energy: MA problem',3737,MA_ERR)
c
c     reassign the handle so that it points to the mirrored GA
c
      d_f1 = d_f1temp
#endif
c
      io_v2=IOGA
c
      if(.not.intorb) THEN  !--------------
        if (read_integrals(2)) then
          call errquit('tce_energy: cannot restart without 2eorb',
     1                  911,GA_ERR)
        endif
c
c       2-e integrals first half
c
        cpu = - util_cpusec()
        wall = - util_wallsec()
        if(fast2e.eq.1) then
          if (nodezero) then
            write(LuOut,*)
            write(LuOut,*) 'tce_ao2e: fast2e=1'
            write(LuOut,*) 'half-transformed integrals in memory'
            call util_flush(LuOut)
          endif
          call tce_ao2e(rtdb,d_ao2e)
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
c        elseif(fast2e.eq.4) then
c          if (nodezero) then
c            write(LuOut,*)
c            write(LuOut,*) 'tce_ao2e: fast2e=4'
c            write(LuOut,*) 'half-transformed integrals in memory'
c            call util_flush(LuOut)
c          endif
c!           call tce_ao2e_new(rtdb,geom,ao_bas_han,tol2e,ipol,d_ao2e)
c          call tce_ao2e_new(rtdb,ao_bas_han,d_ao2e)
c        endif
        cpu = cpu + util_cpusec()
        wall = wall + util_wallsec()
        if (nodezero.and.util_print('time',print_default))
     1    write(LuOut,9020) cpu, wall
        if (nodezero) call util_flush(LuOut)
c
c       2-e integrals second half
c
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
          call tce_mo2e(rtdb,d_ao2e,d_v2,k_v2_offset)
          call reconcilefile(d_v2,size_2e)
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
          write(LuOut,9010) '2-e',filename(1:20)
c          write(LuOut,9090) '2-e',d_v2
          write(LuOut,9020) cpu, wall
          call util_flush(LuOut)
        endif
      elseif(intorb) THEN !----------------------
        if(nodezero) then
          write(6,*)'4-electron integrals stored in orbital form'
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
          call tce_mo2e_offset_intorb(lax_v2_alpha_offset,
     &         kax_v2_alpha_offset,size_2e)
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
c          The GA copy of 2-e integrals
           call tce_filename('v2ga',filename)
           call createfile(filename,d_v2ga,size_2e_ga)
           call ga_zero(d_v2ga)
           d_v2 = d_v2ga
c          The SF copy of 2-e integrals
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
c
c       ============================================================
c
        if((.not.read_integrals(2)).or.(fast2e.eq.9)) then
c
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
c        elseif(fast2e.eq.7) then
c          if(idisk.ne.1) then
c            call errquit('tce_energy: 2emet 7 with idiskx=1 only',
c     2                                             idisk,CALC_ERR)
c          endif
c          call tce_mo2e_zones_4a_disk_4s_new_patch_m3(rtdb,d_v2,
c     1                                  kax_v2_alpha_offset,size_2e)
c        elseif(fast2e.eq.8) then
c          if(idisk.ne.1) then
c            call errquit('tce_energy: 2emet 8 with idiskx=1 only',
c     2                                             idisk,CALC_ERR)
c          endif
c          call tce_mo2e_zones_4a_disk_4s_new_patch_m4(rtdb,d_v2,
c     1                                  kax_v2_alpha_offset,size_2e)
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
c          if (nodezero) write(LuOut,*) 'tce_mo2e_incore_2eorb'
c          call tce_mo2e_incore_2eorb(rtdb,d_v2ga,kax_v2_alpha_offset,
c     1                               size_2e)
c          if (nodezero) then
c            write(LuOut,*) 'tce_mo2e_zones_4a_disk_4s_new_patch'
c          endif
c          call tce_mo2e_zones_4a_disk_4s_new_patch(rtdb,d_v2ga,
c     1                                  kax_v2_alpha_offset,size_2e)
c          if (nodezero) write(LuOut,*) 'tce_mo2e_disk_2eorb'
c          idisk=.true.
c          call tce_mo2e_disk_2eorb(rtdb,d_v2sf,kax_v2_alpha_offset,
c     1                             size_2e)
          if (nodezero) write(LuOut,*) 'tce_mo2e_hybrid_2eorb_split'
          idisk=1
          call tce_mo2e_hybrid_2eorb_split(d_v2ga,d_v2sf,
     1         kax_v2_alpha_offset)
        elseif(fast2e.eq.13) then
          call tce_mo2e_zones_4a_disk_ga_N5(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
        elseif(fast2e.eq.14) then
          call  tce_mo2e_zones_4a_disk_ga_chop_N5(rtdb,d_v2,
     1                                  kax_v2_alpha_offset,size_2e)
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
c
c       ============================================================
c
        elseif (read_integrals(2).and.(.not.(fast2e.eq.9))) then
          if(nodezero) then
            write(LuOut,*) 'Reading in 2-electron integrals now...'
          endif
          call v2_restart_read(d_v2,k_v2_alpha_offset,size_2e,
     &                         handle_v2)
        endif
c
c       ============================================================
c
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
      end if ! intorb -----------------------
c
c     ====================
c     Initial t amplitudes
c     ====================
c
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
      elseif ((model.eq.'ccsdt').or.(model.eq.'cisdt')) then
        needt1 = .true.
        needt2 = .true.
        needt3 = .true.
      elseif (model.eq.'ccsdta') then
        needt1  = .true.
        needt2  = .true.
        needt3a = .true.
      elseif ((model.eq.'ccsdtq').or.(model.eq.'cisdtq')) then
        needt1 = .true.
        needt2 = .true.
        needt3 = .true.
        needt4 = .true.
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
      if(nodezero) write(LuOut,*) 'do_pt = ',do_pt
      if (.not.rtdb_get(rtdb,'tce:do_lam_pt',mt_log,1,do_lam_pt)) then
        do_lam_pt = .false.
      endif
      if(nodezero) write(LuOut,*) 'do_lam_pt = ',do_lam_pt
      if (.not.rtdb_get(rtdb,'tce:do_cr_pt',mt_log,1,do_cr_pt)) then
        do_cr_pt = .false.
      endif
      if(nodezero) write(LuOut,*) 'do_cr_pt = ',do_cr_pt
      if (.not.rtdb_get(rtdb,'tce:do_lcr_pt',mt_log,1,do_lcr_pt)) then
        do_lcr_pt = .false.
      endif
      if(nodezero) write(LuOut,*) 'do_lcr_pt = ',do_lcr_pt
      if (.not.rtdb_get(rtdb,'tce:do_2t_pt',mt_log,1,do_2t_pt)) then
        do_2t_pt = .false.
      endif
      if(nodezero) write(LuOut,*) 'do_2t_pt = ',do_2t_pt
c
c     t1 amplitudes
c
      if (needt1) then
        cpu = - util_cpusec()
        call tce_filename('t1',filename)
        call tce_t1_offset(l_t1_offset,k_t1_offset,size_t1)
c -------
        if(nodezero) then
         write(LuOut,*)'T1-number-of-tasks',int_mb(k_t1_offset)
         call util_flush(LuOut)
        end if
c -------------
#ifdef MIRROR_SING
        call createfile_mirror(filename,d_t1,size_t1)
#else
        call createfile(filename,d_t1,size_t1)
#endif
        if (.not.read_t(1)) then
          call tce_guess_t1(d_t1,k_t1_offset)
#ifdef MIRROR_SING
          call reconcilefile(d_t1)
#endif
        elseif (read_t(1)) then
          call x1_restart_read(d_t1,k_t1_offset,size_t1,0,
     1                         handle_t1,irrep_t)
        endif
cc        call tce_guess_t1_zero(d_t1,k_t1_offset)
        call reconcilefile(d_t1,size_t1)
        cpu = cpu + util_cpusec()
        if (nodezero.and.util_print('t1',print_default)) then
          write(LuOut,*)
          write(LuOut,9000) 't1',size_t1
          write(LuOut,9010) 't1',filename(1:20)
          write(LuOut,9090) 't1',d_t1
c         write(LuOut,9020) cpu
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
        if (lineresp.or.imagresp) then
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
            if(imagresp) then
              call tce_filename(tr1filename(axis+6),filename)
              call createfile(filename,d_tr1(axis+6),
     &             size_tr1(axis))
              if (read_tr(8+1)) then
                call tr1_restart_read(d_tr1(axis+6),k_tr1_offset(axis),
     1            size_tr1(axis),axis+6,handle_tr1(axis+6),irrep_x)
              elseif (read_tr(1)) then
                call tr1_restart_read(d_tr1(axis+6),k_tr1_offset(axis),
     1            size_tr1(axis),axis,handle_tr1(axis+6),irrep_x)
              endif
              call gatoeaf(d_tr1(axis+6))
            endif
            if (nodezero.and.util_print('tr1',print_debug)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
              write(LuOut,9000) 'tr1',size_tr1(axis)
              write(LuOut,9090) 'tr1',d_tr1(axis)
              write(LuOut,9090) 'tr1',d_tr1(axis+3)
              if(imagresp) then
                write(LuOut,9090) 'tr1',d_tr1(axis+6)
              endif
              call util_flush(LuOut)
            endif
          endif
          enddo ! axis
        endif ! lineresp.or.imagresp
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
c
c     t2 amplitudes
c
      if (needt2) then
c--       if(.not.intorb) then ! ----------
        cpu = - util_cpusec()
        call tce_filename('t2',filename)
        call tce_t2_offset(l_t2_offset,k_t2_offset,size_t2)
c ------------
        if(nodezero) then
         write(LuOut,*)'T2-number-of-boxes',int_mb(k_t2_offset)
         call util_flush(LuOut)
        end if
c -------------
        call createfile(filename,d_t2,size_t2)
        if (.not.read_t(2)) then
          call tce_guess_t2(d_v2,k_v2_offset,d_t2,k_t2_offset)
          call reconcilefile(d_t2,size_t2)
c starts for t1 amplitudes t1(in the second order of MBPT)
          if(needt1.and.((model.eq.'ccsd').or.(model.eq.'multi'))) then
            if ((ioalg.eq.2).and.(io_v2.eq.IOGA)) then
              call t1mp2(d_f1,d_t1,d_t2,d_v2,k_f1_offset,k_t1_offset,
     &                                       k_t2_offset,k_v2_offset)
              call reconcilefile(d_t1,size_t1)
              call tce_guess_t1(d_t1,k_t1_offset)
              call reconcilefile(d_t1,size_t1)
            end if
          end if
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
        cpu = cpu + util_cpusec()
        if (nodezero.and.util_print('t2',print_default)) then
          write(LuOut,*)
          write(LuOut,9000) 't2',size_t2
          write(LuOut,9010) 't2',filename(1:20)
          write(LuOut,9090) 't2',d_t2
c         write(LuOut,9020) cpu
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
c--       else  ! ----
c--        cpu = - util_cpusec()
c--        call tce_filename('t2',filename)
c--        call tce_t2_offset(l_t2_offset,k_t2_offset,size_t2)
c--        call createfile(filename,d_t2,size_t2)
c--        if(nodezero) then
c--          write(LuOut,*)'before tce_guess_t2'
c--          call util_flush(LuOut)
c--        end if
c--        call ga_zero(d_t2)
c--       if(nodezero) then
c--         write(LuOut,*)'after tce_guess_t2'
c--         call util_flush(LuOut)
c--       end if
c--        cpu = cpu + util_cpusec()
c--        if (nodezero.and.util_print('t2',print_default)) then
c--          write(LuOut,*)
c--          write(LuOut,9000) 't2',size_t2
c--          write(LuOut,9010) 't2',filename(1:20)
c--c         write(LuOut,9090) 't2',d_t2
c--c         write(LuOut,9020) cpu
c--          call util_flush(LuOut)
c--        endif
c--        if (left) then
c--          call tce_y2_offset(l_l2_offset,k_l2_offset,size_l2)
c--          call tce_filename('lambda2',filename)
c--          call createfile(filename,d_lambda2,size_l2)
c--          call gatoeaf(d_lambda2)
c--        endif
c--       endif !---
        if (lineresp.or.imagresp) then
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
            if (imagresp) then
              call tce_filename(tr2filename(axis+6),filename)
              call createfile(filename,d_tr2(axis+6),
     &             size_tr2(axis))
              if (read_tr(8+2)) then
                call tr2_restart_read(d_tr2(axis+6),k_tr2_offset(axis),
     1            size_tr2(axis),axis+6,handle_tr2(axis+6),irrep_x)
              elseif (read_tr(2)) then
                call tr2_restart_read(d_tr2(axis+6),k_tr2_offset(axis),
     1            size_tr2(axis),axis,handle_tr2(axis+6),irrep_x)
              endif
            endif
            call gatoeaf(d_tr2(axis+6))
            if (nodezero.and.util_print('tr2',print_debug)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
              write(LuOut,9000) 'tr2',size_tr1(axis)
              write(LuOut,9090) 'tr2',d_tr1(axis)
              write(LuOut,9090) 'tr2',d_tr1(axis+3)
              if(imagresp) then
                write(LuOut,9090) 'tr2',d_tr1(axis+6)
              endif
              call util_flush(LuOut)
            endif
          endif
          enddo ! axis
        endif ! lineresp.or.imagresp
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
c
c     t3 amplitudes
c
      if (needt3) then
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
          write(LuOut,9010) 't3',filename(1:20)
c         write(LuOut,9090) 't3',d_t3
c         write(LuOut,9020) cpu
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
        if (lineresp.or.imagresp) then
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
            if (imagresp) then
              call tce_filename(tr3filename(axis+6),filename)
              call createfile(filename,d_tr3(axis+6),
     &             size_tr3(axis))
              if (read_tr(8+3)) then
                call tr3_restart_read(d_tr3(axis+6),k_tr3_offset(axis),
     1            size_tr3(axis),axis+6,handle_tr3(axis+6),irrep_x)
              elseif (read_tr(3)) then
                call tr3_restart_read(d_tr3(axis+6),k_tr3_offset(axis),
     1            size_tr3(axis),axis,handle_tr3(axis+6),irrep_x)
              endif
              call gatoeaf(d_tr3(axis+6))
            endif
          endif
          enddo ! axis
        endif ! lineresp.or.imagresp
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
      endif
c
c     t3a amplitudes
c
      if (needt3a) then
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
          write(LuOut,9010) 't3a',filename(1:20)
c         write(LuOut,9090) 't3a',d_t3
c         write(LuOut,9020) cpu
          call util_flush(LuOut)
        endif
      endif
c
c     t4 amplitudes
c
      if (needt4) then
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
          write(LuOut,9010) 't4',filename(1:20)
c         write(LuOut,9090) 't4',d_t4
c         write(LuOut,9020) cpu
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
        if (lineresp.or.imagresp) then
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
            if (imagresp) then
              call tce_filename(tr4filename(axis+6),filename)
              call createfile(filename,d_tr4(axis+6),
     &             size_tr4(axis))
              if (read_tr(8+4)) then
                call tr4_restart_read(d_tr4(axis+6),k_tr4_offset(axis),
     1            size_tr4(axis),axis+6,handle_tr4(axis+6),irrep_x)
              elseif (read_tr(4)) then
                call tr4_restart_read(d_tr4(axis+6),k_tr4_offset(axis),
     1             size_tr4(axis),axis,handle_tr4(axis+6),irrep_x)
              endif
            endif
            call gatoeaf(d_tr4(axis+6))
          endif
          enddo ! axis
        endif ! lineresp.or.imagresp
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
      endif
c
c     =======================
c     Create residual offsets
c     =======================
c
      call tce_e_offset(l_e_offset,k_e_offset,size_e)
      if (left) then
        call tce_e_offset(l_d0_offset,k_d0_offset,size_d0)
        if (multipole.ge.2) then
          call tce_e_offset(l_q0_offset,k_q0_offset,size_q0)
        endif
        if (multipole.ge.3) then
          call tce_e_offset(l_o0_offset,k_o0_offset,size_o0)
        endif
      endif
      if (lineresp.or.imagresp) then
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
       if (needt2) call tce_t2_offset(l_r2_offset,k_r2_offset,size_r2)
       if (needt3) call tce_t3_offset(l_r3_offset,k_r3_offset,size_r3)
       if (needt3a) call tce_t3a_offset(l_r3_offset,k_r3_offset,size_r3)
       if (needt4) call tce_t4_offset(l_r4_offset,k_r4_offset,size_r4)
      endif
c
c     ======================
c     Create auxiliary files
c     ======================
c
      if (neede2) then
        call tce_filename('e2',filename)
        call tce_e_offset(l_e2_offset,k_e2_offset,size_e2)
        call createfile(filename,d_e2,size_e2)
        call gatoeaf(d_e2)
c       write(LuOut,9090) 'e2',d_e2
      endif
      if (needu1) then
        call tce_filename('u1',filename)
        call tce_t1_offset(l_u1_offset,k_u1_offset,size_u1)
        call createfile(filename,d_u1,size_u1)
        call tce_guess_t1(d_u1,k_u1_offset)
        call reconcilefile(d_u1,size_u1)
c       write(LuOut,9090) 'u2',d_u1
      endif
      if (needu2) then
c--       if(.not.intorb) then
        call tce_filename('u2',filename)
        call tce_t2_offset(l_u2_offset,k_u2_offset,size_u2)
        call createfile(filename,d_u2,size_u2)
        call tce_guess_t2(d_v2,k_v2_offset,d_u2,k_u2_offset)
        call reconcilefile(d_u2,size_u2)
c       write(LuOut,9090) 'u2',d_u2
c--       else
c--        call tce_filename('u2',filename)
c--        call tce_t2_offset(l_u2_offset,k_u2_offset,size_u2)
c--        call createfile(filename,d_u2,size_u2)
c--c        call tce_guess_t2(d_v2,k_v2_offset,d_u2,k_u2_offset)
c--        call ga_zero(d_u2)
c--        call reconcilefile(d_u2,size_u2)
c--c       write(LuOut,9090) 'u2',d_u2
c--       endif
      endif
      if (needw1) then
        call tce_filename('w1',filename)
        call tce_t1_offset(l_w1_offset,k_w1_offset,size_w1)
        call createfile(filename,d_w1,size_w1)
        call gatoeaf(d_w1)
c       write(LuOut,9090) 'w1',d_w1
      endif
      if (needw2) then
        call tce_filename('w2',filename)
        call tce_t2_offset(l_w2_offset,k_w2_offset,size_w2)
        call createfile(filename,d_w2,size_w2)
        call gatoeaf(d_w2)
c       write(LuOut,9090) 'w2',d_w2
      endif
      if (needw3) then
        call tce_filename('w3',filename)
        call tce_t3_offset(l_w3_offset,k_w3_offset,size_w3)
        call createfile(filename,d_w3,size_w3)
        call gatoeaf(d_w3)
c       write(LuOut,9090) 'w3',d_w3
      endif
      if (needw4) then
        call tce_filename('w4',filename)
        call tce_t4_offset(l_w4_offset,k_w4_offset,size_w4)
        call createfile(filename,d_w4,size_w4)
        call gatoeaf(d_w4)
c       write(LuOut,9090) 'w4',d_w4
      endif
c
c     =========================
c     CC / CI / MBPT iterations
c     =========================
c
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
c -------------
c CCD     right
c -------------
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
c -------------
c LCCD    right
c -------------
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
      if ((model.eq.'ccsd').or.do_ccsd) then
c -------------
c CCSD    right
c -------------
#ifdef MIRROR_SING
        call ccsd_energy_noloc(d_e,d_f1,d_v2,d_t1,d_t2,
     1       k_e_offset,k_f1_offset,k_v2_offset,
     2       k_t1_offset,k_t2_offset,
     3       size_t1,size_t2,ref,corr)
#else
        call ccsd_energy_loc(d_e,d_f1,d_v2,d_t1,d_t2,
     1       k_e_offset,k_f1_offset,k_v2_offset,
     2       k_t1_offset,k_t2_offset,
     3       size_t1,size_t2,ref,corr)
#endif
c -------------
c CCSD    left
c -------------
        if (left) then
          call ccsd_lambda(d_d0,d_f1,d_v2,d_d1,d_t1,d_t2,
     1         d_lambda1,d_lambda2,
     2         k_d0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     3         k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset,
     4         size_l1,size_l2,dipole_scf,dipole_cor)
c
c->d3p975
c
c calculated ccsd 1-particle reduced density matrix/ces
c
            if (nroots.eq.0.and.idens) then
               if (.not.intorb)
     1            call ccsd_1prdm(rtdb,d_t1,d_t2,d_lambda1,d_lambda2,
     1                            k_t1_offset,k_t2_offset,
     2                            k_l1_offset,k_l2_offset)
            endif
c
c<-d3p975
c
        endif ! left
c -------------
c CCSD-IR
c -------------
 2234   continue
        if (imagresp) then
          if(nodezero) then
            write(6,*) 'inumfreq = ',inumfreq
            write(6,*) 'Casimir-Polder integration points are:'
            do i = inumfreq,1,-1
              write(6,999) i,ifreq(i)
            enddo
            call util_flush(LuOut)
          endif
  999 format(3x,'ifreq(',i3,') = ',f14.8)
c
c         REAL COMPONENT
c
          do omegacount=inumfreq,1,-1
            omega = ifreq(omegacount)
c            if (omega.eq.(0.0d0)) then
c              call errquit('tce_energy: imag(omega) is zero',0,CALC_ERR)
c            endif
c
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
c
c           offset to the third set of response for real component
            dynaxis = 6
c
            if (nodezero) write(LuOut,9431) omega
c
            if (guess_ir_real.and.(omegacount.eq.1)) then
c
            if (nodezero) write(6,9121) 'Initial guess x = b/Adiag'
c
c           0. Create dummies
c
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
              call tce_zero(d_rr1(axis),size_tr1(axis))
c
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
              call tce_zero(d_rr2(axis),size_tr2(axis))
c
              if (nodezero.and.(util_print('rr1',print_debug)
     1                      .or.util_print('rr2',print_debug))) then
                write(LuOut,*)
                write(LuOut,9440) axisname(axis),irrepname
                write(LuOut,9000) 'rr1',size_tr1(axis)
                write(LuOut,9000) 'rr2',size_tr2(axis)
                write(LuOut,9090) 'rr1',d_rr1(axis)
                write(LuOut,9090) 'rr2',d_rr2(axis)
                call util_flush(LuOut)
              endif
c
c           1. Form b
c
              call ccsd_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     1             k_tr1_offset(axis),k_d1_offset(axis),
     2             k_t1_offset,k_t2_offset)
              call ccsd_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,
     1             k_tr2_offset(axis),k_d1_offset(axis),
     2             k_t1_offset,k_t2_offset)
c
c           2. Hit with preconditioner
c
              call tce_jacobi_tr1(d_rr1(axis),d_tr1(axis+6),
     1                            k_tr1_offset(axis),0.0d0,0.0d0)
              call tce_jacobi_tr2(d_rr2(axis),d_tr2(axis+6),
     1                            k_tr2_offset(axis),0.0d0,0.0d0)
c
c           4. Delete dummies
c
              call deletefile(d_rr1(axis))
              call deletefile(d_rr2(axis))
c
            endif ! guess_ir_real
c
            call tce_diis_init()
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSD-IR (real component)"
c
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
              call tce_zero(d_rr1(axis),size_tr1(axis))
c
              call tce_filename(ir1filename(axis),filename)
              call createfile(filename,d_ir1(axis),size_tr1(axis))
              call tce_zero(d_ir1(axis),size_tr1(axis))
c
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
              call tce_zero(d_rr2(axis),size_tr2(axis))
c
              call tce_filename(ir2filename(axis),filename)
              call createfile(filename,d_ir2(axis),size_tr2(axis))
              call tce_zero(d_ir2(axis),size_tr2(axis))
c
              if (nodezero.and.(util_print('rr1',print_debug)
     1                      .or.util_print('rr2',print_debug))) then
                write(LuOut,*)
                write(LuOut,9440) axisname(axis),irrepname
                write(LuOut,9000) 'rr1',size_tr1(axis)
                write(LuOut,9090) 'rr1',d_rr1(axis)
                write(LuOut,9090) 'ir1',d_ir1(axis)
                write(LuOut,9000) 'rr2',size_tr2(axis)
                write(LuOut,9090) 'rr2',d_rr2(axis)
                write(LuOut,9090) 'ir2',d_ir2(axis)
                call util_flush(LuOut)
              endif
c
              call ccsd_o1(d_ir1(axis),d_d1(axis),d_t1,d_t2,
     1             k_tr1_offset(axis),k_d1_offset(axis),
     2             k_t1_offset,k_t2_offset)
              call eomccsd_x1_old(d_f1,d_ir1(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),k_f1_offset,
     2             k_tr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis))
c
              call ccsd_o2(d_ir2(axis),d_d1(axis),d_t1,d_t2,
     1             k_tr2_offset(axis),k_d1_offset(axis),
     2             k_t1_offset,k_t2_offset)
              call eomccsd_x2_old(d_f1,d_ir2(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             k_f1_offset,k_tr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             size_tr1(axis),size_tr2(axis))
c
              call reconcilefile(d_ir1(axis),size_tr1(axis))
              call reconcilefile(d_ir2(axis),size_tr2(axis))
c
              call eomccsd_x1_old(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1             d_ir1(axis),d_ir2(axis),k_f1_offset,
     2             k_tr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis))
c
              call eomccsd_x2_old(d_f1,d_rr2(axis),d_t1,d_t2,d_v2,
     1             d_ir1(axis),d_ir2(axis),
     2             k_f1_offset,k_tr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             size_tr1(axis),size_tr2(axis))
c
              call reconcilefile(d_rr1(axis),size_tr1(axis))
              call reconcilefile(d_rr2(axis),size_tr2(axis))
c
              call deletefile(d_ir1(axis))
              call deletefile(d_ir2(axis))
c
              call daxpyfile(1,omega*omega,d_tr1(axis+dynaxis),
     1             d_rr1(axis),size_tr1(axis))
              call reconcilefile(d_rr1(axis),size_tr1(axis))
              call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c
              call daxpyfile(1,omega*omega,d_tr2(axis+dynaxis),
     1             d_rr2(axis),size_tr2(axis))
              call reconcilefile(d_rr2(axis),size_tr2(axis))
              call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c
              residual = max(rr1,rr2)
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
                  if (nodezero) then
                    write(LuOut,9082) "T(1) singles",rr1
                    write(LuOut,9082) "T(1) doubles",rr2
                  endif
                endif
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
                goto 2235
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
              call tce_diis3(.false.,iter,.true.,.true.,.false.,.false.,
     1             d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2             size_tr1(axis),
     3             d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4             size_tr2(axis),
     5             dummy,dummy,dummy,dummy,
     6             dummy,dummy,dummy,dummy,omega,2)
              call deletefile(d_rr2(axis))
              call deletefile(d_rr1(axis))
              if (nodezero) call util_flush(LuOut)
            enddo ! iter loop
          call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 2235 continue
c
          endif ! respaxis(axis)
          enddo ! axis loop
c
c         IMAGINARY COMPONENT
c
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
c
          do omegasign=1,1 ! temporary 2
            if (omegasign.eq.1) then
              omega = (1.0d0)*ifreq(omegacount)
              dynaxis = 0
            elseif (omegasign.eq.2) then
              omega = (-1.0d0)*ifreq(omegacount)
              dynaxis = 3
              call daxpyfile(1,-1.0d0,d_tr1(axis+0),
     1             d_tr1(axis+3),size_tr1(axis))
              call daxpyfile(1,-1.0d0,d_tr2(axis+0),
     1             d_tr2(axis+3),size_tr2(axis))
            endif
            if (nodezero) write(LuOut,9431) omega
c
            call tce_diis_init()
c
            if (guess_ir_imag.and.(omegacount.eq.1)) then
c
c           Initial guess x_I = omega*x_R/(Adiag-omega)
           if (nodezero) write(6,9121) 
     1        'Initial guess x_I = omega*x_R/Adiag'
c
c           0. Create dummies
c
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
              call tce_zero(d_rr1(axis),size_tr1(axis))
c
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
              call tce_zero(d_rr2(axis),size_tr2(axis))
c
c           1. Form w_I = omega*x_R
c
              call daxfile(1,omega,d_tr1(axis+6),d_rr1(axis),
     1             size_tr1(axis))
              call daxfile(1,omega,d_tr2(axis+6),d_rr2(axis),
     1             size_tr2(axis))
c
c           3. Hit with preconditioner
c
              call tce_jacobi_tr1(d_rr1(axis),d_tr1(axis+dynaxis),
     1             k_tr1_offset(axis),0.0d0,0.0d0)
              call tce_jacobi_tr2(d_rr2(axis),d_tr2(axis+dynaxis),
     1             k_tr2_offset(axis),0.0d0,0.0d0)
c
c           4. Delete dummies
c
              call deletefile(d_rr1(axis))
              call deletefile(d_rr2(axis))
c
            endif ! guess_ir_imag
c
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSD-IR (imaginary component)"
c
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
c              call tce_zero(d_rr1(axis),size_tr1(axis))
              call daxfile(1,(-1.0d0)*omega,d_tr1(axis+6),
     1             d_rr1(axis),size_tr1(axis))
              call eomccsd_x1_old(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             k_f1_offset,k_tr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis))
              call reconcilefile(d_rr1(axis),size_tr1(axis))
              call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
c              call tce_zero(d_rr2(axis),size_tr2(axis))
              call daxfile(1,(-1.0d0)*omega,d_tr2(axis+6),
     1             d_rr2(axis),size_tr2(axis))
              call eomccsd_x2_old(d_f1,d_rr2(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             k_f1_offset,k_tr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             size_tr1(axis),size_tr2(axis))
              call reconcilefile(d_rr2(axis),size_tr2(axis))
              call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
              residual = max(rr1,rr2)
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
                  if (nodezero) then
                    write(LuOut,9082) "T(1) singles",rr1
                    write(LuOut,9082) "T(1) doubles",rr2
                  endif
                endif
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
                goto 2236
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
              call tce_diis3(.false.,iter,.true.,.true.,.false.,.false.,
     1             d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2             size_tr1(axis),
     3             d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4             size_tr2(axis),
     5             dummy,dummy,dummy,dummy,
     6             dummy,dummy,dummy,dummy,0.0d0,1)
              call deletefile(d_rr2(axis))
              call deletefile(d_rr1(axis))
              if (nodezero) call util_flush(LuOut)
            enddo ! iter loop
          call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 2236 continue
c
          enddo ! omegasign loop
c
          endif ! respaxis(axis)
          enddo ! axis loop
c
c CCSD-IR evaluation step
c
          cpu=-util_cpusec()
          wall=-util_wallsec()
          do axisA = 1, 3
          do axisB = 1, axisA
          alpha(axisA,axisB)=0.0d0
          if (respaxis(axisA).and.respaxis(axisB)) then
            irrep_a=sym_abelian_axis(geom,axisA)
            irrep_b=sym_abelian_axis(geom,axisB)
            irrep_y=irrep_g
c
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,size_a0)
c
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "axisA = ",axisA
            if (nodezero) write(LuOut,*) "axisB = ",axisB
#endif
c
            alpha1=0.0d0
            alpha2=0.0d0
c
            irrep_d=irrep_a
            irrep_tr=irrep_b
c
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha_1 AB"
#endif
            call alpha_1(d_d1(axisA),d_a0,d_t1,d_t2,
     1           d_tr1(axisB+6),d_tr2(axisB+6),
     2           d_lambda1,d_lambda2,
     3           k_d1_offset(axisA),k_a0_offset,
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisB),k_tr2_offset(axisB),
     6           k_l1_offset,k_l2_offset)
c
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha1,1,0)
            call tce_zero(d_a0,size_a0)
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha_1 AB = ",alpha1
#endif
c
            if (axisA.eq.axisB) then
              alpha2=alpha1
            else
              irrep_d=irrep_b
              irrep_tr=irrep_a
#ifdef DEBUG_PRINT
              if (nodezero) write(LuOut,*) "alpha_1 BA"
#endif
              call alpha_1(d_d1(axisB),d_a0,d_t1,d_t2,
     1             d_tr1(axisA+6),d_tr2(axisA+6),
     2             d_lambda1,d_lambda2,
     3             k_d1_offset(axisB),k_a0_offset,
     4             k_t1_offset,k_t2_offset,
     5             k_tr1_offset(axisA),k_tr2_offset(axisA),
     6             k_l1_offset,k_l2_offset)
c
              call reconcilefile(d_a0,size_a0)
              call get_block(d_a0,alpha2,1,0)
              call tce_zero(d_a0,size_a0)
            endif ! axisA.eq.axisB
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha_1 BA = ",alpha2
#endif
c
            alpha(axisA,axisB)=alpha(axisA,axisB)+alpha1+alpha2
c
            alpha1=0.0d0
            alpha2=0.0d0
            alpha3=0.0d0
c
            irrep_tra=irrep_a
            irrep_trb=irrep_b
c
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha_2 A(+i) B(-i)"
#endif
            call alpha_2(d_f1,d_a0,d_t1,d_t2,
     1           d_tr1(axisA),d_tr2(axisA),
c     2           d_tr1(axisB+3),d_tr2(axisB+3), ! temporary
     2           d_tr1(axisB),d_tr2(axisB),
     3           d_v2,d_lambda1,d_lambda2,k_f1_offset,k_a0_offset,
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisA),k_tr2_offset(axisA),
     6           k_tr1_offset(axisB),k_tr2_offset(axisB),k_v2_offset,
     7           k_l1_offset,k_l2_offset)
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha1,1,0)
            call tce_zero(d_a0,size_a0)
            alpha1 = -alpha1 ! temporary
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha_2 A(+i) B(-i) = ",alpha1
#endif
c
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha_2 A(-i) B(+i)"
#endif
            call alpha_2(d_f1,d_a0,d_t1,d_t2,
c     1           d_tr1(axisA+3),d_tr2(axisA+3), ! temporary
     1           d_tr1(axisA),d_tr2(axisA),
     2           d_tr1(axisB),d_tr2(axisB),
     3           d_v2,d_lambda1,d_lambda2,k_f1_offset,k_a0_offset,
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisA),k_tr2_offset(axisA),
     6           k_tr1_offset(axisB),k_tr2_offset(axisB),k_v2_offset,
     7           k_l1_offset,k_l2_offset)
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha2,1,0)
            call tce_zero(d_a0,size_a0)
            alpha2 = -alpha2 ! temporary
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha_2 A(-i) B(+i) = ",alpha2
#endif
c
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha_2 A(r) B(r)"
#endif
            call alpha_2(d_f1,d_a0,d_t1,d_t2,
     1           d_tr1(axisA+6),d_tr2(axisA+6),
     2           d_tr1(axisB+6),d_tr2(axisB+6),
     3           d_v2,d_lambda1,d_lambda2,k_f1_offset,k_a0_offset,
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisA),k_tr2_offset(axisA),
     6           k_tr1_offset(axisB),k_tr2_offset(axisB),k_v2_offset,
     7           k_l1_offset,k_l2_offset)
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha3,1,0)
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha_2 A(r) B(r) = ",alpha3
#endif
            call deletefile(d_a0)
c
            alpha(axisA,axisB)=alpha(axisA,axisB)-alpha1-alpha2+alpha3
            alpha(axisA,axisB)=-1.0d0*alpha(axisA,axisB)
c
#ifdef DEBUG_PRINT
            if (nodezero) write(LuOut,*) "alpha(axisA,axisB) = ",
     1                                    alpha(axisA,axisB)
            call util_flush(LuOut)
#endif
c
          endif ! respaxis(axis)
          enddo ! axisB loop
          enddo ! axisA loop
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
c
          do i = 1, 3
            do j = 1, i
              alphacopy(i,j)=alpha(i,j)
              alphacopy(j,i)=alpha(i,j)
            enddo
          enddo
c
          call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
c
          alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
          alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
          alphaani = dsqrt(0.5d0*alphaani)
c
          if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1        .and.nodezero) write(LuOut,9911)
c
          if (nodezero) write(LuOut,9434) "CCSD Imaginary Response",
     1      ifreq(omegacount),auXnm/ifreq(omegacount),
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
c
          ifreqval(omegacount,1) = alpha(1,1)
          ifreqval(omegacount,2) = alpha(2,2)
          ifreqval(omegacount,3) = alpha(3,3)
          ifreqval(omegacount,4) = alpha(2,1)
          ifreqval(omegacount,5) = alpha(3,1)
          ifreqval(omegacount,6) = alpha(3,2)
          ifreqval(omegacount,7) = alphaiso
          ifreqval(omegacount,8) = alphaani
c
          enddo ! omegacount loop
c
c       Casimir-Polder integration for determination of C6 coefficient
c
          do i = 1, 8
            integral(i) = 0.0d0
            do omegacount = 1,inumfreq-1
              num1 = ifreqval(omegacount+1,i)
              num0 = ifreqval(omegacount,i)
              den1 = ifreq(omegacount+1)
              den0 = ifreq(omegacount)
              num1sq = num1*num1
              num0sq = num0*num0
              numtot = num1sq+num0sq
              dentot = den1 - den0
              integral(i) = integral(i) + 0.5d0*numtot*dentot
            enddo
            integral(i) = integral(i)*3.0d0/pi
          enddo
          if (nodezero) write(LuOut,9435) "CCSD Imaginary Response",
     &       integral(1),integral(2),integral(3),
     &       integral(4),integral(5),integral(6),
     &       integral(7),integral(8)
          call util_flush(LuOut)
c
        endif ! imagresp
 2237   continue
c -------------
c CCSD-LR
c -------------
        if (lineresp) then
          call ccsd_lr_alpha(d_d0,d_a0,d_f1,d_v2,d_d1,
     1         d_t1,d_t2,d_lambda1,d_lambda2,d_tr1,d_tr2,
     2         k_d0_offset,k_a0_offset,
     3         k_f1_offset,k_v2_offset,k_d1_offset,
     4         k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset,
     5         k_tr1_offset,k_tr2_offset,
     6         size_tr1,size_tr2,dipole_scf,dipole_cor)
c
        endif ! lineresp
c -------------
c CCSD-YR
c -------------
        if (leftresp) then
          call ccsd_qr_beta(d_a0,d_b0,d_f1,d_v2,d_d1,
     1         d_t1,d_t2,d_lambda1,d_lambda2,d_yr1,d_yr2,d_tr1,d_tr2,
     2         k_a0_offset,k_b0_offset,
     3         k_f1_offset,k_v2_offset,k_d1_offset,
     4         k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset,
     5         k_tr1_offset,k_tr2_offset,k_yr1_offset,k_yr2_offset,
     6         size_yr1,size_yr2)
c
        endif ! leftresp
c -------------------
c CCSD(T) and CCSD[T]
c -------------------
        if (do_pt.or.((perturbative.eq.'(t)').or.
     1                (perturbative.eq.'[t]'))) then
c
          if (nodezero) write(LuOut,*) 'CCSD(T)'
          if (do_pt) perturbative='(t)'
c
          cpu = - util_cpusec()
          wall = - util_wallsec()
          if (restart_ccsd_t.eq.0) then
            if (fast_ccsd_t) then
              if (nodezero) write(LuOut,*) 'Using tuned CCSD(T) code'
              call ccsd_t_vector(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                           d_v2,k_v2_offset,pt3_2,size_t1)
            else
              if (nodezero) write(LuOut,*) 'Using plain CCSD(T) code'
              call ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                    d_v2,k_v2_offset,pt3_1,pt3_2,size_t1)
            endif
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
            if ((restart_ccsd_t.eq.0).and.(.not.fast_ccsd_t)) then
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
c ---------------------------------
c Lambda-CCSD(T) and Lambda-CCSD[T]
c ---------------------------------
        if (do_lam_pt.or.((perturbative.eq.'lambda(t)').or.
     &                    (perturbative.eq.'lambda[t]'))) then
c
          if (nodezero) write(LuOut,*) 'Lambda-CCSD(T)'
          if (do_lam_pt) perturbative='lambda(t)'
c
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
c -------------------------
c CR-CCSD(T) and CR-CCSD[T]
c -------------------------
        if(do_cr_pt.or.((perturbative.eq.'cr_(t)')
     1             .or. (perturbative.eq.'cr_[t]'))) then
c
          if (nodezero) write(LuOut,*) 'CR-CCSD(T)'
          if (do_cr_pt) perturbative='cr_(t)'
c
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
c ----------------------------------------------------
c LR-CCSD(T) and/or LR-CCSD(T)-1 and/or LR-CCSD(TQ)-1p
c ----------------------------------------------------
        if (do_lcr_pt.or.((perturbative.eq.'lr_(t)').or.
     &                    (perturbative.eq.'lr_(tq1p)').or.
     &                    (perturbative.eq.'lr_(tq1)'))) then
c
          if (nodezero) write(LuOut,*) 'LR-CCSD(T)'
          if (do_lcr_pt) perturbative='lr_(t)'
c
          cpu = - util_cpusec()
          wall = - util_wallsec()
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call lr_ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                   d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                   d_e,k_e_offset,pt3_1,pt3_2,
     3                   pt3_3,pt3_4,pt3_5,pt3_6,size_t1)
          call deletefile(d_e)
c
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
c --------------------------------------
c LR-CCSD(TQ)-1p
c --------------------------------------
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
c --------------------------------------------------------
c LR-CCSD(TQ)-1
c --------------------------------------------------------
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
c ---------------------
c CCSD(2)_T and CCSD(2)
c ---------------------
        if (do_2t_pt.or.((perturbative.eq.'2_t').or.
     1                   (perturbative.eq.'2_tq').or.
     2                   (perturbative.eq.'2_tqf'))) then
c
          if (nodezero) write(LuOut,*) 'CCSD(2)_T'
          if (do_2t_pt) perturbative='2_t'
c
          cpu = - util_cpusec()
          wall = - util_wallsec()
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call lr_ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                   d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                   d_e,k_e_offset,pt3_1,pt3_2,
     3                   pt3_3,pt3_4,pt3_5,pt3_6)
          call deletefile(d_e)
c
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
       endif
      endif
      if ((model.eq.'lccsd').or.do_lccsd) then
c -------------
c LCCSD   right
c -------------
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
c -------------
c CCSDT   right
c -------------
        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CCSDT"
          call tce_filename('e',filename)
          call createfile(filename,d_e,size_e)
          call ccsdt_e(d_f1,d_e,d_t1,d_t2,d_v2,
     1                 k_f1_offset,k_e_offset,
     2                 k_t1_offset,k_t2_offset,k_v2_offset)
          call reconcilefile(d_e,size_e)
          call tce_filename('r1',filename)
          call createfile(filename,d_r1,size_t1)
          call ccsdt_t1(d_f1,d_r1,d_t1,d_t2,d_t3,d_v2,
     1                  k_f1_offset,k_t1_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset)
          call reconcilefile(d_r1,size_t1)
          call tce_filename('r2',filename)
          call createfile(filename,d_r2,size_t2)
          call ccsdt_t2(d_f1,d_r2,d_t1,d_t2,d_t3,d_v2,
     1                  k_f1_offset,k_t2_offset,k_t1_offset,
     2                  k_t2_offset,k_t3_offset,k_v2_offset,size_t2)
          call reconcilefile(d_r2,size_t2)
          call tce_filename('r3',filename)
          call createfile(filename,d_r3,size_t3)
          call ccsdt_t3(d_f1,d_r3,d_t1,d_t2,d_t3,d_v2,
     1                  k_f1_offset,k_t3_offset,k_t1_offset,
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
              write(LuOut,9070) "CCSDT",corr
              write(LuOut,9080) "CCSDT",ref + corr
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
            goto 2327
cccc            ccsdtleft=.true.
cccc            goto 100
cedo            goto 330
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

 2327 continue
c -------------
c CCSDT   left
c -------------
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
ccc            goto 2010
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
c ----------
c CCSDT(2)_Q
c ----------
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
c -------------
c CCSDT-IR
c -------------
        if (imagresp) then
#if !defined(CCSDTLR)
      CALL ERRQUIT('CCSDTLR not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
c
c         REAL COMPONENT
c
          do omegacount=1,inumfreq
            omega = ifreq(omegacount)
c            if (omega.eq.(0.0d0)) then
c              call errquit('tce_energy: imag(omega) is zero',0,CALC_ERR)
c            endif
c
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
c
c           offset to the third set of response for real component
            dynaxis = 6
c
            if (nodezero) write(LuOut,9431) omega
c
            call tce_diis_init()
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSDT-IR (real component)"
c        write(LuOut,*) "singles initialization"
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
              call tce_zero(d_rr1(axis),size_tr1(axis))
              call tce_filename(ir1filename(axis),filename)
              call createfile(filename,d_ir1(axis),size_tr1(axis))
              call tce_zero(d_ir1(axis),size_tr1(axis))
c        write(LuOut,*) "doubles initialization"
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
              call tce_zero(d_rr2(axis),size_tr2(axis))
              call tce_filename(ir2filename(axis),filename)
              call createfile(filename,d_ir2(axis),size_tr2(axis))
              call tce_zero(d_ir2(axis),size_tr2(axis))
c        write(LuOut,*) "triples initialization"
              call tce_filename(rr3filename(axis),filename)
              call createfile(filename,d_rr3(axis),size_tr3(axis))
              call tce_zero(d_rr3(axis),size_tr3(axis))
              call tce_filename(ir3filename(axis),filename)
              call createfile(filename,d_ir3(axis),size_tr3(axis))
              call tce_zero(d_ir3(axis),size_tr3(axis))
c        write(LuOut,*) "ccsdt_o1"
              call ccsdt_o1(d_ir1(axis),d_d1(axis),d_t1,d_t2,
     &             k_tr1_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset)
c        write(LuOut,*) "ccsdt_o2"
              call ccsdt_o2(d_ir2(axis),d_d1(axis),d_t1,d_t2,d_t3,
     &             k_tr2_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset)
c        write(LuOut,*) "ccsdt_o3"
              call ccsdt_o3(d_ir3(axis),d_d1(axis),d_t1,d_t2,d_t3,
     &             k_tr3_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset)
c        write(LuOut,*) "eomccsdt_x1"
              call eomccsdt_x1(d_f1,d_ir1(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
c        write(LuOut,*) "eomccsdt_x2"
              call eomccsdt_x2(d_f1,d_ir2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
c        write(LuOut,*) "eomccsdt_x3"
              call eomccsdt_x3(d_f1,d_ir3(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr3_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
c        write(LuOut,*) "reconcilefile"
              call reconcilefile(d_ir1(axis),size_tr1(axis))
              call reconcilefile(d_ir2(axis),size_tr2(axis))
              call reconcilefile(d_ir3(axis),size_tr3(axis))
c        write(LuOut,*) "eomccsdt_x1"
              call eomccsdt_x1(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1             d_ir1(axis),d_ir2(axis),d_ir3(axis),
     2             k_f1_offset,k_tr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
c        write(LuOut,*) "eomccsdt_x2"
              call eomccsdt_x2(d_f1,d_rr2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_ir1(axis),d_ir2(axis),d_ir3(axis),
     2             k_f1_offset,k_tr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
c        write(LuOut,*) "eomccsdt_x3"
              call eomccsdt_x3(d_f1,d_rr3(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_ir1(axis),d_ir2(axis),d_ir3(axis),
     2             k_f1_offset,k_tr3_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
c        write(LuOut,*) "deletefile"
              call deletefile(d_ir1(axis))
              call deletefile(d_ir2(axis))
              call deletefile(d_ir3(axis))
c        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,omega*omega,d_tr1(axis+dynaxis),
     1             d_rr1(axis),size_tr1(axis))
              call daxpyfile(1,omega*omega,d_tr2(axis+dynaxis),
     1             d_rr2(axis),size_tr2(axis))
              call daxpyfile(1,omega*omega,d_tr3(axis+dynaxis),
     1             d_rr3(axis),size_tr3(axis))
c        write(LuOut,*) "reconcilefile"
              call reconcilefile(d_rr1(axis),size_tr1(axis))
              call reconcilefile(d_rr2(axis),size_tr2(axis))
              call reconcilefile(d_rr3(axis),size_tr3(axis))
c        write(LuOut,*) "tce_residual_tr1"
              call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c        write(LuOut,*) "tce_residual_tr2"
              call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c        write(LuOut,*) "tce_residual_tr3"
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
                goto 2335
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
              call tce_diis3(.false.,iter,.true.,.true.,.true.,.false.,
     1             d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2             size_tr1(axis),
     3             d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4             size_tr2(axis),
     5             d_rr3(axis),d_tr3(axis+dynaxis),k_tr3_offset(axis),
     6             size_tr3(axis),
     7             dummy,dummy,dummy,dummy,omega,1)
              call deletefile(d_rr3(axis))
              call deletefile(d_rr2(axis))
              call deletefile(d_rr1(axis))
              if (nodezero) call util_flush(LuOut)
            enddo ! iter loop
          call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 2335 continue
c
          endif ! respaxis(axis)
          enddo ! axis loop
c
c         IMAGINARY COMPONENT
c
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
c
          do omegasign=1,2
            if (omegasign.eq.1) then
              omega = (1.0d0)*ifreq(omegacount)
              dynaxis = 0
            elseif (omegasign.eq.2) then
              omega = (-1.0d0)*ifreq(omegacount)
              dynaxis = 3
            endif
            if (nodezero) write(LuOut,9431) omega
c
            call tce_diis_init()
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSDT-IR (imaginary component)"
c        write(LuOut,*) "singles initialization"
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
              call tce_zero(d_rr1(axis),size_tr1(axis))
c        write(LuOut,*) "doubles initialization"
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
              call tce_zero(d_rr2(axis),size_tr2(axis))
c        write(LuOut,*) "triples initialization"
              call tce_filename(rr3filename(axis),filename)
              call createfile(filename,d_rr3(axis),size_tr3(axis))
              call tce_zero(d_rr3(axis),size_tr3(axis))
!         write(LuOut,*) "eomccsdt_x1 rr1"
              call eomccsdt_x1(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
!         write(LuOut,*) "eomccsdt_x2 rr2"
              call eomccsdt_x2(d_f1,d_rr2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
!         write(LuOut,*) "eomccsdt_x3 rr3"
              call eomccsdt_x3(d_f1,d_rr3(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr3_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
!         write(LuOut,*) "daxpyfile rr1"
              call daxpyfile(1,(-1.0d0)*omega,d_tr1(axis+6),
     1             d_rr1(axis),size_tr1(axis))
!         write(LuOut,*) "daxpyfile rr2"
              call daxpyfile(1,(-1.0d0)*omega,d_tr2(axis+6),
     1             d_rr2(axis),size_tr2(axis))
!         write(LuOut,*) "daxpyfile rr3"
              call daxpyfile(1,(-1.0d0)*omega,d_tr3(axis+6),
     1             d_rr3(axis),size_tr3(axis))
!         write(LuOut,*) "reconcilefile rr1,rr2"
              call reconcilefile(d_rr1(axis),size_tr1(axis))
              call reconcilefile(d_rr2(axis),size_tr2(axis))
              call reconcilefile(d_rr3(axis),size_tr3(axis))
!         write(LuOut,*) "tce_residual_tr1"
              call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
!         write(LuOut,*) "tce_residual_tr2"
              call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
!         write(LuOut,*) "tce_residual_tr3"
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
                goto 2336
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
              call tce_diis3(.false.,iter,.true.,.true.,.true.,.false.,
     1             d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2             size_tr1(axis),
     3             d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4             size_tr2(axis),
     5             d_rr3(axis),d_tr3(axis+dynaxis),k_tr3_offset(axis),
     6             size_tr3(axis),
     7             dummy,dummy,dummy,dummy,0.0d0,1)
              call deletefile(d_rr3(axis))
              call deletefile(d_rr2(axis))
              call deletefile(d_rr1(axis))
              if (nodezero) call util_flush(LuOut)
            enddo ! iter loop
          call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 2336 continue
c
          enddo ! omegasign loop
c
          endif ! respaxis(axis)
          enddo ! axis loop
c
c CCSDT-IR evaluation step
c
          cpu=-util_cpusec()
          wall=-util_wallsec()
          do axisA = 1, 3
          do axisB = 1, axisA
          alpha(axisA,axisB)=0.0d0
          if (respaxis(axisA).and.respaxis(axisB)) then
            irrep_a=sym_abelian_axis(geom,axisA)
            irrep_b=sym_abelian_axis(geom,axisB)
            irrep_y=irrep_g
c
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,size_a0)
c
c            if (nodezero) write(LuOut,*) "axisA = ",axisA
c            if (nodezero) write(LuOut,*) "axisB = ",axisB
c
            alpha1=0.0d0
            alpha2=0.0d0
c
            irrep_d=irrep_a
            irrep_tr=irrep_b
c
c            write(LuOut,*) "alpha_1 AB"
            call ccsdt_lr_alpha1(d_a0,d_d1(axisA),d_t1,d_t2,d_t3,
     1           d_tr1(axisB+6),d_tr2(axisB+6),d_tr3(axisB+6),
     2           d_lambda1,d_lambda2,d_lambda3,k_a0_offset,
     3           k_d1_offset(axisA),k_t1_offset,k_t2_offset,
     4           k_t3_offset,k_tr1_offset(axisB),
     5           k_tr2_offset(axisB),k_tr3_offset(axisB),
     6           k_l1_offset,k_l2_offset,k_l3_offset)
c
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha1,1,0)
            call tce_zero(d_a0,size_a0)
c            write(LuOut,*) "alpha_1 AB = ",alpha1
c
            if (axisA.eq.axisB) then
              alpha2=alpha1
            else
              irrep_d=irrep_b
              irrep_tr=irrep_a
c              write(LuOut,*) "alpha_1 BA"
            call ccsdt_lr_alpha1(d_a0,d_d1(axisB),d_t1,d_t2,d_t3,
     1           d_tr1(axisA+6),d_tr2(axisA+6),d_tr3(axisA+6),
     2           d_lambda1,d_lambda2,d_lambda3,k_a0_offset,
     3           k_d1_offset(axisB),k_t1_offset,k_t2_offset,
     4           k_t3_offset,k_tr1_offset(axisA),
     5           k_tr2_offset(axisA),k_tr3_offset(axisA),
     6           k_l1_offset,k_l2_offset,k_l3_offset)
c
              call reconcilefile(d_a0,size_a0)
              call get_block(d_a0,alpha2,1,0)
              call tce_zero(d_a0,size_a0)
            endif ! axisA.eq.axisB
c            write(LuOut,*) "alpha_1 BA = ",alpha2
c
            alpha(axisA,axisB)=alpha(axisA,axisB)+alpha1+alpha2
c
            alpha1=0.0d0
            alpha2=0.0d0
            alpha3=0.0d0
c
            irrep_tra=irrep_a
            irrep_trb=irrep_b
c
c            write(LuOut,*) "alpha_2 A(+i) B(-i)"
            call ccsdt_lr_alpha2(d_f1,d_a0,d_t1,d_t2,d_t3,
     1           d_tr1(axisA),d_tr2(axisA),d_tr3(axisA),
     2           d_tr1(axisB+3),d_tr2(axisB+3),d_tr3(axisB+3),
     3           d_v2,d_lambda1,d_lambda2,d_lambda3,
     4           k_f1_offset,k_a0_offset,
     5           k_t1_offset,k_t2_offset,k_t3_offset,
     6           k_tr1_offset(axisA),k_tr2_offset(axisA),
     7           k_tr3_offset(axisA),k_tr1_offset(axisB),
     8           k_tr2_offset(axisB),k_tr3_offset(axisB),
     9           k_v2_offset,k_l1_offset,k_l2_offset,k_l3_offset)
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha1,1,0)
            call tce_zero(d_a0,size_a0)
c            write(LuOut,*) "alpha_2 A(+i) B(-i) = ",alpha1
c
c            write(LuOut,*) "alpha_2 A(-i) B(+i)"
            call ccsdt_lr_alpha2(d_f1,d_a0,d_t1,d_t2,d_t3,
     1           d_tr1(axisA+3),d_tr2(axisA+3),d_tr3(axisA+3),
     2           d_tr1(axisB),d_tr2(axisB),d_tr3(axisB),
     3           d_v2,d_lambda1,d_lambda2,d_lambda3,
     4           k_f1_offset,k_a0_offset,
     5           k_t1_offset,k_t2_offset,k_t3_offset,
     6           k_tr1_offset(axisA),k_tr2_offset(axisA),
     7           k_tr3_offset(axisA),k_tr1_offset(axisB),
     8           k_tr2_offset(axisB),k_tr3_offset(axisB),
     9           k_v2_offset,k_l1_offset,k_l2_offset,k_l3_offset)
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha2,1,0)
            call tce_zero(d_a0,size_a0)
c            write(LuOut,*) "alpha_2 A(-i) B(+i) = ",alpha2
c
c            write(LuOut,*) "alpha_2 A(r) B(r)"
            call ccsdt_lr_alpha2(d_f1,d_a0,d_t1,d_t2,d_t3,
     1           d_tr1(axisA+6),d_tr2(axisA+6),d_tr3(axisA+6),
     2           d_tr1(axisB+6),d_tr2(axisB+6),d_tr3(axisB+6),
     3           d_v2,d_lambda1,d_lambda2,d_lambda3,
     4           k_f1_offset,k_a0_offset,
     5           k_t1_offset,k_t2_offset,k_t3_offset,
     6           k_tr1_offset(axisA),k_tr2_offset(axisA),
     7           k_tr3_offset(axisA),k_tr1_offset(axisB),
     8           k_tr2_offset(axisB),k_tr3_offset(axisB),
     9           k_v2_offset,k_l1_offset,k_l2_offset,k_l3_offset)
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha3,1,0)
c            write(LuOut,*) "alpha_2 A(r) B(r) = ",alpha3
            call deletefile(d_a0)
c
            alpha(axisA,axisB)=alpha(axisA,axisB)-alpha1-alpha2+alpha3
            alpha(axisA,axisB)=-1.0d0*alpha(axisA,axisB)
c
c            write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
c            if (nodezero) write(LuOut,9020) cpu, wall
c
          endif ! respaxis(axis)
          enddo ! axisB loop
          enddo ! axisA loop
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
c
          do i = 1, 3
            do j = 1, i
              alphacopy(i,j)=alpha(i,j)
              alphacopy(j,i)=alpha(i,j)
            enddo
          enddo
c
          call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
c
          alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
          alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
          alphaani = dsqrt(0.5d0*alphaani)
c
          if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1        .and.nodezero) write(LuOut,9911)
c
          if (nodezero) write(LuOut,9434) "CCSDT Imaginary Response",
     1      ifreq(omegacount),auXnm/ifreq(omegacount),
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
c
          enddo ! omegacount loop
c
#endif
        endif ! imagresp
c
c -------------
c CCSDT-LR
c -------------
        if (lineresp) then
#if !defined(CCSDTLR)
      CALL ERRQUIT('CCSDTLR not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
c          if (nodezero) write(LuOut,*) "lineresp is go"
          do omegacount=1,anumfreq
c            if (nodezero) write(LuOut,*) "omegacount is go"
            omega = afreq(omegacount)
            if (omega.eq.(0.0d0)) then
              dynfreq = 1
            else
              dynfreq = 2
            endif
c
          do axis = 1, 3
c            if (nodezero) write(LuOut,*) "axis is go"
          if (respaxis(axis)) then
c            if (nodezero) write(LuOut,*) "respaxis(axis) is go"
            irrep_d=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            if (nodezero.and.util_print('mod1',print_default)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
            endif
            irrep_o=irrep_d
            irrep_x=irrep_d
c
          do omegasign=1,dynfreq
c            if (nodezero) write(LuOut,*) "omegasign is go"
            if (omegasign.eq.1) then
              omega = (1.0d0)*afreq(omegacount)
              dynaxis = 0
            elseif (omegasign.eq.2) then
              omega = (-1.0d0)*afreq(omegacount)
              dynaxis = 3
            endif
            if (nodezero) write(LuOut,9431) omega
c            if (nodezero) write(LuOut,*) "omega = ",omega
c            if (nodezero) write(LuOut,*) "dynaxis = ",dynaxis
c
            call tce_diis_init()
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSDT-LR"
c
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
              call tce_zero(d_rr1(axis),size_tr1(axis))
c        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega,d_tr1(axis+dynaxis),
     1             d_rr1(axis),size_tr1(axis))
c        write(LuOut,*) "eomccsdt_x1"
              call eomccsdt_x1(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
c        write(LuOut,*) "ccsdt_o1"
              call ccsdt_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     &             k_tr1_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset)
c
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
              call tce_zero(d_rr2(axis),size_tr2(axis))
c        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega,d_tr2(axis+dynaxis),
     1             d_rr2(axis),size_tr2(axis))
c        write(LuOut,*) "eomccsdt_x2"
              call eomccsdt_x2(d_f1,d_rr2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
c        write(LuOut,*) "ccsdt_o2"
              call ccsdt_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,d_t3,
     &             k_tr2_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset)
c
              call tce_filename(rr3filename(axis),filename)
              call createfile(filename,d_rr3(axis),size_tr3(axis))
              call tce_zero(d_rr3(axis),size_tr3(axis))
c        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega,d_tr3(axis+dynaxis),
     1             d_rr3(axis),size_tr3(axis))
c        write(LuOut,*) "eomccsdt_x3"
              call eomccsdt_x3(d_f1,d_rr3(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,k_tr3_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis))
c        write(LuOut,*) "ccsdt_o3"
              call ccsdt_o3(d_rr3(axis),d_d1(axis),d_t1,d_t2,d_t3,
     &             k_tr3_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset)
c
c        write(LuOut,*) "reconcilefile"
              call reconcilefile(d_rr1(axis),size_tr1(axis))
              call reconcilefile(d_rr2(axis),size_tr2(axis))
              call reconcilefile(d_rr3(axis),size_tr3(axis))
c        write(LuOut,*) "tce_residual_tr1"
              call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c        write(LuOut,*) "tce_residual_tr2"
              call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c        write(LuOut,*) "tce_residual_tr3"
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
c
          enddo ! omegasign loop
          endif ! respaxis(axis)
          enddo ! axis loop
c
c         Evaluate Dipole Polarizability
c
          do axisA = 1, 3
          do axisB = 1, axisA
          alpha(axisA,axisB)=0.0d0
          if (respaxis(axisA).and.respaxis(axisB)) then
            irrep_a=sym_abelian_axis(geom,axisA)
            irrep_b=sym_abelian_axis(geom,axisB)
            irrep_y=irrep_g
c
c            if (nodezero) write(LuOut,*) "axisA = ",axisA
c            if (nodezero) write(LuOut,*) "axisB = ",axisB
c
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,size_a0)
c
              alpha1=0.0d0
              alpha2=0.0d0
c
              if (omega.ne.0.0d0) then
                dynaxis = 3
              else
                dynaxis = 0
              endif
c
              irrep_c=irrep_a
              irrep_oa=irrep_a
              irrep_ob=irrep_b
c              if (nodezero) write(LuOut,*) "ccsdt_lr_alpha_offdiag #1"
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
c
              call reconcilefile(d_a0,size_a0)
              call get_block(d_a0,alpha1,1,0)
              call tce_zero(d_a0,size_a0)
c              write(LuOut,*) "alpha1 = ",alpha1
c
              if ((omega.eq.0.0d0).or.(axisA.eq.axisB)) then
                alpha2=alpha1
              else
c                if (nodezero) write(LuOut,*) "ccsdt_lr_alpha_offdiag #2"
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
c
                call reconcilefile(d_a0,size_a0)
                call get_block(d_a0,alpha2,1,0)
                call tce_zero(d_a0,size_a0)
              endif ! axisA.eq.axisB
c              write(LuOut,*) "alpha2 = ",alpha2
c
              alpha(axisA,axisB)=(-0.5d0)*(alpha1+alpha2)
c
            call deletefile(d_a0)
c
c            write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
c            if (nodezero) write(LuOut,9020) cpu, wall
c
          endif ! respaxis(axis)
          enddo ! axisB loop
          enddo ! axisA loop
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
c
          do i = 1, 3
            do j = 1, i
              alphacopy(i,j)=alpha(i,j)
              alphacopy(j,i)=alpha(i,j)
            enddo
          enddo
c
          call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
c
          alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
          alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
          alphaani = dsqrt(0.5d0*alphaani)
c
          if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1        .and.nodezero) write(LuOut,9911)
c
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
c
          enddo ! omegacount loop
c
#endif
        endif ! lineresp
c
c -------------
c CCSDT-YR
c -------------
        if (leftresp) then
#if !defined(CCSDTLR)
      CALL ERRQUIT('CCSDTLR not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
c
          do omegacount=1,anumfreq
            omega1 = afreq(omegacount)
            if (omega1.eq.(0.0d0)) then
              dynfreq = 1
            else
              dynfreq = 2
            endif
c
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
c
          do omegasign=1,dynfreq
            if (omegasign.eq.1) then
              omega1 = (1.0d0)*afreq(omegacount)
              dynaxis = 0
            elseif (omegasign.eq.2) then
              omega1 = (-1.0d0)*afreq(omegacount)
              dynaxis = 3
            endif
            if (nodezero) write(LuOut,9431) omega1
c            if (nodezero) write(LuOut,*) "omega1 = ",omega1
c            if (nodezero) write(LuOut,*) "dynaxis = ",dynaxis
c
            call tce_diis_init()
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSDT-YR (Static)"
c        write(LuOut,*) "tce_filename"
              call tce_filename(yrr1filename(axis),filename)
c        write(LuOut,*) "createfile"
              call createfile(filename,d_yrr1(axis),size_yr1(axis))
c        write(LuOut,*) "tce_zero"
              call tce_zero(d_yrr1(axis),size_yr1(axis))
              call daxpyfile(1,(1.0d0)*omega1,d_yr1(axis+dynaxis),
     1             d_yrr1(axis),size_yr1(axis))
c        write(LuOut,*) "ccsdt_y_o1"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_lambda_o1(d_yrr1(axis),d_d1(axis),d_t1,d_t2,
     1             d_t3,d_lambda1,d_lambda2,d_lambda3,
     2             k_yr1_offset(axis),k_d1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,
     3             k_l1_offset,k_l2_offset,k_l3_offset)
c              call ga_print(d_yrr1(axis))
c        write(LuOut,*) "ccsdt_y_tr1"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_y_tr1(d_f1,d_yrr1(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_lambda1,d_lambda2,d_lambda3,
     2             k_f1_offset,k_yr1_offset(axis),k_t1_offset,
     3             k_t2_offset,k_v2_offset,k_tr1_offset(axis),
     4             k_tr2_offset(axis),k_tr3_offset(axis),
     5             k_l1_offset,k_l2_offset,k_l3_offset)
c              call ga_print(d_yrr1(axis))
c        write(LuOut,*) "eomccsdt_y1"
              irrep_y=irrep_d ! using Lambda(1) amplitudes
              call eomccsdt_y1(d_f1,d_yrr1(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_yr1(axis+dynaxis),d_yr2(axis+dynaxis),
     2             d_yr3(axis+dynaxis),k_f1_offset,k_yr1_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_yr1_offset(axis),k_yr2_offset(axis),
     5             k_yr3_offset(axis))
c              call ga_print(d_yrr1(axis))
c        write(LuOut,*) "tce_filename"
              call tce_filename(yrr2filename(axis),filename)
c        write(LuOut,*) "createfile"
              call createfile(filename,d_yrr2(axis),size_yr2(axis))
c        write(LuOut,*) "tce_zero"
              call tce_zero(d_yrr2(axis),size_yr2(axis))
              call daxpyfile(1,(1.0d0)*omega1,d_yr2(axis+dynaxis),
     1             d_yrr2(axis),size_yr2(axis))
c        write(LuOut,*) "ccsdt_y_o2"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_lambda_o2(d_yrr2(axis),d_d1(axis),d_t1,d_t2,
     1             d_lambda1,d_lambda2,d_lambda3,k_yr2_offset(axis),
     2             k_d1_offset(axis),k_t1_offset,k_t2_offset,
     3             k_l1_offset,k_l2_offset,k_l3_offset)
c        write(LuOut,*) "ccsdt_y_tr2"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_y_tr2(d_f1,d_yrr2(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_lambda1,d_lambda2,
     2             d_lambda3,k_f1_offset,k_yr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_tr3_offset(axis),k_l1_offset,k_l2_offset,
     6             k_l3_offset)
c        write(LuOut,*) "eomccsdt_y2"
              irrep_y=irrep_d ! using Lambda(1) amplitudes
              call eomccsdt_y2(d_f1,d_yrr2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_yr1(axis+dynaxis),d_yr2(axis+dynaxis),
     2             d_yr3(axis+dynaxis),k_f1_offset,k_yr2_offset(axis),
     3             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     4             k_yr1_offset(axis),k_yr2_offset(axis),
     5             k_yr3_offset(axis))
c        write(LuOut,*) "tce_filename"
              call tce_filename(yrr3filename(axis),filename)
c        write(LuOut,*) "createfile"
              call createfile(filename,d_yrr3(axis),size_yr3(axis))
c        write(LuOut,*) "tce_zero"
              call tce_zero(d_yrr2(axis),size_yr2(axis))
              call daxpyfile(1,(1.0d0)*omega1,d_yr2(axis+dynaxis),
     1             d_yrr2(axis),size_yr2(axis))
c        write(LuOut,*) "ccsdt_y_o3"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_lambda_o3(d_yrr3(axis),d_d1(axis),d_t1,
     1             d_lambda2,d_lambda3,k_yr3_offset(axis),
     2             k_d1_offset(axis),k_t1_offset,
     3             k_l2_offset,k_l3_offset)
c        write(LuOut,*) "ccsdt_y_tr3"
              irrep_y=irrep_g ! using Lambda(0) amplitudes
              call ccsdt_y_tr3(d_f1,d_yrr3(axis),d_t1,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),d_lambda2,
     2             d_lambda3,k_f1_offset,k_yr3_offset(axis),
     3             k_t1_offset,k_v2_offset,
     4             k_tr1_offset(axis),k_tr2_offset(axis),
     5             k_l2_offset,k_l3_offset)
c        write(LuOut,*) "eomccsdt_y3"
              irrep_y=irrep_d ! using Lambda(1) amplitudes
              call eomccsdt_y3(d_f1,d_yrr3(axis),d_t1,d_t2,d_v2,
     1             d_yr1(axis+dynaxis),d_yr2(axis+dynaxis),
     2             d_yr3(axis+dynaxis),k_f1_offset,k_yr3_offset(axis),
     3             k_t1_offset,k_t2_offset,k_v2_offset,
     4             k_yr1_offset(axis),k_yr2_offset(axis),
     5             k_yr3_offset(axis))
c        write(LuOut,*) "reconcilefile"
              call reconcilefile(d_yrr1(axis),size_yr1(axis))
              call reconcilefile(d_yrr2(axis),size_yr2(axis))
              call reconcilefile(d_yrr3(axis),size_yr3(axis))
c        write(LuOut,*) "ddotfile"
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
c
c                call ga_print(d_yrr1(axis))
c                call ga_print(d_yrr2(axis))
c                call ga_print(d_yrr3(axis))
c
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
c       write(LuOut,*) "before tce_diis2"
              irrep_y=irrep_d
              call tce_diis2(.true.,iter,.true.,.true.,.false.,.false.,
     1             d_yrr1(axis),d_yr1(axis+dynaxis),k_yr1_offset(axis),
     2             size_yr1(axis),
     3             d_yrr2(axis),d_yr2(axis+dynaxis),k_yr2_offset(axis),
     4             size_yr2(axis),
     5             d_yrr3(axis),d_yr3(axis+dynaxis),k_yr3_offset(axis),
     6             size_yr3(axis),
     6             dummy,dummy,dummy,dummy,omega)
c       write(LuOut,*) "after tce_diis2"
              call deletefile(d_yrr3(axis))
              call deletefile(d_yrr2(axis))
              call deletefile(d_yrr1(axis))
              if (nodezero) call util_flush(LuOut)
            enddo ! iter loop
          call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 2339 continue
c
          enddo ! omegasign loop
c
          endif ! respaxis(axis)
          enddo ! axis loop
          call errquit('tce_energy: manual stop',911,CALC_ERR)
c
c         EVALUATE POLARIZABILITY
c
          cpu=-util_cpusec()
          wall=-util_wallsec()
          do axisA = 1, 3
          irrep_a=sym_abelian_axis(geom,axisA)
          do axisB = 1, 3
          irrep_b=sym_abelian_axis(geom,axisB)
          alpha(axisA,axisB)=0.0d0
          if ( respaxis(axisA).and.respaxis(axisB) ) then
c
            alpha1=0.0d0
            alpha2=0.0d0
            alpha3=0.0d0
            alpha4=0.0d0
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,size_a0)
c
            irrep_y=irrep_g
            irrep_d=irrep_a
            irrep_tr=irrep_b
c
            call alpha_1(d_d1(axisA),d_a0,d_t1,d_t2,
     1           d_tr1(axisB),d_tr2(axisB),
     2           d_lambda1,d_lambda2,k_d1_offset(axisA),k_a0_offset,
     4           k_t1_offset,k_t2_offset,k_tr1_offset(axisB),
     5           k_tr2_offset(axisB),k_l1_offset,k_l2_offset)
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha1,1,0)
            call tce_zero(d_a0,size_a0)
c
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
c
            irrep_o=irrep_a
            irrep_y=irrep_b
c
            call alpha_2_yr(d_a0,d_d1(axisA),d_t1,d_t2,
     1           d_yr1(axisB),d_yr2(axisB),
     2           k_a0_offset,k_d1_offset(axisA),k_t1_offset,
     3           k_t2_offset,k_yr1_offset(axisB),k_yr2_offset(axisB))
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha3,1,0)
            call tce_zero(d_a0,size_a0)
c
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
c
c            write(LuOut,*) "alpha1 = ",alpha1
c            write(LuOut,*) "alpha2 = ",alpha2
c            write(LuOut,*) "alpha3 = ",alpha3
c            write(LuOut,*) "alpha4 = ",alpha4
            alpha(axisA,axisB)=(-0.5d0)*(alpha1+alpha2+alpha3+alpha4)
            call deletefile(d_a0)
c
          endif ! respaxis(axisA,axisB)
          enddo ! axisB loop
          enddo ! axisA loop
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
c
          do j = 1, 3
            do i = 1, 3
              alphacopy(i,j)=alpha(i,j)
            enddo
          enddo
c
          call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
c
          alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
          alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
          alphaani = dsqrt(0.5d0*alphaani)
c
          if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1        .and.nodezero) write(LuOut,9911)
c
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
c
c
c         EVALUATE HYPERPOLARIZABILITY
c
          cpu=-util_cpusec()
          wall=-util_wallsec()
          do axisA = 1, 3
          do axisB = 1, 3
          do axisC = 1, 3
            beta(axisA,axisB,axisC)=0.0d0
          enddo
          enddo
          enddo
c
          do axisA = 1, 3
            irrep_a=sym_abelian_axis(geom,axisA)
          do axisB = 1, 3
            irrep_b=sym_abelian_axis(geom,axisB)
          do axisC = 1, 3
            irrep_c=sym_abelian_axis(geom,axisC)
          if ( ( respaxis(axisA).and.respaxis(axisB) )
     &                          .and.respaxis(axisC) ) then
c
c            write(LuOut,*) "axisA = ",axisA
c            write(LuOut,*) "axisB = ",axisB
c            write(LuOut,*) "axisC = ",axisC
c            call util_flush(LuOut)
c
            beta1=0.0d0
            beta2=0.0d0
            beta3=0.0d0
            beta4=0.0d0
            call tce_filename('b0',filename)
            call createfile(filename,d_b0,size_b0)
c
            irrep_y=irrep_g
            irrep_tra=irrep_a
            irrep_trb=irrep_b
            irrep_trc=irrep_c
c            write(LuOut,*) "beta_1"
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
c
            irrep_y=irrep_g
            irrep_d=irrep_a
            irrep_tra=irrep_b
            irrep_trb=irrep_c
c            write(LuOut,*) "beta_2"
            call beta_2(d_d1(axisA),d_b0,d_tr1(axisB),d_tr2(axisB),
     1           d_tr1(axisC),d_tr2(axisC),d_lambda1,d_lambda2,
     2           k_d1_offset(axisA),k_b0_offset,
     3           k_tr1_offset(axisB),k_tr2_offset(axisB),
     4           k_tr1_offset(axisC),k_tr2_offset(axisC),
     5           k_l1_offset,k_l2_offset)
            call reconcilefile(d_b0,size_b0)
            call get_block(d_b0,beta2,1,0)
            call tce_zero(d_b0,size_b0)
c
            irrep_y=irrep_a
            irrep_tra=irrep_b
            irrep_trb=irrep_c
c            write(LuOut,*) "beta_3"
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
c
            irrep_y=irrep_a
            irrep_d=irrep_b
            irrep_tr=irrep_c
c            write(LuOut,*) "beta_4"
            call beta_4(d_d1(axisB),d_b0,d_t1,d_t2,
     1           d_tr1(axisC),d_tr2(axisC),d_yr1(axisA),d_yr2(axisA),
     2           k_d1_offset(axisB),k_b0_offset,k_t1_offset,k_t2_offset,
     3           k_tr1_offset(axisC),k_tr2_offset(axisC),
     4           k_yr1_offset(axisA),k_yr2_offset(axisA))
            call reconcilefile(d_b0,size_b0)
            call get_block(d_b0,beta4,1,0)
            call tce_zero(d_b0,size_b0)
c
c            write(LuOut,*) "beta1 = ",beta1
c            write(LuOut,*) "beta2 = ",beta2
c            write(LuOut,*) "beta3 = ",beta3
c            write(LuOut,*) "beta4 = ",beta4
c
c           This term contributes to only one component
c
            beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta1
c
c           These terms both contribute to three components
c
            beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta2
            beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta2
            beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta2
c
            beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta3
            beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta3
            beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta3
c
c           This term contributes to six components
c
            beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta4
            beta(axisA,axisC,axisB)=beta(axisA,axisC,axisB)-beta4
            beta(axisB,axisA,axisC)=beta(axisB,axisA,axisC)-beta4
            beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta4
            beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta4
            beta(axisC,axisB,axisA)=beta(axisC,axisB,axisA)-beta4
c
            call deletefile(d_b0)
c
          endif ! respaxis(axisA,axisB,axisC)
          enddo ! axisB loop
          enddo ! axisB loop
          enddo ! axisA loop
c
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
c
          enddo ! omegacount loop
c
#endif
        endif ! leftresp
c
      endif
      if (model.eq.'ccsdta') then
c -----------------------------
c CCSDt (CCSDT-active)   right
c -----------------------------
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
cedo            goto 330
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
      endif
      if (model.eq.'ccsdtq') then
#if !defined(CCSDTQ)
      CALL ERRQUIT('CCSDTQ not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
c -------------
c CCSDTQ  right
c -------------
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
c -------------
c CCSDTQ  left
c -------------
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
c            goto 100
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
c -------------
c CCSDTQ-LR
c -------------
 2437   continue
        if (lineresp) then
c
c          if (nodezero) write(LuOut,*) "lineresp is go"
          do omegacount=1,anumfreq
c            if (nodezero) write(LuOut,*) "omegacount is go"
            omega1 = afreq(omegacount)
            if (omega1.eq.(0.0d0)) then
              dynfreq = 1
            else
              dynfreq = 2
            endif
c
          do axis = 1, 3
c            if (nodezero) write(LuOut,*) "axis is go"
          if (respaxis(axis)) then
c            if (nodezero) write(LuOut,*) "respaxis(axis) is go"
            irrep_d=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            if (nodezero.and.util_print('mod1',print_default)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
            endif
            irrep_o=irrep_d
            irrep_x=irrep_d
c
          do omegasign=1,dynfreq
c            if (nodezero) write(LuOut,*) "omegasign is go"
            if (omegasign.eq.1) then
              omega1 = (1.0d0)*afreq(omegacount)
              dynaxis = 0
            elseif (omegasign.eq.2) then
              omega1 = (-1.0d0)*afreq(omegacount)
              dynaxis = 3
            endif
            if (nodezero) write(LuOut,9431) omega1
c            if (nodezero) write(LuOut,*) "omega1 = ",omega1
c            if (nodezero) write(LuOut,*) "dynaxis = ",dynaxis
c
            call tce_diis_init()
            do iter=1,maxiter
              cpu=-util_cpusec()
              wall=-util_wallsec()
              if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSDTQ-LR"
c
              call tce_filename(rr1filename(axis),filename)
              call createfile(filename,d_rr1(axis),size_tr1(axis))
              call tce_zero(d_rr1(axis),size_tr1(axis))
c        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega1,d_tr1(axis+dynaxis),
     1             d_rr1(axis),size_tr1(axis))
c        write(LuOut,*) "eomccsdtq_x1"
              call eomccsdtq_x1(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,
     3             k_tr1_offset(axis),
     4             k_t1_offset,k_t2_offset,k_v2_offset,
     5             k_tr1_offset(axis),k_tr2_offset(axis),
     6             k_tr3_offset(axis))
c        write(LuOut,*) "ccsdtq1_o1"
              call ccsdtq_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     &             k_tr1_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset)
c
              call tce_filename(rr2filename(axis),filename)
              call createfile(filename,d_rr2(axis),size_tr2(axis))
              call tce_zero(d_rr2(axis),size_tr2(axis))
c        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega1,d_tr2(axis+dynaxis),
     1             d_rr2(axis),size_tr2(axis))
c        write(LuOut,*) "eomccsdtq_x2"
              call eomccsdtq_x2(d_f1,d_rr2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_tr4(axis+dynaxis),
     3             k_f1_offset,k_tr2_offset(axis),
     4             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     5             k_tr1_offset(axis),k_tr2_offset(axis),
     6             k_tr3_offset(axis),k_tr4_offset(axis))
c        write(LuOut,*) "ccsdtq_o2"
              call ccsdtq_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,d_t3,
     &             k_tr2_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset)
c
              call tce_filename(rr3filename(axis),filename)
              call createfile(filename,d_rr3(axis),size_tr3(axis))
              call tce_zero(d_rr3(axis),size_tr3(axis))
c        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega1,d_tr3(axis+dynaxis),
     1             d_rr3(axis),size_tr3(axis))
c        write(LuOut,*) "eomccsdtq_x3"
              call eomccsdtq_x3(d_f1,d_rr3(axis),d_t1,d_t2,d_t3,d_t4,
     1             d_v2,d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_tr4(axis+dynaxis),
     3             k_f1_offset,k_tr3_offset(axis),
     4             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5             k_v2_offset,
     6             k_tr1_offset(axis),k_tr2_offset(axis),
     7             k_tr3_offset(axis),k_tr4_offset(axis))
c        write(LuOut,*) "ccsdtq_o3"
              call ccsdtq_o3(d_rr3(axis),d_d1(axis),d_t1,d_t2,d_t3,d_t4,
     &             k_tr3_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset)
c
              call tce_filename(rr4filename(axis),filename)
              call createfile(filename,d_rr4(axis),size_tr4(axis))
              call tce_zero(d_rr4(axis),size_tr4(axis))
c        write(LuOut,*) "daxpyfile"
              call daxpyfile(1,(-1.0d0)*omega1,d_tr4(axis+dynaxis),
     1             d_rr4(axis),size_tr4(axis))
c        write(LuOut,*) "eomccsdtq_x4"
              call eomccsdtq_x4(d_f1,d_rr4(axis),d_t1,d_t2,d_t3,d_t4,
     1             d_v2,d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_tr4(axis+dynaxis),
     3             k_f1_offset,k_tr4_offset(axis),
     4             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5             k_v2_offset,
     6             k_tr1_offset(axis),k_tr2_offset(axis),
     7             k_tr3_offset(axis),k_tr4_offset(axis))
c        write(LuOut,*) "ccsdtq_o4"
              call ccsdtq_o4(d_rr4(axis),d_d1(axis),d_t1,d_t2,d_t3,d_t4,
     &             k_tr4_offset(axis),k_d1_offset(axis),
     &             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset)
c
c        write(LuOut,*) "reconcilefile"
              call reconcilefile(d_rr1(axis),size_tr1(axis))
              call reconcilefile(d_rr2(axis),size_tr2(axis))
              call reconcilefile(d_rr3(axis),size_tr3(axis))
              call reconcilefile(d_rr4(axis),size_tr4(axis))
c        write(LuOut,*) "tce_residual_tr1"
              call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c        write(LuOut,*) "tce_residual_tr2"
              call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c        write(LuOut,*) "tce_residual_tr3"
              call tce_residual_tr3(d_rr3(axis),k_tr3_offset(axis),rr3)
c        write(LuOut,*) "tce_residual_tr4"
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
c
          enddo ! omegasign loop
          endif ! respaxis(axis)
          enddo ! axis loop
c
c         Evaluate Dipole Polarizability
c
          do axisA = 1, 3
          do axisB = 1, axisA
          alpha(axisA,axisB)=0.0d0
          if (respaxis(axisA).and.respaxis(axisB)) then
            irrep_a=sym_abelian_axis(geom,axisA)
            irrep_b=sym_abelian_axis(geom,axisB)
            irrep_y=irrep_g
c
c            if (nodezero) write(LuOut,*) "axisA = ",axisA
c            if (nodezero) write(LuOut,*) "axisB = ",axisB
c
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,size_a0)
c
              alpha1=0.0d0
              alpha2=0.0d0
c
              if (omega1.ne.0.0d0) then
                dynaxis = 3
              else
                dynaxis = 0
              endif
c
              irrep_c=irrep_a
              irrep_oa=irrep_a
              irrep_ob=irrep_b
              irrep_tra=irrep_a
              irrep_trb=irrep_b
c              if (nodezero) write(LuOut,*) "ccsdtq_lr_alpha"
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
c
              call reconcilefile(d_a0,size_a0)
              call get_block(d_a0,alpha1,1,0)
              call tce_zero(d_a0,size_a0)
c              write(LuOut,*) "alpha1 = ",alpha1
c
              if (omega1.eq.0.0d0) then
                alpha2=alpha1
              else
c                if (nodezero) write(LuOut,*) "ccsdtq_lr_alpha"
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
c
                call reconcilefile(d_a0,size_a0)
                call get_block(d_a0,alpha2,1,0)
                call tce_zero(d_a0,size_a0)
              endif ! axisA.eq.axisB
c              write(LuOut,*) "alpha2 = ",alpha2
c
              alpha(axisA,axisB)=(-0.5d0)*(alpha1+alpha2)
c
            call deletefile(d_a0)
c
c            write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
c            if (nodezero) write(LuOut,9020) cpu, wall
c
          endif ! respaxis(axis)
          enddo ! axisB loop
          enddo ! axisA loop
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
c
          do i = 1, 3
            do j = 1, i
              alphacopy(i,j)=alpha(i,j)
              alphacopy(j,i)=alpha(i,j)
            enddo
          enddo
c
          call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
c
          alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
          alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
          alphaani = dsqrt(0.5d0*alphaani)
c
          if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1        .and.nodezero) write(LuOut,9911)
c
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
c
          enddo ! omegacount loop
c
        endif ! lineresp
c
#endif
      endif
      if ((model.eq.'qcisd').or.do_qcisd) then
c -------------
c QCISD   right
c -------------
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
c            goto 100
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
c -------------------
c QCISD(T) and QCISD[T]
c -------------------
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
c -------------------------
c CR-QCISD(T)
c -------------------------
        if (do_cr_pt.or.(perturbative.eq.'cr_(t)')) then
          cpu = - util_cpusec()
          wall = - util_wallsec()
c          write(LuOut,*) 'tce_energy: cr_qcisd_t'
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
c -------------
c CISD    right
c -------------
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
c -------------
c CISDT   right
c -------------
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
c -------------
c CISDTQ  right
c -------------
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
c -------------
c MBPT(2) right
c -------------
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
c -------------
c MBPT(3) right
c -------------
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
c->d3p975
            call mbpt2_t2(d_f1,d_r2,d_u1,d_u2,d_w2,d_w3,d_v2,
     1                    k_f1_offset,k_r2_offset,k_u1_offset,
     2                    k_u2_offset,k_w2_offset,
     3                    k_w3_offset,k_v2_offset)
c<-d3p975
            call reconcilefile(d_r2,size_r2)
            call tce_filename('r3',filename)
            call createfile(filename,d_r3,size_r3)
c->d3p975
            call mbpt2_t3(d_f1,d_r3,d_u1,d_u2,d_w2,d_w3,d_v2,
     1                    k_f1_offset,k_r3_offset,k_u1_offset,
     2                    k_u2_offset,k_w2_offset,k_w3_offset,
     3                    k_v2_offset)
c<-d3p975
            call reconcilefile(d_r3,size_r3)
c->d3p975
          elseif ((model.eq.'mbpt4sdq').or.
     1             (model.eq.'mbpt4sdq_t')) then
            call tce_filename('r2',filename)
            call createfile(filename,d_r2,size_r2)
            call mbpt2_t2_sd(d_f1,d_r2,d_u1,d_u2,d_w2,d_v2,
     1                       k_f1_offset,k_r2_offset,k_u1_offset,
     1                       k_u2_offset,k_w2_offset,k_v2_offset)
c
            call reconcilefile(d_r2,size_r2)
c<-d3p975
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
c->d3p975
            call tce_diis(.false.,iter,.true.,.true.,.true.,.false.,
     1                    d_r1,d_w1,k_w1_offset,size_w1,
     2                    d_r2,d_w2,k_w2_offset,size_w2,
     3                    d_r3,d_w3,k_w3_offset,size_w3,
     4                    dummy,dummy,dummy,dummy)
c<-d3p975
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
c -------------
c MBPT(4) right
c -------------
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
c
            if (model.eq.'mbpt4sdq_t') then
              cpu=-util_cpusec()
              wall=-util_wallsec()
              call ccsd_t(d_u1,k_u1_offset,d_u2,k_u2_offset,
     1                    d_v2,k_v2_offset,pt3_1,pt3_2,size_t1)
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
            endif
c
            if (nodezero) then
              if(model.eq.'mbpt4') then
                write(LuOut,9060)
                write(LuOut,9070) "MBPT(4)",mbpt4
                write(LuOut,9080) "MBPT(4)",ref + corr
              endif
c
              if (model.eq.'mbpt4sdq') then
                write(LuOut,9060)
                write(LuOut,9070) "MBPT(4)(SDQ)",mbpt4
                write(LuOut,9080) "MBPT(4)(SDQ)",ref + corr
              endif
c
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
c            goto 2010
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
c
      if (nodezero) call util_flush(LuOut)
c
c     Delete files for energy calculation
c
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
     1    call errquit("tce_energy: MA problem",27,MA_ERR)
      endif
      if (needt3) then
        if (.not.ma_pop_stack(l_r3_offset))
     1    call errquit("tce_energy: MA problem",21,MA_ERR)
      endif
      if (needt3a) then
        if (.not.ma_pop_stack(l_r3_offset))
     1    call errquit("tce_energy: MA problem",21,MA_ERR)
      endif
      if (needt2) then
        if (.not.ma_pop_stack(l_r2_offset))
     1    call errquit("tce_energy: MA problem",22,MA_ERR)
      endif
      if (needt1) then
        if (.not.ma_pop_stack(l_r1_offset))
     1    call errquit("tce_energy: MA problem",23,MA_ERR)
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
      if (lineresp.or.imagresp) then
        if (.not.ma_pop_stack(l_a0_offset))
     1    call errquit("tce_energy: MA problem",26,MA_ERR)
      endif
      if (left) then
        if (multipole.ge.3) then
          if (.not.ma_pop_stack(l_o0_offset))
     1      call errquit("tce_energy: MA problem",300,MA_ERR)
        endif
        if (multipole.ge.2) then
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
c link to QM/MM
        if (.not.rtdb_put(rtdb,'tce:enegr',mt_dbl,1,ref+corr))
     1  call errquit('tce_energy: RTDB problem',0,MA_ERR)
c
        tce_energy=.true.
      else
c
c     ============================
c     Excited-state CI, CC, & MBPT
c     ============================
c
        irrep_g = 0
c
c       Restricted closed shell: irrep_g is always zero.
c
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
c
          if(nodezero) then
           write(LuOut,9251) hbard
          end if
c new EOMCC solver
       if(.not.((model.ne.'ccsd').and.(eom_solver.eq.2))) then
       if (.not.ma_push_get(mt_dbl,hbard*hbard,'hbar',
     1  l_hbar,k_hbar)) call errquit('tce_eom_xdiagon: MA problem',0,
     2  MA_ERR)
        do ihbar=1,hbard*hbard
         dbl_mb(k_hbar+ihbar-1)=0.0d0
        enddo
        do ihbar=1,hbard
         dbl_mb(k_hbar+(ihbar-1)*hbard+ihbar-1)=1.0d+8
        enddo
        end if
c
c
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
c
          if (needt1) then
            call tce_x1_offset(l_x1_offset,k_x1_offset,size_x1)
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
            call tce_x2_offset(l_x2_offset,k_x2_offset,size_x2)
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
cccx            if (left) then
cccx              call tce_y3_offset(l_y3_offset,k_y3_offset,size_y3)
cccx              call tce_filename('ry3',filename)
cccx              call createfile(filename,d_ry3,size_y3)
cccx              call gatoeaf(d_ry3)
cccx            endif
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
c
c         ------------------------------
c         Generate initial trial vectors
c         ------------------------------
c
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (model.eq.'ccsd') then


      if(eom_solver.eq.1) then !  ------- EOM-SOLVER ----- (solver 1)
c
c
      if(nodezero) then
       write(LuOut,*)'EOMCCSD SOLVER TYPE 1'
       call util_flush(LuOut)
      end if
c
      if(ccsd_var.ne.'cc2') then
      ister=0
c rec_s_1
       CALL OFFSET_eomccsd_x1_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs1',-1,MA_ERR)
c rec_s_2
       CALL OFFSET_eomccsd_x1_2_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_2_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs2',-1,MA_ERR)
c rec_s_3
       CALL OFFSET_eomccsd_x1_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs3',-1,MA_ERR)
c rec_s_4
       CALL OFFSET_eomccsd_x1_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs4',-1,MA_ERR)
c rec_s_5
       CALL OFFSET_eomccsd_x1_7_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_7_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_s_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs5',-1,MA_ERR)
c rec_d_1
       CALL OFFSET_eomccsd_x2_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd1',-1,MA_ERR)
c rec_d_2
       CALL OFFSET_eomccsd_x2_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_3_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd2',-1,MA_ERR)
c rec_d_3
       CALL OFFSET_eomccsd_x2_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd3',-1,MA_ERR)
c rec_d_4
       CALL OFFSET_eomccsd_x2_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd4',-1,MA_ERR)
c rec_d_5
       CALL OFFSET_eomccsd_x2_6_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_6_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd5',-1,MA_ERR)
c rec_d_6
       CALL OFFSET_eomccsd_x2_8_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_6,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd6',-1,MA_ERR)
c rec_d_7
       CALL OFFSET_eomccsd_x2_8_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_3_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_7,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd7',-1,MA_ERR)
c rec_d_8
       CALL OFFSET_eomccsd_x2_8_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_4_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_8,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd8',-1,MA_ERR)
c rec_d_9
      CALL OFFSET_eomccsd_x2_10_1_1(l_aux_off,k_aux_off,size_aux)
      CALL TCE_FILENAME('eomccsd_x2_10_1_1_i2',filename)
      CALL CREATEFILE(filename,rec_d_9,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd9',-1,MA_ERR)
      end if 
c
      if(ccsd_var.eq.'lr-ccsd') then
            call tce_x1_offset(l_j1_offset,k_j1_offset,size_j1)
            call tce_filename('j1',filename)
            call createfile(filename,d_j1,size_j1)
cc            call gatoeaf(d_j1)
            call tce_zero(d_j1,size_j1)
c
            call tce_x2_offset(l_j2_offset,k_j2_offset,size_j2)
            call tce_filename('j2',filename)
            call createfile(filename,d_j2,size_j2)
cc            call gatoeaf(d_j2)
            call tce_zero(d_j2,size_j2)
c
            call tce_loc_j1_one(d_j1,k_j1_offset)
            call tce_loc_j2_one(d_j2,k_j2_offset)
c
            call tce_loc_j1(d_j1,k_j1_offset,d_t1,k_t1_offset)
            call tce_loc_j2(d_j2,k_j2_offset,d_t2,k_t2_offset)
c
            call tce_loc_j1_inv(d_j1,k_j1_offset)
            call tce_loc_j2_inv(d_j2,k_j2_offset)
      end if
c
            call tce_eom_xguess(needt1,needt2,needt3,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
            if (nxtrials.eq.0) goto 200
ccxsg
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
ccxsg
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
c
c           call hbar_init
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
c - T1/X1 LOCALIZATION ----------
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,size_x1,'x1_local',
     1      l_x1_local,k_x1_local))
     1      call errquit('x1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
        call ma_zero(dbl_mb(k_x1_local),size_x1)
c    copy d_t1 ==> l_t1_local
c    copy x1(ivec) ==> l_x1_local
cc        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
cc        call ga_get(x1(ivec),1,size_x1,1,1,dbl_mb(k_x1_local),1)
      call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
      call get_block(x1(ivec),dbl_mb(k_x1_local),size_x1,0)
ccx       end if
c -------------------------------
         if(ccsd_var.eq.'cc2') then
          call cc2_x1(d_f1,xp1(ivec),k_t1_local,d_t2,d_v2,
     &      k_x1_local,x2(ivec),
     &      k_f1_offset,k_x1_offset,k_t1_offset,k_t2_offset,
     &      k_v2_offset,k_x1_offset,k_x2_offset)
         else
                  call eomccsd_x1_n(d_f1,xp1(ivec),k_t1_local,d_t2,d_v2,
     1              k_x1_local,x2(ivec),k_f1_offset,k_x1_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,
     3              k_x1_offset,k_x2_offset,
     4 ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9)
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
                  call eomccsd_x2_n(d_f1,xp2(ivec),k_t1_local,d_t2,d_v2,
     1              k_x1_local,x2(ivec),k_f1_offset,k_x2_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,
     3              k_x1_offset,k_x2_offset,size_x1,size_x2,
     4     ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5             rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6             rec_d_6,rec_d_7,rec_d_8,rec_d_9)
c
           ister=1
c
         end if
                  call reconcilefile(xp2(ivec),size_x2)
c - T1/X1 LOCALIZATION --
         if(.not.MA_POP_STACK(l_x1_local))
     &      call errquit('l_x1_local',4,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c -----------------------
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
              call tce_eom_xdiagon(needt1,needt2,needt3,needt4,
     1          size_x1,size_x2,size_x3,size_x4,
     2          k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
     3          d_rx1,d_rx2,d_rx3,d_rx4,
     4          dbl_mb(k_omegax),dbl_mb(k_residual),k_hbar,iter)
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
c               call hbar_tidy
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
c ---
            if(nodezero) then
             write(6,*)'largest EOMCCSD amplitudes: R1 and R2'
             call util_flush(6)
            end if
             do jvec=1,nroots_reduced
              call tce_print_x1(xc1(jvec),k_x1_offset,0.10d0,irrep_x)
              call tce_print_x2(xc2(jvec),k_x2_offset,0.10d0,irrep_x)
              call util_flush(LuOut)
             enddo
c ---
                goto 200
              endif
            enddo
c
           if(ccsd_var.eq.'lr-ccsd')then
            call deletefile(d_j2)
            if (.not.ma_pop_stack(l_j2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
            call deletefile(d_j1)
            if (.not.ma_pop_stack(l_j1_offset))
     1        call errquit("tce_energy: MA problem",37,MA_ERR)
           end if
c
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
c
         call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
c
c
c
      else !  ------- EOM-SOLVER --------------- (solver 2)
c
c
      if(nodezero) then
       write(LuOut,*)'EOMCCSD SOLVER TYPE 2'
       call util_flush(LuOut)
      end if
c
       modelname = "CIS-like initial guess for EOM-CCSD"
        if (nodezero) write(LuOut,9220)
     1  modelname(1:inp_strlen(modelname))
c
        call cis(LuOut,size_x1,k_x1_offset,d_f1,k_f1_offset,d_v2,
     1           k_v2_offset,k_omegax,l_omegax,d_rx1)
c
c
            if(nodezero) then
             write(LuOut,*)'INITIAL STARTS -----'
             call util_flush(6)
            end if
            do iroot=1,nroots_reduced
             call tce_print_x1(x1_ini(iroot),k_x1_offset,0.10d0,irrep_x)
            enddo
            if(nodezero) then
             write(LuOut,*)'                    '
             write(LuOut,*)'INITIAL STARTS -----'
             call util_flush(6)
            end if
c
c
            if (nxtrials.eq.0) goto 200
            call tce_eom_init
c
      if(ccsd_var.ne.'cc2') then
      ister=0
c rec_s_1
       CALL OFFSET_eomccsd_x1_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs1',-1,MA_ERR)
c rec_s_2
       CALL OFFSET_eomccsd_x1_2_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_2_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs2',-1,MA_ERR)
c rec_s_3
       CALL OFFSET_eomccsd_x1_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs3',-1,MA_ERR)
c rec_s_4
       CALL OFFSET_eomccsd_x1_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_s_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs4',-1,MA_ERR)
c rec_s_5
       CALL OFFSET_eomccsd_x1_7_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x1_7_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_s_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rs5',-1,MA_ERR)
c rec_d_1
       CALL OFFSET_eomccsd_x2_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_1_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_1,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd1',-1,MA_ERR)
c rec_d_2
       CALL OFFSET_eomccsd_x2_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_3_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_2,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd2',-1,MA_ERR)
c rec_d_3
       CALL OFFSET_eomccsd_x2_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_4_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_3,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd3',-1,MA_ERR)
c rec_d_4
       CALL OFFSET_eomccsd_x2_5_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_5_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_4,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd4',-1,MA_ERR)
c rec_d_5
       CALL OFFSET_eomccsd_x2_6_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_6_1_i1',filename)
       CALL CREATEFILE(filename,rec_d_5,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd5',-1,MA_ERR)
c rec_d_6
       CALL OFFSET_eomccsd_x2_8_1_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_1_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_6,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd6',-1,MA_ERR)
c rec_d_7
       CALL OFFSET_eomccsd_x2_8_3_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_3_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_7,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd7',-1,MA_ERR)
c rec_d_8
       CALL OFFSET_eomccsd_x2_8_4_1(l_aux_off,k_aux_off,size_aux)
       CALL TCE_FILENAME('eomccsd_x2_8_4_1_i2',filename)
       CALL CREATEFILE(filename,rec_d_8,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd8',-1,MA_ERR)
c rec_d_9
      CALL OFFSET_eomccsd_x2_10_1_1(l_aux_off,k_aux_off,size_aux)
      CALL TCE_FILENAME('eomccsd_x2_10_1_1_i2',filename)
      CALL CREATEFILE(filename,rec_d_9,size_aux)
       IF (.not.MA_POP_STACK(l_aux_off)) CALL ERRQUIT('rd9',-1,MA_ERR)
      end if 
c
      if(ccsd_var.eq.'lr-ccsd') then
            call tce_x1_offset(l_j1_offset,k_j1_offset,size_j1)
            call tce_filename('j1',filename)
            call createfile(filename,d_j1,size_j1)
cc            call gatoeaf(d_j1)
            call tce_zero(d_j1,size_j1)
c
            call tce_x2_offset(l_j2_offset,k_j2_offset,size_j2)
            call tce_filename('j2',filename)
            call createfile(filename,d_j2,size_j2)
cc            call gatoeaf(d_j2)
            call tce_zero(d_j2,size_j2)
c
            call tce_loc_j1_one(d_j1,k_j1_offset)
            call tce_loc_j2_one(d_j2,k_j2_offset)
c
            call tce_loc_j1(d_j1,k_j1_offset,d_t1,k_t1_offset)
            call tce_loc_j2(d_j2,k_j2_offset,d_t2,k_t2_offset)
c
            call tce_loc_j1_inv(d_j1,k_j1_offset)
            call tce_loc_j2_inv(d_j2,k_j2_offset)
      end if
c
c
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
c
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
c
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
c
c open k_h_eff  and k_c_eff
c
        if (.not.MA_PUSH_GET(mt_dbl,
     1      xdiis*xdiis*nroots_reduced*nroots_reduced,'h_eff',
     1      l_h_eff,k_h_eff)) call errquit('h_eff',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,xdiis*nroots_reduced*nroots_reduced,
     1      'c_eff',
     1      l_c_eff,k_c_eff)) call errquit('c_eff',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,nroots_reduced,'residual',
     1          l_residual,k_residual)) call errquit('l_res',1,MA_ERR)
c
            iter=1
c
             if(iter.eq.1) then
              do iroot=1,nroots_reduced
               call copyfile(x1_ini(iroot),xc1(iroot),size_x1)
               call tce_zero(xc2(iroot),size_x2)
              enddo
             endif 
c
 111         continue  !iter continuation
c
             imicro = 1
c
c
            do iroot=1,nroots_reduced
             call copyfile(xc1(iroot),x1(iroot),size_x1)
             call copyfile(xc2(iroot),x2(iroot),size_x2)
            enddo
c
            call norm_x_vector_mr(k_x1_offset,k_x2_offset,
     1                            size_x1,size_x2,nroots_reduced)
c
            do iroot=1,nroots_reduced
             call tce_zero(xp1(iroot),size_x1)
             call tce_zero(xp2(iroot),size_x2)
            enddo
c
             imicro = 1
c
             do i=1,xdiis*nroots_reduced*nroots_reduced
              dbl_mb(k_c_eff+i-1)=0.0d0
             enddo
             do i=1,nroots_reduced
             dbl_mb(k_c_eff+(i-1)*xdiis*nroots_reduced+i-1)=1.0d0
             enddo
             do i=1,xdiis*xdiis*nroots_reduced*nroots_reduced
               dbl_mb(k_h_eff+i-1)=0.0d0
             enddo
c
c    declare h_eff c_eff - MA 
c
      cpu=-util_cpusec()
      wall=-util_wallsec()
c
      do iroot=1,nroots_reduced ! x-x-x-x-x-x-x-x-x
c
c - T1/X1 LOCALIZATION ----------
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,size_x1,'x1_local',
     1      l_x1_local,k_x1_local))
     1      call errquit('x1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
        call ma_zero(dbl_mb(k_x1_local),size_x1)
cc        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
cc        call ga_get(x1(iroot),1,size_x1,1,1,dbl_mb(k_x1_local),1)
        call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
        call get_block(x1(iroot),dbl_mb(k_x1_local),size_x1,0)
c -------------------------------
         if(ccsd_var.eq.'cc2') then
          call cc2_x1(d_f1,xp1(iroot),k_t1_local,d_t2,d_v2,
     &      k_x1_local,x2(iroot),
     &      k_f1_offset,k_x1_offset,k_t1_offset,k_t2_offset,
     &      k_v2_offset,k_x1_offset,k_x2_offset)
         else
          call eomccsd_x1_n(d_f1,xp1(iroot),k_t1_local,d_t2,d_v2,
     1      k_x1_local,x2(iroot),k_f1_offset,k_x1_offset,
     2      k_t1_offset,k_t2_offset,k_v2_offset,
     3      k_x1_offset,k_x2_offset,
     4      ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5      rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6      rec_d_6,rec_d_7,rec_d_8,rec_d_9)
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
           call eomccsd_x2_n(d_f1,xp2(iroot),k_t1_local,d_t2,d_v2,
     1       k_x1_local,x2(iroot),k_f1_offset,k_x2_offset,
     2       k_t1_offset,k_t2_offset,k_v2_offset,
     3       k_x1_offset,k_x2_offset,size_x1,size_x2,
     4       ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9)
         end if
         if(ccsd_var.eq.'lr-ccsd') then
           call tce_j1_x1_scaling(xp1(iroot),k_x1_offset,
     1          d_j1,k_j1_offset)
         end if
c
         ister=1
c - T1/X1 LOCALIZATION --
         if(.not.MA_POP_STACK(l_x1_local))
     &      call errquit('l_x1_local',4,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c -----------------------------
         enddo ! x-x-x-x-x-x-x-x-x
c   calculate  H-bar delta x
c
           call form_heff_mr(k_x1_offset,k_x2_offset,size_x1,size_x2,
     &        dbl_mb(k_c_eff),dbl_mb(k_h_eff),dbl_mb(k_omegax),
     &        imicro,xdiis,iter)
c
c
          do imicro=2,xdiis !microloop(imicro)
c
          do iroot=1,nroots_reduced
           ipos=(imicro-1)*nroots_reduced+iroot
           call tce_zero(x1(ipos),size_x1)
           call tce_zero(x2(ipos),size_x2)
           call tce_zero(xp1(ipos),size_x1)
           call tce_zero(xp2(ipos),size_x2)
          enddo
c
          do iroot=1,nroots_reduced
           dbl_mb(k_residual+iroot-1)=0.0d0
          enddo
c
          do iroot=1,nroots_reduced
           call form_correction_mr(k_x1_offset,k_x2_offset,
     1                  size_x1,size_x2,
     1                  dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     1                  dbl_mb(k_omegax+iroot-1),
     1                  imicro,xdiis,dbl_mb(k_residual+iroot-1),iroot)
          enddo
c
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
c
           converged=.true.
           do iroot=1,nroots_reduced
             if(dabs(dbl_mb(k_residual+iroot-1)).gt.thresh) 
     &             converged=.false. 
           enddo 
c
            if(converged) then
             do iroot=1,nroots_reduced
              call form_best_guess_mr(k_x1_offset,k_x2_offset,
     &                size_x1,size_x2,
     &                dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     &                imicro-1,
     &                xdiis,iroot)
             enddo
c
               do iroot=1,nroots_reduced
                if(nodezero)
     &                  write(LuOut,*)'EOMCCSD root nr. ',iroot
                call tce_print_x1(xc1(iroot),k_x1_offset,0.10d0,irrep_x)
                call tce_print_x2(xc2(iroot),k_x2_offset,0.10d0,irrep_x)
               enddo
              go to 112
            end if
c all together
            call deltax_ort_norm_mr(k_x1_offset,k_x2_offset,
     1                              size_x1,size_x2,
     1                              imicro,xdiis)
c
c
c
         do iroot=1,nroots_reduced ! --- iroot
c
          ipos=(imicro-1)*nroots_reduced+iroot
c
c - T1/X1 LOCALIZATION ----------
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,size_x1,'x1_local',
     1      l_x1_local,k_x1_local))
     1      call errquit('x1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
        call ma_zero(dbl_mb(k_x1_local),size_x1)
cc        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
cc        call ga_get(x1(ipos),1,size_x1,1,1,dbl_mb(k_x1_local),1)
        call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
        call get_block(x1(ipos),dbl_mb(k_x1_local),size_x1,0)
c -------------------------------
         if(ccsd_var.eq.'cc2') then
          call cc2_x1(d_f1,xp1(ipos),k_t1_local,d_t2,d_v2,
     &      k_x1_local,x2(ipos),
     &      k_f1_offset,k_x1_offset,k_t1_offset,k_t2_offset,
     &      k_v2_offset,k_x1_offset,k_x2_offset)
         else
          call eomccsd_x1_n(d_f1,xp1(ipos),k_t1_local,d_t2,d_v2,
     1      k_x1_local,x2(ipos),k_f1_offset,k_x1_offset,
     2      k_t1_offset,k_t2_offset,k_v2_offset,
     3      k_x1_offset,k_x2_offset,
     4      ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5      rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6      rec_d_6,rec_d_7,rec_d_8,rec_d_9)
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
           call eomccsd_x2_n(d_f1,xp2(ipos),k_t1_local,d_t2,d_v2,
     1       k_x1_local,x2(ipos),k_f1_offset,k_x2_offset,
     2       k_t1_offset,k_t2_offset,k_v2_offset,
     3       k_x1_offset,k_x2_offset,size_x1,size_x2,
     4       ister,rec_s_1,rec_s_2,rec_s_3,rec_s_4,rec_s_5,
     5       rec_d_1,rec_d_2,rec_d_3,rec_d_4,rec_d_5,
     6       rec_d_6,rec_d_7,rec_d_8,rec_d_9)
         end if
         if(ccsd_var.eq.'lr-ccsd') then
           call tce_j1_x1_scaling(xp1(ipos),k_x1_offset,
     1          d_j1,k_j1_offset)
         end if
c
         ister=1
c - T1/X1 LOCALIZATION --
         if(.not.MA_POP_STACK(l_x1_local))
     &      call errquit('l_x1_local',4,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c -----------------------
c
        enddo ! --- iroot 
c
c             calculate  H-bar delta x
c
           call form_heff_mr(k_x1_offset,k_x2_offset,size_x1,size_x2,
     &        dbl_mb(k_c_eff),dbl_mb(k_h_eff),dbl_mb(k_omegax),
     &        imicro,xdiis,iter)
c
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
c
          enddo   !microloop(imicro)
c
             do iroot=1,nroots_reduced
              call form_best_guess_mr(k_x1_offset,k_x2_offset,
     &                size_x1,size_x2,
     &                dbl_mb(k_c_eff+(iroot-1)*xdiis*nroots_reduced),
     &                xdiis,
     &                xdiis,iroot)
             enddo
c
               go to 111
c
 112          continue !from converged
c

         if(.not.MA_POP_STACK(l_residual))
     &      call errquit('l_residual',4,MA_ERR)
         if(.not.MA_POP_STACK(l_c_eff))
     &      call errquit('l_c_eff',4,MA_ERR)
         if(.not.MA_POP_STACK(l_h_eff))
     &      call errquit('l_h_eff',4,MA_ERR)
c
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
c
c
c
c xc1 xc2 are deleted by tce_eom_tidy (later on)
c
c
           if(ccsd_var.eq.'lr-ccsd')then
            call deletefile(d_j2)
            if (.not.ma_pop_stack(l_j2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
            call deletefile(d_j1)
            if (.not.ma_pop_stack(l_j1_offset))
     1        call errquit("tce_energy: MA problem",37,MA_ERR)
           end if
c
c
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
c
         do ivec=1,maxtrials
           if(x1_ini_exist(ivec))  then
            call deletefile(x1_ini(ivec))
            x1_ini_exist(ivec) = .false.
           end if
         enddo
c
c
c
         end if ! ------- EOM-SOLVER ---------------
c
c
c
c -------------
c EOM-CCSD    left
c -------------
  200       continue
c
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
c
            if (left) then
            call tce_eom_yguess(needt1,needt2,needt3,needt4,
     1        size_y1,size_y2,size_y3,size_y4,
     2        k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset)
            if (nytrials.eq.0) goto 400
            modelname = "EOM-CCSD left-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
c           call hbar_init
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
c               call hbar_tidy
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
                    if (nodezero) write(LuOut,9430) "EOM-CCSD",
     1                      dipole_exc(1),dipole_exc(1)*debye,
     2                      dipole_exc(2),dipole_exc(2)*debye,
     3                      dipole_exc(3),dipole_exc(3)*debye
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
c                    irrep_x = irrep   !! Caution !!
c                    call tce_print_x1(xc1(ivec),k_x1_offset,printtol)
c                    call tce_print_x2(xc2(ivec),k_x2_offset,printtol)
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
c
c->d3p975
c
c   eom-ccsd 1-particle reduced density matrix/ces
c
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
c 
c<-d3p975
c
                enddo
                modelname = "EOM-CCSD"
                goto 400
              endif
            enddo
            call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
            endif


c -------------------------
c CR-EOMCCSD(T)
c -------------------------
        if ((perturbative.eq.'creom_(t)').or.
     &      (perturbative.eq.'emb1').or.(perturbative.eq.'emb2')) then
        irrep_c = irrep_x  !important
        go to 8922
ccccx      if(irrep_x.eq.irrep_g) then !ground-state calculations
c      if(irrep_x.eq.irrep_t) then !ground-state calculations
c            call tce_x2_offset(l_c2_offset,k_c2_offset,size_c2)
c            call tce_filename('c2',filename)
c            call createfile(filename,d_c2,size_c2)
ccc            call gatoeaf(d_c2)
c            call tce_zero(d_c2,size_c2)
c       call tce_e_offset(l_e_offset,k_e_offset,size_e)
c       call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
c       call tce_e_offset(l_ex2_offset,k_ex2_offset,size_ex2)
c          call tce_filename('e',filename)
c          call createfile(filename,d_e,size_e)
c          call tce_zero(d_e,size_e)
c          call tce_filename('ex1',filename)
c          call createfile(filename,d_ex1,size_ex1)
c          call tce_zero(d_ex1,size_ex1)
c          call tce_filename('ex2',filename)
c          call createfile(filename,d_ex2,size_ex2)
c          call tce_zero(d_ex2,size_ex2)
c      pt3_6=0.0d0
c      call cr_eomccsd_t_gr(d_t1,k_t1_offset,d_t2,k_t2_offset,
c     2                     d_f1,k_f1_offset,d_v2,k_v2_offset,
c     3                     d_e,k_e_offset,
c     4                     d_ex1,k_ex1_offset,size_ex1,
c     5                     d_ex2,k_ex2_offset,size_ex2,
c     7                     d_c2,k_c2_offset,size_c2,
c     8                     pt3_6)
c       energy_grcr=ref+corr+pt3_6
c       if (nodezero) then
c       write(LuOut,9080) "ground-state CR-CCSD(T)IA",ref+corr+pt3_6
c       write(LuOut,*)
c       call util_flush(LuOut)
c       end if
cc
c          call deletefile(d_ex2)
c          call deletefile(d_ex1)
c          call deletefile(d_e)
c          if (.not.ma_pop_stack(l_ex2_offset))
c     1     call errquit("tce_energy: MA problem",24,MA_ERR)
c          if (.not.ma_pop_stack(l_ex1_offset))
c     1     call errquit("tce_energy: MA problem",24,MA_ERR)
c          if (.not.ma_pop_stack(l_e_offset))
c     1     call errquit("tce_energy: MA problem",24,MA_ERR)
cc
c            call deletefile(d_c2)
c            if (.not.ma_pop_stack(l_c2_offset))
c     1        call errquit("tce_energy: MA problem",36,MA_ERR)
cc
c      end if
c UWAGA
c ACHTUNG
c BHUMAHE
c ATTENTION
 8922  continue
c
        DO 3443,ivec = 1,nroots_reduced
          cpu = - util_cpusec()
          wall = - util_wallsec()
c d_e,d_ex1,d_ex2 global arrays
c
            call tce_x1_offset(l_c1_offset,k_c1_offset,size_c1)
            call tce_filename('c1',filename)
            call createfile(filename,d_c1,size_c1)
cc            call gatoeaf(d_c1)
            call tce_zero(d_c1,size_c1)
c
            call tce_x2_offset(l_c2_offset,k_c2_offset,size_c2)
            call tce_filename('c2',filename)
            call createfile(filename,d_c2,size_c2)
cc            call gatoeaf(d_c2)
            call tce_zero(d_c2,size_c2)
c
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
      call cr_eomccsd_t(xc1(ivec),k_x1_offset,xc2(ivec),k_x2_offset,
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
c
          call deletefile(d_ex2)
          call deletefile(d_ex1)
          call deletefile(d_e)
          if (.not.ma_pop_stack(l_ex2_offset))
     1     call errquit("tce_energy: MA problem",24,MA_ERR)
          if (.not.ma_pop_stack(l_ex1_offset))
     1     call errquit("tce_energy: MA problem",24,MA_ERR)
          if (.not.ma_pop_stack(l_e_offset))
     1     call errquit("tce_energy: MA problem",24,MA_ERR)
c
          cpu = cpu + util_cpusec()
          wall = wall + util_wallsec()
       if(perturbative.eq.'creom_(t)') then
          if (nodezero) then
           write(LuOut,9080)"CR-EOMCCSD(T)IA",ref+corr+excit+pt3_6
c UWAGA
c ACHTUNG
c BHUMAHE
c ATTENTION
c           write(LuOut,9085)"CR-EOMCCSD(T)IA - CR-CCSD(T)IA",
c     &                      (ref+corr+excit+pt3_6-energy_grcr)*au2ev
           write(LuOut,9085)"delta-CR-EOMCCSD(T),IA        ",
     &                      (excit+pt3_6)*au2ev
           write(LuOut,*)
           call util_flush(LuOut)
          endif
         end if
       if(perturbative.eq.'emb1') then
          if (nodezero) then
           write(LuOut,9080)"r-CR-EOMCCSD(T)-I",ref+corr+excit+pt3_6
c UWAGA
c ACHTUNG
c BHUMAHE
c ATTENTION
           write(LuOut,9085)"delta-r-CR-EOMCCSD(T)-I        ",
     &                      (excit+pt3_6)*au2ev
           write(LuOut,*)
           call util_flush(LuOut)
          endif
         end if
       if(perturbative.eq.'emb2') then
          if (nodezero) then
           write(LuOut,9080)"r-CR-EOMCCSD(T)-II",ref+corr+excit+pt3_6
c UWAGA
c ACHTUNG
c BHUMAHE
c ATTENTION
           write(LuOut,9085)"delta-r-CR-EOMCCSD(T)-II        ",
     &                      (excit+pt3_6)*au2ev
           write(LuOut,*)
           call util_flush(LuOut)
          endif
         end if

c redefining the vector with excitation energies
          dbl_mb(k_omegax+ivec-1)=dbl_mb(k_omegax+ivec-1)+pt3_6
c
            call deletefile(d_c2)
            if (.not.ma_pop_stack(l_c2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
            call deletefile(d_c1)
            if (.not.ma_pop_stack(l_c1_offset))
     1        call errquit("tce_energy: MA problem",37,MA_ERR)
c
 3443   CONTINUE
        end if !cr-eom
c
c
          elseif (model.eq.'ccsdt') then
c -------------
c EOM-CCSDT   right
c -------------
            call tce_eom_xguess(needt1,needt2,needt3,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
cedo            if (nxtrials.eq.0) goto 210
            if (nxtrials.eq.0) then
               ccsdtleft=.true.
               goto 400
            endif
ccxsg
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
ccxsg
            modelname = "EOM-CCSDT right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
c           call hbar_init
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
     4          dbl_mb(k_omegax),dbl_mb(k_residual),k_hbar,iter)
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
c               call hbar_tidy
                if (nodezero) then
                  write(LuOut,9240)
                  call util_flush(LuOut)
                endif
c ---
            if(nodezero) then
             write(6,*)'largest EOMCCSDT amplitudes: R1 and R2'
             call util_flush(6)
            end if
             do jvec=1,nroots_reduced
c              call tce_print_x1(xc1(jvec),k_x1_offset,0.10d0)
c              call tce_print_x2(xc2(jvec),k_x2_offset,0.10d0)
              call tce_print_x1(xc1(jvec),k_x1_offset,printtol,irrep_x)
              call tce_print_x2(xc2(jvec),k_x2_offset,printtol,irrep_x)
              call util_flush(LuOut)
             enddo
c ---
                ccsdtleft=.true.
                goto 400
cedo                goto 210
              endif
            enddo
            call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
c
c
c
c
          elseif (model.eq.'ccsdta') then
c -----------------------
c EOM-CCSDt  right
c -----------------------
            call tce_eom_xguess(needt1,needt2,needt3a,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
cccx            if (nxtrials.eq.0) goto 210
ccxsg
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
ccxsg
            modelname = "EOM-CCSDt right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
c           call hbar_init
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
     4          dbl_mb(k_omegax),dbl_mb(k_residual),k_hbar,iter)
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
c               call hbar_tidy
                if (nodezero) then
                  write(LuOut,9240)
                  call util_flush(LuOut)
                endif
c ---
            if(nodezero) then
             write(6,*)'largest EOMCCSDt amplitudes: R1 and R2 '
             call util_flush(6)
            end if
             do jvec=1,nroots_reduced
c              call tce_print_x1(xc1(jvec),k_x1_offset,0.10d0)
c              call tce_print_x2(xc2(jvec),k_x2_offset,0.10d0)
              call tce_print_x1(xc1(jvec),k_x1_offset,printtol,irrep_x)
              call tce_print_x2(xc2(jvec),k_x2_offset,printtol,irrep_x)
              call util_flush(LuOut)
             enddo
c ---
                ccsdtleft=.true.
                goto 400
cedo                goto 210
              endif
            enddo
            call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
c
c
c
c
          elseif (model.eq.'ccsdtq') then
#if !defined(CCSDTQ)
      CALL ERRQUIT('CCSDTQ not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
c -------------
c EOM-CCSDTQ  right
c -------------
            call tce_eom_xguess(needt1,needt2,needt3,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
            if (nxtrials.eq.0) goto 220
ccxsg
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
ccxsg
            modelname = "EOM-CCSDTQ right-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
c           call hbar_init
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
     4          dbl_mb(k_omegax),dbl_mb(k_residual),k_hbar,iter)
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
c               call hbar_tidy
		if (nodezero) then
		  write(LuOut,9240)
		  call util_flush(LuOut)
		endif
		goto 220
	      endif
	    enddo
	    call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
c -------------
c EOM-CCSDTQ  left
c -------------
  220       continue
	    if (left) then
	    call tce_eom_yguess(needt1,needt2,needt3,needt4,
     1        size_y1,size_y2,size_y3,size_y4,
     2        k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset)
	    if (nytrials.eq.0) goto 400
	    modelname = "EOM-CCSDTQ left-hand side"
	    if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
c           call hbar_init
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
cccx                  call eomccsdtq_y1(d_f1,yp1(ivec),d_t1,d_t2,d_t3,d_t4,
cccx     1              d_v2,y1(ivec),y2(ivec),y3(ivec),y4(ivec),
cccx     2              k_f1_offset,k_y1_offset,k_t1_offset,k_t2_offset,
cccx     3              k_t3_offset,k_t4_offset,k_v2_offset,k_y1_offset,
cccx     4              k_y2_offset,k_y3_offset,k_y4_offset)
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
cccx                  call eomccsdtq_y2(d_f1,yp2(ivec),d_t1,d_t2,d_t3,d_t4,
cccx     1              d_v2,y1(ivec),y2(ivec),y3(ivec),y4(ivec),
cccx     2              k_f1_offset,k_y2_offset,k_t1_offset,k_t2_offset,
cccx     3              k_t3_offset,k_t4_offset,k_v2_offset,k_y1_offset,
cccx     4              k_y2_offset,k_y3_offset,k_y4_offset)
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
cccx                    call eomccsdtq_y3(d_f1,yp3(ivec),d_t1,d_t2,d_t3,
cccx     1              d_v2,y1(ivec),y2(ivec),y3(ivec),y4(ivec),
cccx     2              k_f1_offset,k_y3_offset,k_t1_offset,k_t2_offset,
cccx     3              k_t3_offset,k_v2_offset,k_y1_offset,
cccx     4              k_y2_offset,k_y3_offset,k_y4_offset)
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
cccx                  call eomccsdtq_y4(d_f1,yp4(ivec),d_t1,d_t2,
cccx     1              d_v2,y2(ivec),y3(ivec),y4(ivec),
cccx     2              k_f1_offset,k_y4_offset,k_t1_offset,k_t2_offset,
cccx     3              k_v2_offset,k_y2_offset,k_y3_offset,k_y4_offset)
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
c               call hbar_tidy
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
cccx                    call eomccsdtq_denominator(d_d0,xc1(ivec),xc2(ivec),
cccx     1                xc3(ivec),xc4(ivec),yc1(jvec),yc2(jvec),yc3(jvec),
cccx     2                yc4(jvec),k_d0_offset,k_x1_offset,k_x2_offset,
cccx     3                k_x3_offset,k_x4_offset,k_y1_offset,k_y2_offset,
cccx     4                k_y3_offset,k_y4_offset)
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
cccx                      call eomccsdtq_density1(d_d1(axis),d_d0,d_t1,d_t2,
cccx     1                  d_t3,d_t4,d_x0,xc1(ivec),xc2(ivec),xc3(ivec),
cccx     2                  xc4(ivec),d_y0,yc1(jvec),yc2(jvec),yc3(jvec),
cccx     3                  yc4(jvec),k_d1_offset(axis),k_d0_offset,
cccx     4                  k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
cccx     5                  k_x0_offset,k_x1_offset,k_x2_offset,k_x3_offset,
cccx     6                  k_x4_offset,k_y0_offset,k_y1_offset,k_y2_offset,
cccx     7                  k_y3_offset,k_y4_offset)
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
cccx                      call eomccsdtq_density1(d_d1(axis),d_d0,d_t1,d_t2,
cccx     1                  d_t3,d_t4,d_x0,xc1(ivec),xc2(ivec),xc3(ivec),
cccx     2                  xc4(ivec),d_y0,d_lambda1,d_lambda2,d_lambda3,
cccx     3                  d_lambda4,k_d1_offset(axis),k_d0_offset,
cccx     4                  k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
cccx     5                  k_x0_offset,k_x1_offset,k_x2_offset,k_x3_offset,
cccx     6                  k_x4_offset,k_y0_offset,k_l1_offset,k_l2_offset,
cccx     7                  k_l3_offset,k_l4_offset)
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
cccx                      call eomccsdtq_density1(d_d1(axis),d_d0,d_t1,d_t2,
cccx     1                  d_t3,d_t4,d_x0,d_r1,d_r2,d_r3,
cccx     2                  d_r4,d_y0,yc1(jvec),yc2(jvec),yc3(jvec),
cccx     3                  yc4(jvec),k_d1_offset(axis),k_d0_offset,
cccx     4                  k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
cccx     5                  k_x0_offset,k_r1_offset,k_r2_offset,k_r3_offset,
cccx     6                  k_r4_offset,k_y0_offset,k_y1_offset,k_y2_offset,
cccx     7                  k_y3_offset,k_y4_offset)
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
                    if (nodezero) write(LuOut,9430) "EOM-CCSDTQ",
     1                      dipole_exc(1),dipole_exc(1)*debye,
     2                      dipole_exc(2),dipole_exc(2)*debye,
     3                      dipole_exc(3),dipole_exc(3)*debye
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
                    irrep_x = irrep   !! Caution !!
c                    call tce_print_x1(xc1(ivec),k_x1_offset,printtol)
c                    call tce_print_x2(xc2(ivec),k_x2_offset,printtol)
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
c
c         --------------------
c         Write energy to RTDB
c         --------------------
c
  400     continue
c -------------
c EOM-CCSDT   left
c -------------
          if(ccsdtleft) then
            if (left) then
            call tce_eom_yguess(needt1,needt2,needt3,needt4,
     1        size_y1,size_y2,size_y3,size_y4,
     2        k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset)
            if (nytrials.eq.0) goto 4400
            modelname = "EOM-CCSDT left-hand side"
            if (nodezero) write(LuOut,9220)
     1        modelname(1:inp_strlen(modelname))
c           call hbar_init
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
c               call hbar_tidy
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
                    if (nodezero) write(LuOut,9430) "EOM-CCSDT",
     1                      dipole_exc(1),dipole_exc(1)*debye,
     2                      dipole_exc(2),dipole_exc(2)*debye,
     3                      dipole_exc(3),dipole_exc(3)*debye
                    call dratoga(xc1(ivec))
                    call dratoga(xc2(ivec))
                    irrep_x = irrep   !! Caution !!
c                    call tce_print_x1(xc1(ivec),k_x1_offset,printtol)
c                    call tce_print_x2(xc2(ivec),k_x2_offset,printtol)
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
         endif
c     done with 2nd ccsdtleft
 4400    ccsdtleft=.false.
          call tce_eom_tidy
          if ((targetsym.eq.irrepname)
     1      .and.(target.le.nroots_reduced)) then
            if (nodezero.and.util_print('excited state',print_low)) then
              write(LuOut,*)
              write(LuOut,9300) target
              write(LuOut,9340) targetsym
              write(LuOut,*)
              write(LuOut,9310) modelname(1:inp_strlen(modelname)),
     1          ref+corr
              write(LuOut,9320) modelname(1:inp_strlen(modelname)),
     1          dbl_mb(k_omegax+target-1)
              write(LuOut,9330) modelname(1:inp_strlen(modelname)),
     1          ref+corr+dbl_mb(k_omegax+target-1)
            endiF
            if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,
     1        ref+corr+dbl_mb(k_omegax+target-1)))
     2        call errquit('tce_energy: RTDB problem',0,RTDB_ERR)
            tce_energy=.true.
          elseif ((targetsym.eq.'none').and.(irrep_x.eq.0)
     1      .and.(target.le.nroots_reduced)) then
            if (nodezero.and.util_print('excited state',print_low)) then
ccx              write(LuOut,*)
ccx              write(LuOut,9310) modelname(1:inp_strlen(modelname)),
ccx     1          ref+corr
ccx              write(LuOut,9320) modelname(1:inp_strlen(modelname)),
ccx     1          dbl_mb(k_omegax+target-1)
ccx              write(LuOut,9330) modelname(1:inp_strlen(modelname)),
ccx     1          ref+corr+dbl_mb(k_omegax+target-1)
            endif
            if (.not.rtdb_put(rtdb,'tce:omega',mt_dbl,1,
     1        dbl_mb(k_omegax+target-1)*au2ev))
     2        call errquit('tce_energy: RTDB problem',0,RTDB_ERR)
            if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,
     1        ref+corr+dbl_mb(k_omegax+target-1)))
     2        call errquit('tce_energy: RTDB problem',0,RTDB_ERR)
            tce_energy=.true.
          endif
c link to QM/MM (no symmetry OR a-symmetry) -------
          if(irrep.eq.0) then
            if (.not.rtdb_put(rtdb,'tce:eneex',mt_dbl,nroots_reduced,
     1        dbl_mb(k_omegax)))
     2        call errquit('tce_energy: RTDB problem',0,RTDB_ERR)
          end if
c -------------------------------------------------
c
c         ------------------------------------
c         Terminate excited state calculations
c         ------------------------------------
c
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
ccc
          if (needt3a) then
cccx           if (left) then
cccx              call deletefile(d_ry3)
cccx              if (.not.ma_pop_stack(l_y3_offset))
cccx     1          call errquit("tce_energy: MA problem",35,MA_ERR)
cccx            endif
            call deletefile(d_rx3)
            if (.not.ma_pop_stack(l_x3_offset))
     1        call errquit("tce_energy: MA problem",35,MA_ERR)
          endif
ccc
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
c new EOMCC solver
       if(.not.((model.ne.'ccsd').and.(eom_solver.eq.2))) then
      if (.not.ma_pop_stack(l_hbar))
     1  call errquit('tce_eom_xdiagon: MA problem',12,MA_ERR)
      end if
          endif
        enddo
      endif
c
c     ===========================
c     End of all TCE calculations
c     ===========================
c
c->d3p975
      do i = 1, ipol
        if (.not.ga_destroy(g_movecs(i)))
     1    call errquit('tce_energy: GA problem',0,GA_ERR)
      enddo
c<-d3p975
c
      if (io_v2.eq.IOSF) then
        call deletefile_sf(d_v2)
      else
        call deletefile(d_v2)
      endif
      if (left) then
        if (multipole.ge.3) then
          do axis=10,1,-1
            call deletefile(d_o1(axis))
          enddo
        endif
        if (multipole.ge.2) then
          do axis=6,1,-1
            call deletefile(d_q1(axis))
          enddo
        endif
        do axis=3,1,-1
          call deletefile(d_d1(axis))
        enddo
      endif
      call deletefile(d_f1)
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
        if (lineresp.or.imagresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_tr4(axis))
            call deletefile(d_tr4(axis+3))
            if (imagresp) call deletefile(d_tr4(axis+6))
            if (.not.ma_pop_stack(l_tr4_offset(axis)))
     1        call errquit("tce_energy: MA problem",89,MA_ERR)
          endif
          enddo ! axis
        endif ! lineresp.or.imagresp
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
        if (lineresp.or.imagresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_tr3(axis))
            call deletefile(d_tr3(axis+3))
            if (imagresp) call deletefile(d_tr3(axis+6))
            if (.not.ma_pop_stack(l_tr3_offset(axis)))
     1        call errquit("tce_energy: MA problem",88,MA_ERR)
          endif
          enddo ! axis
        endif ! lineresp.or.imagresp
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
cccx        if (left) then
cccx          call deletefile(d_lambda3)
cccx          if (.not.ma_pop_stack(l_l3_offset))
cccx     1      call errquit("tce_energy: MA problem",25,MA_ERR)
cccx        endif
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
        if (lineresp.or.imagresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_tr2(axis))
            call deletefile(d_tr2(axis+3))
            if (imagresp) call deletefile(d_tr2(axis+6))
            if (.not.ma_pop_stack(l_tr2_offset(axis)))
     1        call errquit("tce_energy: MA problem",87,MA_ERR)
          endif
          enddo ! axis
        endif ! lineresp.or.imagresp
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
        if (lineresp.or.imagresp) then
          do axis = 3, 1, -1
          if (respaxis(axis)) then
            call deletefile(d_tr1(axis))
            call deletefile(d_tr1(axis+3))
            if (imagresp) call deletefile(d_tr1(axis+6))
            if (.not.ma_pop_stack(l_tr1_offset(axis)))
     1        call errquit("tce_energy: MA problem",86,MA_ERR)
          endif
          enddo ! axis
        endif ! lineresp.or.imagresp
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
       if (.not.ma_pop_stack(l_v2spin_offset))
     1  call errquit("tce_energy: MA problem",59,MA_ERR)
c--       if (.not.ma_pop_stack(l_o2ot_offset))
c--     1  call errquit("tce_energy: MA problem",58,MA_ERR)
c--       if (.not.ma_pop_stack(l_o2o_offset))
c--     1  call errquit("tce_energy: MA problem",57,MA_ERR)
       if (.not.ma_pop_stack(l_v2_alpha_offset))
     1  call errquit("tce_energy: MA problem",56,MA_ERR)
      end if
      if(.not.intorb) then
      if (.not.ma_pop_stack(l_v2_offset))
     1  call errquit("tce_energy: MA problem",2,MA_ERR)
      end if
      if (.not.ma_pop_stack(l_f1_offset))
     1  call errquit("tce_energy: MA problem",3,MA_ERR)
      if (left) then
        if (multipole.ge.3) then
          do axis=10,1,-1
            if (.not.ma_pop_stack(l_o1_offset(axis)))
     1        call errquit("tce_energy: MA problem",3030+axis,MA_ERR)
          enddo
        endif
        if (multipole.ge.2) then
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
      end if
c ---
      if(activecalc.or.(.not.intorb)) THEN
      if (.not.ma_pop_stack(l_active))
     1  call errquit("tce_energy: MA problem",34,MA_ERR)
      end if
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
c
c     ===============
c     Destroy a mutex
c     ===============
c
      if (.not.ga_destroy_mutexes())
     1  call errquit('tce_energy: GA problem',1,GA_ERR)
c
c     =========
c     Terminate
c     =========
c
      call tce_tidy(rtdb)
      call util_print_pop
      return
c
c     ======
c     Format
c     ======
c
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
c 9126 format(1x,'imaxsize = ',2x,i6)
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
     5  1x,'------------------------------------')
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
c
c     Ordering in mptval (from defNxyz):
c     xx, xy, xz, yy, yz, zz
 9432 format(/,1x,A,' quadrupole moments / hartree & Buckingham',/,
     1  1x,'--------------------------------',/
     2  1x,'XX ',2f15.7,/
     3  1x,'XY ',2f15.7,/
     4  1x,'XZ ',2f15.7,/
     5  1x,'YY ',2f15.7,/
     6  1x,'YZ ',2f15.7,/
     7  1x,'ZZ ',2f15.7,/
     8  1x,'--------------------------------')
c
c     Ordering in mptval (from defNxyz):
c     xxx, xxy, xxz, yyx, xyz, zzx, yyy, zzy, zzz
 9433 format(/,1x,A,' octupole moments / hartree & Octos',/,
     1  1x,'--------------------------------',/
     2  1x,'XXX ',2f15.7,/
     3  1x,'YYY ',2f15.7,/
     4  1x,'ZZZ ',2f15.7,/
     5  1x,'XXY ',2f15.7,/
     6  1x,'XXZ ',2f15.7,/
     7  1x,'YYX ',2f15.7,/
     8  1x,'YYZ ',2f15.7,/
     9  1x,'ZZX ',2f15.7,/
     1  1x,'ZZY ',2f15.7,/
     2  1x,'XYZ ',2f15.7,/
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
c Use 9431 for frequency
 9502 format(1x,'-----------------------------------------------')
 9503 format(1x,'beta(',a1,',',a1,',',a1,') = ',f18.9,a)
 9504 format(1x,'Static Hyperpolarizability')
 9505 format(1x,'Second Harmonic Generation (SHG)')
 9506 format(1x,'Optical Rectification (OR)')
 9507 format(1x,'Electro-Optic Pockels Effect (EOPE)')
 9508 format(1x,'beta(',f1.5,',',f1.5,',',f1.5,')')
      end
