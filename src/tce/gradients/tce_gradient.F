      logical function tce_gradient(rtdb)
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
c     
c     CI, CC, & MBPT
c     
      integer rtdb              ! Run-time database
      logical nodezero          ! True if node 0
      logical recompf           ! True if recompute Fock
      double precision cpu      ! CPU sec counter
      double precision wall     ! WALL sec counter
      integer g_ao1e(2)         ! GA handle for AO Fock matrices
      integer d_ao2e            ! SF handle for AO 2e integrals
      integer d_f1              ! SF handle for MO 1e integrals
      integer l_f1_offset       ! Offset for 1e integral file
      integer k_f1_offset       ! Offset for 1e integral file
      integer size_1e           ! File size in doubles
      integer d_v2              ! SF handle for MO 2e integrals
      integer l_v2_offset       ! Offset for 2e integral file
      integer k_v2_offset       ! Offset for 2e integral file
      integer size_2e           ! File size in doubles
c     --- relevant if intorb = .true. ----
      integer lax_v2_alpha_offset ! Offset for 2e integrals if intorb=.true.
      integer kax_v2_alpha_offset ! Offset for 2e integrals if intorb=.true.
      integer lx_o2o_offset,kx_o2o_offset
      integer lx_o2ot_offset,kx_o2ot_offset
      integer lx_v2spin_offset,kx_v2spin_offset
      integer length
c     ------------------------------------
      integer d_t1              ! SF handle for t1 amplitudes
      integer l_t1_offset       ! Offset for t1 file
      integer k_t1_offset       ! Offset for t1 file
      integer size_t1           ! File size in doubles
      integer d_t2              ! SF handle for t2 amplitudes
      integer l_t2_offset       ! Offset for t2 file
      integer k_t2_offset       ! Offset for t2 file
      integer size_t2           ! File size in doubles
      integer d_t3              ! SF handle for t3 amplitudes
      integer l_t3_offset       ! Offset for t3 file
      integer k_t3_offset       ! Offset for t3 file
      integer size_t3           ! File size in doubles
      integer d_t4              ! SF handle for t4 amplitudes
      integer l_t4_offset       ! Offset for t4 file
      integer k_t4_offset       ! Offset for t4 file
      integer size_t4           ! File size in doubles
      integer d_r1              ! SF handle for r1 amplitudes
      integer l_r1_offset       ! Offset for r1 file
      integer k_r1_offset       ! Offset for r1 file
      integer size_r1           ! File size in doubles
      integer d_r2              ! SF handle for r2 amplitudes
      integer l_r2_offset       ! Offset for r2 file
      integer k_r2_offset       ! Offset for r2 file
      integer size_r2           ! File size in doubles
      integer d_r3              ! SF handle for r3 amplitudes
      integer l_r3_offset       ! Offset for r3 file
      integer k_r3_offset       ! Offset for r3 file
      integer size_r3           ! File size in doubles
      integer d_r4              ! SF handle for r4 amplitudes
      integer l_r4_offset       ! Offset for r4 file
      integer k_r4_offset       ! Offset for r4 file
      integer size_r4           ! File size in doubles
      integer d_e               ! SF handle for e file
      integer l_e_offset        ! Offset for e file
      integer k_e_offset        ! Offset for e file
      integer size_e            ! File size in doubles
      integer d_u1              ! SF handle for u1 file
      integer l_u1_offset       ! Offset for u1 file
      integer k_u1_offset       ! Offset for u1 file
      integer size_u1           ! File size in doubles
      integer d_u2              ! SF handle for u2 file
      integer l_u2_offset       ! Offset for u2 file
      integer k_u2_offset       ! Offset for u2 file
      integer size_u2           ! File size in doubles
      integer d_w1              ! SF handle for w1 file
      integer l_w1_offset       ! Offset for w1 file
      integer k_w1_offset       ! Offset for w1 file
      integer size_w1           ! File size in doubles
      integer d_w2              ! SF handle for w2 file
      integer l_w2_offset       ! Offset for w2 file
      integer k_w2_offset       ! Offset for w2 file
      integer size_w2           ! File size in doubles
      integer d_w3              ! SF handle for w3 file
      integer l_w3_offset       ! Offset for w3 file
      integer k_w3_offset       ! Offset for w3 file
      integer size_w3           ! File size in doubles
      integer d_w4              ! SF handle for w4 file
      integer l_w4_offset       ! Offset for w4 file
      integer k_w4_offset       ! Offset for w4 file
      integer size_w4           ! File size in doubles
      integer d_e2              ! SF handle for e2 file
      integer l_e2_offset       ! Offset for e2 file
      integer k_e2_offset       ! Offset for e2 file
      integer size_e2           ! File size in doubles
      double precision ref      ! Ground state energy
      double precision corr     ! Correlation energy
      double precision mbpt2    ! Correlation energy
      double precision mbpt3    ! Correlation energy
      double precision mbpt4    ! Correlation energy
      double precision r1       ! Residual t1
      double precision r2       ! Residual t2
      double precision r3       ! Residual t3
      double precision r4       ! Residual t4
      double precision residual ! Largest residual
c      integer iter              ! Iteration counter
      integer dummy             ! Dummy argument for DIIS
c     FOR CR-EOMCCSD(T) ---------------------
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
c     
c     ->pengdong.fan
c     
c     for ccsd/eom-ccsd 1-e reduced density matrices
c     
c     integer d_hh, size_hh
c     integer l_ph_offset,k_ph_offset
c     integer d_ph, size_ph
c     integer l_pp_offset,k_pp_offset
c     integer d_pp, size_pp
c     integer l_hp_offset,k_hp_offset
c     integer d_hp, size_hp
c     integer dim_rdm_ao, l_rdm_ao, k_rdm_ao
c     integer nh,np
c     integer dim_rdm_mo_hh,dim_rdm_mo_hp
c     integer dim_rdm_mo_ph,dim_rdm_mo_pp
c     integer l_rdm_mo_hh,k_rdm_mo_hh
c     integer l_rdm_mo_hp,k_rdm_mo_hp
c     integer l_rdm_mo_ph,k_rdm_mo_ph
c     integer l_rdm_mo_pp,k_rdm_mo_pp
c     integer l_mo_p,k_mo_p,l_mo_h,k_mo_h
c     integer dim_mu_h,dim_mu_p,dim_mo_h,dim_mo_p
c     integer l_mu_h,k_mu_h,l_mu_p,k_mu_p
c     integer hole,particle,unitno
c     logical ao_rdm_write
c     external ao_rdm_write
c     double precision res,a(10,7),b(10,10),c(7,10)
c     double precision res_x0,res_y0
      double precision r0
c     
c     <-pengdong.fan      
c     
c     ---------------------------------------
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
      character*255 filename
      character*4 irrepname
      logical dft_energy
      external dft_energy
      logical scf
      external scf
      integer i
      integer j
c     --- debug ---
      integer dup1,dup2,dup3,dup4,dup5,dup6
c     -------------

c     
c     CC Lambda equation
c     
      integer g_aod1(3)         ! GA handle for AO dipole matrices
      integer d_d1(3)           ! MO dipole integral file
      integer l_d1_offset(3)    ! Offset for d1 integral file
      integer k_d1_offset(3)    ! Offset for d1 integral file
      integer size_d1(3)        ! File size in doubles
      character*3 dfilename(3)  ! File name stub
      data dfilename/'d1x','d1y','d1z'/
      character*1 axisname(3)   ! Axis
      data axisname/'X','Y','Z'/
      double precision debye    ! Hartree to Debye
      parameter (debye=-2.541766d0)
      integer d_lambda1         ! Lambda 1 file
      integer d_lambda2         ! Lambda 2 file
      integer d_lambda3         ! Lambda 3 file
      integer d_lambda4         ! Lambda 4 file
      integer d_lr1             ! Lambda 1 residual file
      integer d_lr2             ! Lambda 2 residual file
      integer d_lr3             ! Lambda 3 residual file
      integer d_lr4             ! Lambda 4 residual file
      integer l_l1_offset       ! Offset for Lambda 1 file
      integer k_l1_offset       ! Offset for Lambda 1 file
      integer size_l1           ! File size in doubles
      integer l_l2_offset       ! Offset for Lambda 2 file
      integer k_l2_offset       ! Offset for Lambda 2 file
      integer size_l2           ! File size in doubles
      integer l_l3_offset       ! Offset for Lambda 3 file
      integer k_l3_offset       ! Offset for Lambda 3 file
      integer size_l3           ! File size in doubles
      integer l_l4_offset       ! Offset for Lambda 4 file
      integer k_l4_offset       ! Offset for Lambda 4 file
      integer size_l4           ! File size in doubles
      integer d_d0              ! Dipole moment
      integer l_d0_offset       ! Offset for d0 file
      integer k_d0_offset       ! Offset for d0 file
      integer size_d0           ! File size in doubles
      integer axis              ! X=1, Y=2, & Z=3
      double precision dipole_scf(3) ! Dipole moment
      double precision dipole_cor(3) ! Dipole moment
      double precision dipole_exc(3) ! Dipole moment
      double precision cdipole  ! Dipole moment correlation
      double precision ddotfile
      external ddotfile
      integer sym_abelian_axis  ! Move this and the following to "sym.fh"
      external sym_abelian_axis
c     
c     Excited-state CI, CC, & MBPT
c     
      integer irrep             ! Symmetry loop index
      integer irrep_g           ! Ground state symmetry
      integer d_rx1             ! RHS residual file
      integer d_rx2             ! RHS residual file
      integer d_rx3             ! RHS residual file
      integer d_rx4             ! RHS residual file
      integer d_ry1             ! LHS residual file
      integer d_ry2             ! LHS residual file
      integer d_ry3             ! LHS residual file
      integer d_ry4             ! LHS residual file
      integer d_x0              ! Offset for x0 file
      integer l_x0_offset       ! Offset for x0 file
      integer k_x0_offset       ! Offset for x0 file
      integer size_x0           ! File size in doubles
      integer l_x1_offset       ! Offset for x1 file
      integer k_x1_offset       ! Offset for x1 file
      integer size_x1           ! File size in doubles
      integer l_x2_offset       ! Offset for x2 file
      integer k_x2_offset       ! Offset for x2 file
      integer size_x2           ! File size in doubles
      integer l_x3_offset       ! Offset for x3 file
      integer k_x3_offset       ! Offset for x3 file
      integer size_x3           ! File size in doubles
      integer l_x4_offset       ! Offset for x4 file
      integer k_x4_offset       ! Offset for x4 file
      integer size_x4           ! File size in doubles
      integer d_y0              ! Offset for y0 file
      integer l_y0_offset       ! Offset for y0 file
      integer k_y0_offset       ! Offset for y0 file
      integer size_y0           ! File size in doubles
      integer l_y1_offset       ! Offset for y1 file
      integer k_y1_offset       ! Offset for y1 file
      integer size_y1           ! File size in doubles
      integer l_y2_offset       ! Offset for y2 file
      integer k_y2_offset       ! Offset for y2 file
      integer size_y2           ! File size in doubles
      integer l_y3_offset       ! Offset for y3 file
      integer k_y3_offset       ! Offset for y3 file
      integer size_y3           ! File size in doubles
      integer l_y4_offset       ! Offset for y4 file
      integer k_y4_offset       ! Offset for y4 file
      integer size_y4           ! File size in doubles
      integer ivec,jvec         ! Current trial vector
      integer l_omegax          ! Excitation energy
      integer k_omegax          ! Excitation energy
      integer l_omegay          ! Excitation energy
      integer k_omegay          ! Excitation energy
      integer l_residual        ! Residual
      integer k_residual        ! Residual
c
      double precision omegai
c
      logical converged
      double precision pt3_1    ! [T] energy
      double precision pt3_2    ! (T)-[T] energy
c     LR-CCSD(T) method
      double precision pt3_3,pt3_4,pt3_5,pt3_6
c     --------------------------------------------
c     LR-CCSD(TQ) method
      double precision pt4_lr
c     --------------------------------------------
c     
      double precision au2ev    ! Conversion factor from a.u. to eV
      double precision numerator ! <L|D exp(T)|R>
      double precision denominator ! <L|exp(T)|R>
      double precision transition(3) ! Transition moments
      parameter (au2ev=27.2113961d0)
      character*255 modelname
c     
      nodezero=(ga_nodeid().eq.0)
c     
c     =========================
c     Ground-state HF/DFT first
c     =========================
c     
      tce_gradient = .false.
      if (.not.rtdb_get(rtdb,'tce:reference',mt_int,1,reference)) then
         reference = 1
         if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1        call errquit('tce_gradient: failed to write reference',0,
     2        RTDB_ERR)
      endif
      if (reference.eq.0) then
         if (.not.dft_energy(rtdb)) return
         if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,ref))
     1        call errquit('tce_gradient: failed to get dft energy',0,
     2        RTDB_ERR)
      else if (reference.eq.1) then
         if (.not.scf(rtdb)) return
         if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,ref))
     1        call errquit('tce_gradient: failed to get scf energy',0,
     2        RTDB_ERR)
      else
         call errquit('tce_gradient: illegal reference',reference,
     1        UNKNOWN_ERR)
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
      if (nodezero) call util_flush(LuOut)
c     
c     =================
c     Totally symmetric
c     =================
c     
      irrep_x = 0
      irrep_y = 0
      irrep_c = 0
c     
c     ==============
c     Create a mutex
c     ==============
c     
      if (.not.ga_create_mutexes(1))
     1     call errquit('tce_gradient: GA problem',0,GA_ERR)
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
c     1-e integrals
c     
      cpu = - util_cpusec()
      wall = - util_wallsec()
      if (.not.rtdb_get(rtdb,'tce:recompf',mt_log,1,recompf))
     1     call errquit('tce_gradient: failed reading from rtdb',0,
     2     RTDB_ERR)
      if (recompf) then
         call tce_ao1e_fock2e(rtdb,g_ao1e)
      endif
      do i = 1, ipol
         if (.not.ga_destroy(g_movecs(i)))
     1        call errquit('tce_gradient: GA problem',0,GA_ERR)
      enddo
      call tce_filename('f1',filename)
      call tce_mo1e_offset(l_f1_offset,k_f1_offset,size_1e)
      call createfile(filename,d_f1,size_1e)
      if (recompf) then
         call tce_mo1e(g_ao1e,d_f1,k_f1_offset)
      else
         call tce_mo1e_epsilon(d_f1)
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
         write(LuOut,9090) '1-e',d_f1
         write(LuOut,9020) cpu, wall
         call util_flush(LuOut)
      endif
c     
      if (.not.intorb) then
c     
c     2-e integrals first half
c     
         cpu = - util_cpusec()
         wall = - util_wallsec()
         call tce_ao2e(rtdb,d_ao2e)
         cpu = cpu + util_cpusec()
         wall = wall + util_wallsec()
         if (nodezero.and.util_print('time',print_default))
     1        write(LuOut,9020) cpu, wall
         if (nodezero) call util_flush(LuOut)
c     
c     2-e integrals second half
c     
         cpu = - util_cpusec()
         wall = - util_wallsec()
         call tce_filename('v2',filename)
         call tce_mo2e_offset(l_v2_offset,k_v2_offset,size_2e)
         call createfile(filename,d_v2,size_2e)
         call tce_mo2e(rtdb,d_ao2e,d_v2,k_v2_offset)
         call reconcilefile(d_v2,size_2e)
         cpu = cpu + util_cpusec()
         wall = wall + util_wallsec()
         if (nodezero.and.util_print('mo2e',print_default)) then
            write(LuOut,*)
            write(LuOut,9000) '2-e',size_2e
            write(LuOut,9010) '2-e',filename(1:20)
            write(LuOut,9090) '2-e',d_v2
            write(LuOut,9020) cpu, wall
            call util_flush(LuOut)
         endif
c     
      end if
c     
      if(intorb) then
c     
         if(nodezero) then
            write(6,*)'4-electron integrals stored in orbital form'
         end if
         cpu = - util_cpusec()
         wall = - util_wallsec()
         call tce_filename('v2',filename)
         length = 0 
         call tce_mo2e_offset_intorb(lax_v2_alpha_offset,
     &        kax_v2_alpha_offset,size_2e)
c     --- debug ---
         call ga_sync()
c     -------------
         l_v2_alpha_offset=lax_v2_alpha_offset 
         k_v2_alpha_offset=kax_v2_alpha_offset 
         call createfile(filename,d_v2,size_2e)
         call ga_zero(d_v2)
         if (nodezero) then
            write(LuOut,*)
            write(LuOut,9000) 'v2   ',size_2e
            call util_flush(LuOut)
         endif
         d_v2orb = d_v2
         cpu = - util_cpusec()
         wall = - util_wallsec()
         call tce_mo2e_trans(rtdb,d_v2,kax_v2_alpha_offset,size_2e)
         cpu = cpu + util_cpusec()
         wall = wall + util_wallsec()
         if (nodezero) then
            write(LuOut,9020) cpu, wall
            call util_flush(LuOut)
         endif
         call tce_mo2e_offset_size(length)
c     
         call e2_offset_2_offset(lx_o2o_offset,kx_o2o_offset,
     &        lx_o2ot_offset,kx_o2ot_offset,
     &        lx_v2spin_offset,kx_v2spin_offset,length) 
c     
         l_o2o_offset    = lx_o2o_offset
         k_o2o_offset    = kx_o2o_offset
         l_o2ot_offset   = lx_o2ot_offset
         k_o2ot_offset   = kx_o2ot_offset
         l_v2spin_offset = lx_v2spin_offset
         k_v2spin_offset = kx_v2spin_offset
c     
         k_v2_offset = k_v2spin_offset
c     
      end if
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
      if (model.eq.'ccd') then
         needt2 = .true.
      else if (model.eq.'lccd') then
         needt2 = .true.
      else if (model.eq.'ccsd') then
         needt1 = .true.
         needt2 = .true.
      else if (model.eq.'lccsd') then
         needt1 = .true.
         needt2 = .true.
      else if (model.eq.'ccsdt') then
         needt1 = .true.
         needt2 = .true.
         needt3 = .true.
      else if (model.eq.'ccsdta') then
         needt1  = .true.
         needt2  = .true.
         needt3a = .true.
      else if (model.eq.'ccsdtq') then
         needt1 = .true.
         needt2 = .true.
         needt3 = .true.
         needt4 = .true.
      else if (model.eq.'qcisd') then
         needt1 = .true.
         needt2 = .true.
      else if (model.eq.'cisd') then
         needt1 = .true.
         needt2 = .true.
      else if (model.eq.'cisdt') then
         needt1 = .true.
         needt2 = .true.
         needt3 = .true.
      else if (model.eq.'cisdtq') then
         needt1 = .true.
         needt2 = .true.
         needt3 = .true.
         needt4 = .true.
      else if (model.eq.'mbpt2') then
         needu1 = .true.
         needu2 = .true.
         needt1 = .true.
         needt2 = .true.
      else if (model.eq.'mbpt3') then
         needu1 = .true.
         needu2 = .true.
         needw1 = .true.
         needw2 = .true.
         needt1 = .true.
         needt2 = .true.
      else if (model.eq.'mbpt4') then
         neede2 = .true.
         needu1 = .true.
         needu2 = .true.
         needw1 = .true.
         needw2 = .true.
         needw3 = .true.
         needt1 = .true.
         needt2 = .true.
         needt3 = .true.
      endif
c     
c     t1 amplitudes
c     
      left = .true.
      if (needt1) then
         cpu = - util_cpusec()
         call tce_filename('t1',filename)
         call tce_t1_offset(l_t1_offset,k_t1_offset,size_t1)
         call createfile(filename,d_t1,size_t1)
         call tce_guess_t1(d_t1,k_t1_offset)
         call reconcilefile(d_t1,size_t1)
         cpu = cpu + util_cpusec()
         if (nodezero.and.util_print('t1',print_default)) then
            write(LuOut,*)
            write(LuOut,9000) 't1',size_t1
            write(LuOut,9010) 't1',filename(1:20)
            write(LuOut,9090) 't1',d_t1
            write(LuOut,9020) cpu
            call util_flush(LuOut)
         endif
         if (left) then
            call tce_y1_offset(l_l1_offset,k_l1_offset,size_l1)
            call tce_filename('lambda1',filename)
            call createfile(filename,d_lambda1,size_l1)
            call gatoeaf(d_lambda1)
         endif
      endif
c     
c     t2 amplitudes
c     
      if (needt2) then
         cpu = - util_cpusec()
         call tce_filename('t2',filename)
         call tce_t2_offset(l_t2_offset,k_t2_offset,size_t2)
         call createfile(filename,d_t2,size_t2)
         if(nodezero) then
            write(LuOut,*)'before tce_guess_t2'
            call util_flush(LuOut)
         end if
         call tce_guess_t2(d_v2,k_v2_offset,d_t2,k_t2_offset)
         call reconcilefile(d_t2,size_t2)
         if(nodezero) then
            write(LuOut,*)'after tce_guess_t2'
            call util_flush(LuOut)
         end if
         cpu = cpu + util_cpusec()
         if (nodezero.and.util_print('t2',print_default)) then
            write(LuOut,*)
            write(LuOut,9000) 't2',size_t2
            write(LuOut,9010) 't2',filename(1:20)
            write(LuOut,9090) 't2',d_t2
            write(LuOut,9020) cpu
            call util_flush(LuOut)
         endif
         if (left) then
            call tce_y2_offset(l_l2_offset,k_l2_offset,size_l2)
            call tce_filename('lambda2',filename)
            call createfile(filename,d_lambda2,size_l2)
            call gatoeaf(d_lambda2)
         endif
      endif
c     
c     t3 amplitudes
c     
      if (needt3) then
         cpu = - util_cpusec()
         call tce_filename('t3',filename)
         call tce_t3_offset(l_t3_offset,k_t3_offset,size_t3)
         call createfile(filename,d_t3,size_t3)
         call gatoeaf(d_t3)
         cpu = cpu + util_cpusec()
         if (nodezero.and.util_print('t3',print_default)) then
            write(LuOut,*)
            write(LuOut,9000) 't3',size_t3
            write(LuOut,9010) 't3',filename(1:20)
            write(LuOut,9090) 't3',d_t3
            write(LuOut,9020) cpu
            call util_flush(LuOut)
         endif
         if (left) then
            call tce_y3_offset(l_l3_offset,k_l3_offset,size_l3)
            call tce_filename('lambda3',filename)
            call createfile(filename,d_lambda3,size_l3)
            call gatoeaf(d_lambda3)
         endif
      endif
c     
c     t3a amplitudes
c     
      if (needt3a) then
         cpu = - util_cpusec()
         call tce_filename('t3',filename)
         call tce_t3a_offset(l_t3_offset,k_t3_offset,size_t3)
         call createfile(filename,d_t3,size_t3)
         call gatoeaf(d_t3)
         cpu = cpu + util_cpusec()
         if (nodezero.and.util_print('t3',print_default)) then
            write(LuOut,*)
            write(LuOut,9000) 't3a',size_t3
            write(LuOut,9010) 't3a',filename(1:20)
            write(LuOut,9090) 't3a',d_t3
            write(LuOut,9020) cpu
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
         call gatoeaf(d_t4)
         cpu = cpu + util_cpusec()
         if (nodezero.and.util_print('t4',print_default)) then
            write(LuOut,*)
            write(LuOut,9000) 't4',size_t4
            write(LuOut,9010) 't4',filename(1:20)
            write(LuOut,9090) 't4',d_t4
            write(LuOut,9020) cpu
            call util_flush(LuOut)
         endif
         if (left) then
            call tce_y4_offset(l_l4_offset,k_l4_offset,size_l4)
            call tce_filename('lambda4',filename)
            call createfile(filename,d_lambda4,size_l4)
            call gatoeaf(d_lambda4)
         endif
      endif
c     
c     =======================
c     Create residual offsets
c     =======================
c     
      call tce_e_offset(l_e_offset,k_e_offset,size_e)
      if (needt1) call tce_t1_offset(l_r1_offset,k_r1_offset,size_r1)
      if (needt2) call tce_t2_offset(l_r2_offset,k_r2_offset,size_r2)
      if (needt3) call tce_t3_offset(l_r3_offset,k_r3_offset,size_r3)
      if (needt3a) call tce_t3a_offset(l_r3_offset,k_r3_offset,size_r3)
      if (needt4) call tce_t4_offset(l_r4_offset,k_r4_offset,size_r4)
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
         write(LuOut,9090) 'e2',d_e2
      endif
      if (needu1) then
         call tce_filename('u1',filename)
         call tce_t1_offset(l_u1_offset,k_u1_offset,size_u1)
         call createfile(filename,d_u1,size_u1)
         call tce_guess_t1(d_u1,k_u1_offset)
         call reconcilefile(d_u1,size_u1)
         write(LuOut,9090) 'u1',d_u1
      endif
      if (needu2) then
         call tce_filename('u2',filename)
         call tce_t2_offset(l_u2_offset,k_u2_offset,size_u2)
         call createfile(filename,d_u2,size_u2)
         call tce_guess_t2(d_v2,k_v2_offset,d_u2,k_u2_offset)
         call reconcilefile(d_u2,size_u2)
         write(LuOut,9090) 'u2',d_u2
      endif
      if (needw1) then
         call tce_filename('w1',filename)
         call tce_t1_offset(l_w1_offset,k_w1_offset,size_w1)
         call createfile(filename,d_w1,size_w1)
         call gatoeaf(d_w1)
         write(LuOut,9090) 'w1',d_w1
      endif
      if (needw2) then
         call tce_filename('w2',filename)
         call tce_t2_offset(l_w2_offset,k_w2_offset,size_w2)
         call createfile(filename,d_w2,size_w2)
         call gatoeaf(d_w2)
         write(LuOut,9090) 'w2',d_w2
      endif
      if (needw3) then
         call tce_filename('w3',filename)
         call tce_t3_offset(l_w3_offset,k_w3_offset,size_w3)
         call createfile(filename,d_w3,size_w3)
         call gatoeaf(d_w3)
         write(LuOut,9090) 'w3',d_w3
      endif
      if (needw4) then
         call tce_filename('w4',filename)
         call tce_t4_offset(l_w4_offset,k_w4_offset,size_w4)
         call createfile(filename,d_w4,size_w4)
         call gatoeaf(d_w4)
         write(LuOut,9090) 'w4',d_w4
      endif
c     
c     =========================
c     CC / CI / MBPT iterations
c     =========================
c     
      if (model.eq.'ccd') then
c     -------------
c     CCD     right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if (model.eq.'lccd') then
c     -------------
c     LCCD    right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if (model.eq.'ccsd') then
c     -------------
c     CCSD    right
c     -------------
         modelname = "CCSD"
c
c         call ccsd_energy(d_f1,d_e,d_t1,d_t2,d_v2,d_r1,d_r2,
c     1        k_f1_offset,k_e_offset,k_t1_offset,
c     2        k_t2_offset,k_v2_offset,k_r1_offset,
c     3        k_r2_offset,size_e,size_t1,size_t2,
c     4        size_r1,size_r2,ref,corr)
c
        call ccsd_energy_loc(d_e,d_f1,d_v2,d_t1,d_t2,
     1       k_e_offset,k_f1_offset,k_v2_offset,
     2       k_t1_offset,k_t2_offset,
     3       size_t1,size_t2,ref,corr)
c
         if(nroots.eq.0) then
         call ccsd_left(d_f1,d_lr1,d_lr2,d_t1,d_t2,d_v2,d_lambda1,
     1        d_lambda2,k_f1_offset,k_l1_offset,k_t1_offset,
     2        k_t2_offset,k_v2_offset,k_l2_offset,size_l1,size_l2,
     3        iter,nodezero,residual)
c
         call ccsd_gradient(d_f1,k_f1_offset,d_v2,k_v2_offset,d_t1,
     1           k_t1_offset,d_t2,k_t2_offset,d_lambda1,k_l1_offset,
     2           d_lambda2,k_l2_offset,
     3           rtdb)
         endif
      else if (model.eq.'lccsd') then
c     -------------
c     LCCSD   right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if (model.eq.'ccsdt') then
c     -------------
c     CCSDT   right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if (model.eq.'ccsdta') then 
c     -----------------------------
c     CCSDt (CCSDT-active)   right
c     -----------------------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if (model.eq.'ccsdtq') then
c     -------------
c     CCSDTQ  right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if (model.eq.'qcisd') then
c     -------------
c     QCISD   right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if (model.eq.'cisd') then
c     -------------
c     CISD    right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if (model.eq.'cisdt') then
c     -------------
c     CISDT   right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if (model.eq.'cisdtq') then
c     -------------
c     CISDTQ  right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      else if ((model.eq.'mbpt2').or.
     1         (model.eq.'mbpt3').or.
     2         (model.eq.'mbpt4')) then
c     -------------
c     MBPT(2) right
c     -------------
         call errquit("Not yet implemented",0,CAPMIS_ERR)
c
      endif
      if (nodezero) call util_flush(LuOut)
c     
c     Delete files for energy calculation
c     
      if (needw4) then
         call deletefile(d_w4)
         if (.not.ma_pop_stack(l_w4_offset))
     1        call errquit("tce_gradient: MA problem",1,MA_ERR)
      endif
      if (needw3) then
         call deletefile(d_w3)
         if (.not.ma_pop_stack(l_w3_offset))
     1        call errquit("tce_gradient: MA problem",2,MA_ERR)
      endif
      if (needw2) then
         call deletefile(d_w2)
         if (.not.ma_pop_stack(l_w2_offset))
     1        call errquit("tce_gradient: MA problem",3,MA_ERR)
      endif
      if (needw1) then
         call deletefile(d_w1)
         if (.not.ma_pop_stack(l_w1_offset))
     1        call errquit("tce_gradient: MA problem",4,MA_ERR)
      endif
      if (needu2) then
         call deletefile(d_u2)
         if (.not.ma_pop_stack(l_u2_offset))
     1        call errquit("tce_gradient: MA problem",5,MA_ERR)
      endif
      if (needu1) then
         call deletefile(d_u1)
         if (.not.ma_pop_stack(l_u1_offset))
     1        call errquit("tce_gradient: MA problem",6,MA_ERR)
      endif
      if (neede2) then
         call deletefile(d_e2)
         if (.not.ma_pop_stack(l_e2_offset))
     1        call errquit("tce_gradient: MA problem",7,MA_ERR)
      endif
      if (needt4) then
         if (.not.ma_pop_stack(l_r4_offset))
     1        call errquit("tce_gradient: MA problem",8,MA_ERR)
      endif
      if (needt3) then
         if (.not.ma_pop_stack(l_r3_offset))
     1        call errquit("tce_gradient: MA problem",9,MA_ERR)
      endif
      if (needt3a) then
         if (.not.ma_pop_stack(l_r3_offset))
     1        call errquit("tce_gradient: MA problem",10,MA_ERR)
      endif
      if (needt2) then
         if (.not.ma_pop_stack(l_r2_offset))
     1        call errquit("tce_gradient: MA problem",11,MA_ERR)
      endif
      if (needt1) then
         if (.not.ma_pop_stack(l_r1_offset))
     1        call errquit("tce_gradient: MA problem",12,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_e_offset))
     1     call errquit("tce_gradient: MA problem",13,MA_ERR)
      if(nroots.eq.0) then
         if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,ref+corr))
     1        call errquit('tce_gradient: RTDB problem',0,MA_ERR)
c     link to QM/MM
         if (.not.rtdb_put(rtdb,'tce:enegr',mt_dbl,1,ref+corr))
     1        call errquit('tce_gradient: RTDB problem',0,MA_ERR)
c     
         tce_gradient=.true.
      else
c     
c     ============================
c     Excited-state CI, CC, & MBPT
c     ============================
c     
         irrep_g = 0
c     
c     Restricted closed shell: irrep_g is always zero.
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
               if (.not.ma_push_get(mt_dbl,maxtrials,'omegax',
     1              l_omegax,k_omegax)) 
     2              call errquit('tce_gradient: MA problem',14,
     3              MA_ERR)
               if (left) then
                  if (.not.ma_push_get(mt_dbl,maxtrials,'omegay',
     1                 l_omegay,k_omegay)) 
     2                 call errquit('tce_gradient: MA problem',15,
     3                 MA_ERR)
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
               cpu=-util_cpusec()
               wall=-util_wallsec()
               if (model.eq.'ccsd') then
c     -------------
c     CCSD    right
c     -------------
c                  call errquit("Not yet implemented",0,CAPMIS_ERR)
c
                  modelname="EOM-CCSD"
                  call eomccsd_energy(d_rx1,d_rx2,size_x1,size_x2,
     1                 k_x1_offset,k_x2_offset,
     2                 d_f1,d_v2,d_t1,d_t2,
     3                 k_f1_offset,k_v2_offset,
     4                 k_t1_offset,k_t2_offset,k_omegax)
                  call eomccsd_left(d_ry1,d_ry2,size_y1,size_y2,
     1                 k_y1_offset,k_y2_offset,
     2                 d_f1,d_v2,d_t1,d_t2,
     3                 k_f1_offset,k_v2_offset,
     4                 k_t1_offset,k_t2_offset,k_omegay)
c
                  ivec=nroots_reduced
                  call dratoga(xc1(ivec))
                  call dratoga(xc2(ivec))
                  call dratoga(yc1(ivec))
                  call dratoga(yc2(ivec))
                  omegai=dbl_mb(k_omegax+ivec-1)
c                  call tce_print_x1(xc1(ivec),k_x1_offset,1.0d-1)
c                  call tce_print_x2(xc2(ivec),k_x2_offset,1.0d-1)
                  call eomccsd_gradients(rtdb,
     1                 d_f1,k_f1_offset,d_v2,k_v2_offset,
     1                 d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                 xc1(ivec),k_x1_offset,size_x1,
     1                 xc2(ivec),k_x2_offset,size_x2,
     1                 yc1(ivec),k_y1_offset,size_y1,
     1                 yc2(ivec),k_y2_offset,size_y2,
     1                 omegai)
                  call gatodra(yc2(ivec))
                  call gatodra(yc1(ivec))
                  call gatodra(xc2(ivec))
                  call gatodra(xc1(ivec))
c
               else if (model.eq.'ccsdt') then
c     -------------
c     CCSDT   right
c     -------------
                  call errquit("Not yet implemented",0,CAPMIS_ERR)
c
               else if (model.eq.'ccsdta') then
c     -----------------------
c     CCSDt (EOMCCSDt)  right
c     -----------------------
                  call errquit("Not yet implemented",0,CAPMIS_ERR)
c
               else if (model.eq.'ccsdtq') then
c     -------------
c     CCSDTQ  right
c     -------------
                  call errquit("Not yet implemented",0,CAPMIS_ERR)
c
               else
                  call errquit("Not yet implemented",0,CAPMIS_ERR)
               endif
c     
c     --------------------
c     Write energy to RTDB
c     --------------------
c     
               call tce_eom_tidy
               if ((targetsym.eq.irrepname)
     1              .and.(target.le.nroots_reduced)) then
                  if (nodezero.and.util_print('excited state',print_low)
     1               ) then
                     write(LuOut,*)
                     write(LuOut,9300) target
                     write(LuOut,9340) targetsym
                     write(LuOut,*)
                     write(LuOut,9310) 
     1                    modelname(1:inp_strlen(modelname)),
     1                    ref+corr
                     write(LuOut,9320) 
     1                    modelname(1:inp_strlen(modelname)),
     1                    dbl_mb(k_omegax+target-1)
                     write(LuOut,9330) 
     1                    modelname(1:inp_strlen(modelname)),
     1                    ref+corr+dbl_mb(k_omegax+target-1)
                  endiF
                  if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,
     1                 ref+corr+dbl_mb(k_omegax+target-1)))
     2                 call errquit('tce_gradient: RTDB problem',0,
     3                 RTDB_ERR)
                  tce_gradient=.true.
               else if ((targetsym.eq.'none').and.(irrep_x.eq.0)
     1                 .and.(target.le.nroots_reduced)) then
                  if (.not.rtdb_put(rtdb,'tce:omega',mt_dbl,1,
     1                 dbl_mb(k_omegax+target-1)*au2ev))
     2                 call errquit('tce_gradient: RTDB problem',0,
     3                 RTDB_ERR)
                  if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,
     1                 ref+corr+dbl_mb(k_omegax+target-1)))
     2                 call errquit('tce_gradient: RTDB problem',0,
     3                 RTDB_ERR)
                  tce_gradient=.true.
               endif
c     link to QM/MM (no symmetry OR a-symmetry) -------
               if(irrep.eq.0) then
                  if (.not.rtdb_put(rtdb,'tce:eneex',mt_dbl,
     1                 nroots_reduced,dbl_mb(k_omegax)))
     2                 call errquit('tce_gradient: RTDB problem',0,
     3                 RTDB_ERR)
               end if
c     
c     ------------------------------------
c     Terminate excited state calculations
c     ------------------------------------
c     
               if (needt4) then
                  if (left) then
                     call deletefile(d_ry4)
                     if (.not.ma_pop_stack(l_y4_offset))
     1                    call errquit("tce_gradient: MA problem",16,
     2                    MA_ERR)
                  endif
                  call deletefile(d_rx4)
                  if (.not.ma_pop_stack(l_x4_offset))
     1                 call errquit("tce_gradient: MA problem",17,
     2                 MA_ERR)
               endif
               if (needt3) then
                  if (left) then
                     call deletefile(d_ry3)
                     if (.not.ma_pop_stack(l_y3_offset))
     1                    call errquit("tce_gradient: MA problem",18,
     2                    MA_ERR)
                  endif
                  call deletefile(d_rx3)
                  if (.not.ma_pop_stack(l_x3_offset))
     1                 call errquit("tce_gradient: MA problem",19,
     2                 MA_ERR)
               endif
               if (needt3a) then
                  call deletefile(d_rx3)
                  if (.not.ma_pop_stack(l_x3_offset))
     1                 call errquit("tce_gradient: MA problem",20,
     2                 MA_ERR)
               endif
               if (needt2) then
                  if (left) then
                     call deletefile(d_ry2)
                     if (.not.ma_pop_stack(l_y2_offset))
     1                    call errquit("tce_gradient: MA problem",21,
     2                    MA_ERR)
                  endif
                  call deletefile(d_rx2)
                  if (.not.ma_pop_stack(l_x2_offset))
     1                 call errquit("tce_gradient: MA problem",22,
     2                 MA_ERR)
               endif
               if (needt1) then
                  if (left) then
                     call deletefile(d_ry1)
                     if (.not.ma_pop_stack(l_y1_offset))
     1                    call errquit("tce_gradient: MA problem",23,
     2                    MA_ERR)
                  endif
                  call deletefile(d_rx1)
                  if (.not.ma_pop_stack(l_x1_offset))
     1                 call errquit("tce_gradient: MA problem",24,
     2                 MA_ERR)
               endif
               if (left) then
                  if (.not.ma_pop_stack(l_y0_offset))
     1                 call errquit("tce_gradient: MA problem",25,
     2                 MA_ERR)
                  if (.not.ma_pop_stack(l_x0_offset))
     1                 call errquit("tce_gradient: MA problem",26,
     2                 MA_ERR)
                  if (.not.ma_pop_stack(l_d0_offset))
     1                 call errquit("tce_gradient: MA problem",27,
     2                 MA_ERR)
                  if (.not.ma_pop_stack(l_omegay))
     1                 call errquit("tce_gradient: MA problem",28,
     2                 MA_ERR)
               endif
               if (.not.ma_pop_stack(l_omegax))
     1              call errquit("tce_gradient: MA problem",29,
     2              MA_ERR)
            endif
         enddo
      endif
c     
c     ===========================
c     End of all TCE calculations
c     ===========================
c     
      call deletefile(d_v2)
      call deletefile(d_f1)
      if (needt4) then
         if (left) then
            call deletefile(d_lambda4)
            if (.not.ma_pop_stack(l_l4_offset))
     1           call errquit("tce_gradient: MA problem",30,MA_ERR)
         endif
         call deletefile(d_t4)
         if (.not.ma_pop_stack(l_t4_offset))
     1        call errquit("tce_gradient: MA problem",31,MA_ERR)
      endif
      if (needt3) then
         if (left) then
            call deletefile(d_lambda3)
            if (.not.ma_pop_stack(l_l3_offset))
     1           call errquit("tce_gradient: MA problem",32,MA_ERR)
         endif
         call deletefile(d_t3)
         if (.not.ma_pop_stack(l_t3_offset))
     1        call errquit("tce_gradient: MA problem",33,MA_ERR)
      endif
      if (needt3a) then
         call deletefile(d_t3)
         if (.not.ma_pop_stack(l_t3_offset))
     1        call errquit("tce_gradient: MA problem",34,MA_ERR)
      endif
      if (needt2) then
         if (left) then
            call deletefile(d_lambda2)
            if (.not.ma_pop_stack(l_l2_offset))
     1           call errquit("tce_gradient: MA problem",35,MA_ERR)
         endif
         call deletefile(d_t2)
         if (.not.ma_pop_stack(l_t2_offset))
     1        call errquit("tce_gradient: MA problem",36,MA_ERR)
      endif
      if (needt1) then
         if (left) then
            call deletefile(d_lambda1)
            if (.not.ma_pop_stack(l_l1_offset))
     1           call errquit("tce_gradient: MA problem",37,MA_ERR)
         endif
         call deletefile(d_t1)
         if (.not.ma_pop_stack(l_t1_offset))
     1        call errquit("tce_gradient: MA problem",38,MA_ERR)
      endif
      if(intorb) then
         if (.not.ma_pop_stack(l_v2spin_offset))
     1        call errquit("tce_gradient: MA problem",39,MA_ERR)
         if (.not.ma_pop_stack(l_v2_alpha_offset))
     1        call errquit("tce_gradient: MA problem",40,MA_ERR)
      end if
      if(.not.intorb) then
         if (.not.ma_pop_stack(l_v2_offset))
     1        call errquit("tce_gradient: MA problem",41,MA_ERR)
      end if
      if (.not.ma_pop_stack(l_f1_offset))
     1     call errquit("tce_gradient: MA problem",42,MA_ERR)
      if(intorb) then
         if (.not.ma_pop_stack(l_offset_alpha))
     1        call errquit("tce_gradient: MA problem",43,MA_ERR)
         if (.not.ma_pop_stack(l_range_alpha))
     1        call errquit("tce_gradient: MA problem",44,MA_ERR)
         if (.not.ma_pop_stack(l_sym_alpha))
     1        call errquit("tce_gradient: MA problem",45,MA_ERR)
         if (.not.ma_pop_stack(l_spin_alpha))
     1        call errquit("tce_gradient: MA problem",46,MA_ERR)
         if (.not.ma_pop_stack(l_b2am))
     1        call errquit("tce_gradient: MA problem",47,MA_ERR)
      end if
c     ---
      if (activecalc.or.(.not.intorb)) then
         if (.not.ma_pop_stack(l_active))
     1        call errquit("tce_gradient: MA problem",48,MA_ERR)
      end if
      if (.not.ma_pop_stack(l_alpha))
     1     call errquit("tce_gradient: MA problem",49,MA_ERR)
      if (.not.ma_pop_stack(l_offset))
     1     call errquit("tce_gradient: MA problem",50,MA_ERR)
      if (.not.ma_pop_stack(l_range))
     1     call errquit("tce_gradient: MA problem",51,MA_ERR)
      if (.not.ma_pop_stack(l_sym))
     1     call errquit("tce_gradient: MA problem",52,MA_ERR)
      if (.not.ma_pop_stack(l_spin))
     1     call errquit("tce_gradient: MA problem",53,MA_ERR)
      if(.not.ma_pop_stack(l_mo_index))
     1     call errquit("tce_gradient: ma problem",54,ma_err)
      if (.not.ma_pop_stack(l_evl_sorted))
     1     call errquit("tce_gradient: MA problem",55,MA_ERR)
      if (.not.ma_pop_stack(l_irs_sorted))
     1     call errquit("tce_gradient: MA problem",56,MA_ERR)
      if (.not.ma_pop_stack(l_spin_sorted))
     1     call errquit("tce_gradient: MA problem",57,MA_ERR)
      if (.not.ma_pop_stack(l_movecs_sorted))
     1     call errquit("tce_gradient: MA problem",58,MA_ERR)
c     
c     ===============
c     Destroy a mutex
c     ===============
c     
      if (.not.ga_destroy_mutexes())
     1     call errquit('tce_gradient: GA problem',1,GA_ERR)
c     
c     =========
c     Terminate
c     =========
c     
      call tce_grad_tidy(rtdb)
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
     1   1x,'--------------------------------------------------------',/
     2   1x,'Iter          Residuum       Correlation     Cpu    Wall',/
     3   1x,'--------------------------------------------------------')
 9060 format(
     1   1x,'--------------------------------------------------------',/
     2   1x,'Iterations converged')
 9070 format(1x,A,' correlation energy / hartree = ',f25.15)
 9080 format(1x,A,' total energy / hartree       = ',f25.15)
 9085 format(1x,A,' excitation energy (eV)       = ',f12.5)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9250 format(1x,'Ground-state symmetry is ',A4)
 9200 format(1x,'=========================================',/,
     1     1x,'Excited-state calculation ( ',A4,'symmetry)',/,
     2     1x,'=========================================')
 9210 format(/,1x,'Iteration ',i3,' using ',i4,' trial vectors')
 9220 format(/,1x,A,' iterations',/,1x,
     1  '--------------------------------------------------------------'
     2  ,/,1x,
     3  '     Residuum       Omega / hartree  Omega / eV    Cpu    Wall'
     4  ,/,1x,
     5 '--------------------------------------------------------------')
 9230 format(1x,f17.13,f18.13,f11.5,2f8.1)
 9240 format(1x,
     1  '--------------------------------------------------------------'
     2  ,/,1x,'Iterations converged')
 9300 format(' Target root     =',i3)
 9340 format(' Target symmetry =',1x,a4)
 9310 format(1x,A,' ground state energy / hartree  =',f25.15)
 9320 format(1x,A,' excitation energy / hartree    =',f25.15)
 9330 format(1x,A,' excited state energy / hartree =',f25.15)
 9400 format(/,1x,A,' iterations',/,
     1     1x,'--------------------------------------',/
     2     1x,'Iter          Residuum     Cpu    Wall',/
     3     1x,'--------------------------------------')
 9410 format(
     1     1x,'--------------------------------------',/
     2     1x,'Iterations converged')
 9420 format(1x,i4,f18.13,2f8.1)
 9430 format(/,1x,A,' dipole moments / hartree & Debye',/,
     1     1x,'--------------------------------',/
     2     1x,'X ',2f15.7,/
     3     1x,'Y ',2f15.7,/
     4     1x,'Z ',2f15.7,/
     5     1x,'--------------------------------')
 9440 format(1x,A1,' axis ( ',A4,'symmetry)')
 9450 format(1x,'(T) & [T] code has been provided by ',
     1     'A.A.Auer (Waterloo)')
 9460 format(/,1x,'Excited state root',i3,/
     1     1x,'Excitation energy / hartree =',f25.15,/
     2     1x,'                  / eV      =',f25.15)
 9470 format(/,1x,A,' transition moments / hartree',/
     1     1x,'--------------------------------------------',/
     2     1x,'X',f13.7,' Y',f13.7,' Z',f13.7,/
     3     1x,'Oscillator Strength            ',f13.7,/
     4     1x,'--------------------------------------------')
      end
c $Id$
