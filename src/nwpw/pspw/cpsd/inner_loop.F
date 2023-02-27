* $Id$
*

*  ************************************************************
*  *                MPI cpsd routine                          *
*  *                                                          *
*  *  This is a developing cpsdv3 parallel code wrtten in     *
*  *  Fortran and MPI.                                        *
*  *                                                          *
*  *    + mpl message passing library used                    *
*  *                                                          *
*  *    + ngp is used instead of nfft in this proceudure      *
*  *                                                          *
*  *    + error checking is based on aimd.h parameters        *
*  *      then control file                                   *
*  ************************************************************
#ifdef USE_OPENMP
#define PSI_OMP
#endif

  
      subroutine inner_loop(ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi1,psi2,dn,
     >                      it_in,E,deltae,deltac,deltar,
     >                      hml,lmd,lmd1,first_iteration,
     >                      psi_r,Hpsi,
     >                      fractional,occ1,occ2)
#ifdef PSI_OMP
      use omp_lib
#endif
      implicit none
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax
      complex*16 psi1(npack1,nemax)
      complex*16 psi2(npack1,nemax)
      real*8     dn(2*nfft3d,2)
      integer    it_in
      real*8     E(*)
      real*8     deltae,deltac,deltar
      real*8     hml(*),lmd(*),lmd1(*)
      logical    first_iteration

*     **** very big workspace variables ****
      real*8     psi_r(2*nfft3d,nemax)
      complex*16 Hpsi(npack1,nemax)

      logical fractional
      real*8 occ1(*),occ2(*)


#include "bafdecls.fh"
#include "errquit.fh"
cccccccc#include "frac_occ.fh"

*     **** local variables ****
      logical move,sic,hfx
      integer n2ft3d,np_i,np_j,np
      integer i,n,n1(2),n2(2),it,ms,ierr
      integer nx,ny,nz
      integer gga
      real*8  sum,Eold,eorbit,eion,ehartr,eke,elocal,enlocal
      real*8  exc,exc2,pxc,pxc2,dte,scal1,scal2,dv,dt
      real*8  ehsic,phsic,exsic,pxsic,ehfx,phfx
      real*8  deltamm
      real*8  e_lj, e_q, e_spring,Eapc,Papc
      real*8  ehartree_atom,ecmp_cmp,ecmp_pw,exc_atom,pxc_atom
      !real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      real*8 dipole(3)



*     **** MA local variables ****
      logical value,field_exist,cosmo_on,cosmo1_on
      logical V_APC_on
*     real*8     tmp_L(8*nemax*nemax)
*     complex*16 tmp1(nfft3d)
*     complex*16 tmp2(nfft3d)
c     complex*16  vl(nfft3d)
c     complex*16  vc(nfft3d)
c     complex*16 dng(nfft3d)
c     real*8     xcp(2*nfft3d,2)
c     real*8     xce(2*nfft3d,2)
c     real*8     fion(3,natmx)
      integer tmp_L(2)
      integer tmp1(2),tmp2(2)
      integer vl(2),vc(2),dng(2)
      integer vlr_l(2),rho(2),r_grid(2)
      integer v_field(2)
      integer xcp(2),xce(2),dnall(2)
      integer natmx,fion(2),ftest(2)
      integer sumi(2)
      integer npack0
      integer nion1

*     ***** external functions ****
      logical  control_gram_schmidt
      external control_gram_schmidt
      logical  control_move,psp_semicore,psp_pawexist,control_Fei
      logical  pspw_SIC,pspw_SIC_relaxed
      logical  pspw_HFX,pspw_HFX_relaxed
      integer  ion_nion,control_gga
      real*8   control_time_step,control_fake_mass,ion_dti
      real*8   lattice_omega,coulomb_e,ewald_e
      external control_move,psp_semicore,psp_pawexist,control_Fei
      external pspw_SIC,pspw_SIC_relaxed
      external pspw_HFX,pspw_HFX_relaxed
      external ion_nion,control_gga
      external control_time_step,control_fake_mass,ion_dti
      external lattice_omega,coulomb_e,ewald_e
      integer  control_version
      external control_version
      real*8   ion_ion_e,E_vnonlocal
      external ion_ion_e,E_vnonlocal

      double precision Dneall_m_trace
      external         Dneall_m_trace
      logical  Dneall_m_push_get_block,Dneall_m_pop_stack
      external Dneall_m_push_get_block,Dneall_m_pop_stack

*     ***** QM/MM external functions ****
      logical  pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      real*8   pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix
      external pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found,pspw_bqext,pspw_Efield_found
      external pspw_charge_found,pspw_bqext,pspw_Efield_found
      integer  pspw_charge_nion,pspw_Efield_type
      external pspw_charge_nion,pspw_Efield_type
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge
      !real*8   electron_psi_v_field_ave,pspw_Efield_Energy_ion
      !external electron_psi_v_field_ave,pspw_Efield_Energy_ion
      real*8   pspw_Efield_Energy_ion
      external pspw_Efield_Energy_ion
      real*8   pspw_Efield_Energy_dipole
      external pspw_Efield_Energy_dipole
      real*8   dipole_Efield_e,dipole_Efield_p
      external dipole_Efield_e,dipole_Efield_p

      logical  psp_U_psputerm,meta_found
      external psp_U_psputerm,meta_found
      logical  nwpw_meta_gga_on,ion_disp_on
      external nwpw_meta_gga_on,ion_disp_on
      real*8   nwpw_meta_gga_pxc,ion_disp_energy
      external nwpw_meta_gga_pxc,ion_disp_energy


*     **** paw external functions ****
      real*8   psp_kinetic_atom,psp_kinetic_core,psp_valence_core_atom
      real*8   psp_ncmp_vloc,psp_vloc_atom,psp_ion_core
      real*8   psp_hartree_atom,psp_hartree_cmp_cmp,psp_hartree_cmp_pw
      external psp_kinetic_atom,psp_kinetic_core,psp_valence_core_atom
      external psp_ncmp_vloc,psp_vloc_atom,psp_ion_core
      external psp_hartree_atom,psp_hartree_cmp_cmp,psp_hartree_cmp_pw

*     **** cosmo external functions ****
      logical  nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      external nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      real*8   nwpw_cosmo_Eqq,nwpw_cosmo_EQionq,nwpw_cosmo_EQelcq
      external nwpw_cosmo_Eqq,nwpw_cosmo_EQionq,nwpw_cosmo_EQelcq


      integer tid,nthreads
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

#ifdef PSI_OMP
      integer nida,nidb
      integer  mthr

      real*8  adiff
      logical notgram

      integer MASTER,taskid
      parameter (MASTER=0)
      integer  Parallel_maxthreads
      external Parallel_maxthreads
      integer thrlmbda(2)

      INTEGER(kind=omp_nest_lock_kind) reduce_lock1
      INTEGER(kind=omp_nest_lock_kind) reduce_lock2
      INTEGER(kind=omp_nest_lock_kind) reduce_lock3
      common / reduce_ffm / reduce_lock1,reduce_lock2,reduce_lock3

*     **** matrix_blocking common block ****
      integer mblock(2),nblock(2),algorithm(2)
      common /matrix_blocking/ mblock,nblock,algorithm
#endif

      call nwpw_timing_start(42)


      call Pack_npack(0,npack0)
      n2ft3d = 2*nfft3d
      deltamm = 0.0d0

      field_exist = pspw_charge_found().or.pspw_Efield_found()
      sic        = pspw_SIC()
      hfx        = pspw_HFX()
      gga = control_gga()
      cosmo_on  = nwpw_cosmo_on()
      cosmo1_on = nwpw_cosmo1_on()
      V_APC_on  = pspw_V_APC_on()


      call nwpw_timing_start(12)
*     **** allocate MA local variables ****
      value = Dneall_m_push_get_block(1,8,tmp_L)
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))

      if (control_version().eq.3) then
       value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vc',  vc(2),  vc(1))
      end if

      if (control_version().eq.4) then
       value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'vc',vc(2),vc(1))
       value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'vlr_l',vlr_l(2),vlr_l(1))
      end if

      if ((field_exist).or.(control_version().eq.4)) then
       value = value.and.
     >    BA_push_get(mt_dbl,(3*n2ft3d),'r_grid',r_grid(2),r_grid(1))
      end if

       value = value.and.
     >  BA_push_get(mt_dbl,(n2ft3d),'v_field',v_field(2),v_field(1))

      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vloc', vl(2), vl(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'rho',rho(2),rho(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'dng',dng(2), dng(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xcp',xcp(2), xcp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xce',xce(2), xce(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'dnall',dnall(2),dnall(1))
      natmx = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) natmx = natmx + pspw_charge_nion()
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'fion',fion(2),fion(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'ftest',ftest(2),ftest(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemax),'sumi',sumi(2),sumi(1))

      call Parallel_shared_vector_zero(.false.,3*natmx,dbl_mb(fion(1)))
      call Parallel_shared_vector_zero(.true.,3*natmx,dbl_mb(ftest(1)))
#ifdef PSI_OMP
      mthr = Parallel_maxthreads()
      mthr = max(mthr,4)

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nida = 0
      if (taskid.eq.MASTER) nida = 1
      nidb = npack1-nida
      
      !write(*,*) mthr,ne(1)

c      value = value .and. MA_set_auto_verify(.true.)
      value = value.and.
     >        BA_push_get(mt_dbl,(mthr*ne(1)*ne(1)*8),'thrlmbda',
     >                                     thrlmbda(2),thrlmbda(1))

      call omp_init_nest_lock(reduce_lock1)
      call omp_init_nest_lock(reduce_lock2)
      call omp_init_nest_lock(reduce_lock3)
*     **** define blocking for dgemm matrix multiply ****
      algorithm(1) = -1
      algorithm(2) = -1
      do ms=1,ispin
         if (ne(ms).gt.(mthr*np)) then
            algorithm(ms) = 1
            call Parallel_matrixblocking(mthr*np,ne(ms),ne(ms),
     >                                mblock(ms),nblock(ms))
         else if (ne(ms).gt.(mthr)) then
            algorithm(ms) = 0
            call Parallel_matrixblocking(mthr,ne(ms),ne(ms),
     >                                mblock(ms),nblock(ms))
         else
            algorithm(ms) = -1
         end if
      end do

#endif

      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      call Parallel_shared_vector_zero(.false.,
     >                           4*nfft3d,dbl_mb(dnall(1)))
      call Parallel_shared_vector_zero(.false.,4*nfft3d,dbl_mb(xcp(1)))
      call Parallel_shared_vector_zero(.true.,4*nfft3d,dbl_mb(xce(1)))
      !call D3dB_r_nZero(1,2,dbl_mb(dnall(1)))
      !call D3dB_r_nZero(1,2,dbl_mb(xcp(1)))
      !call D3dB_r_nZero(1,2,dbl_mb(xce(1)))

      call nwpw_timing_end(12)
     
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      move = control_move()
      

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      dt = control_time_step()
      dte = dt/dsqrt(control_fake_mass())
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()

      if ((control_version().eq.4).or.(field_exist))
     >   call lattice_r_grid(dbl_mb(r_grid(1)))

      call pspw_Lin_HFX_reset()






*     ******************************************
*     ****                                  ****
*     ****   Start of steepest descent loop ****
*     ****                                  ****
*     ******************************************
!$OMP PARALLEL private(it,tid,nthreads)
      tid      = Parallel_threadid()
      nthreads = Parallel_nthreads()
!      write(*,*) 'num threads: ',nthreads
      do it=1,it_in


        call Pack_c_nCopy(1,nemax,psi2,psi1)

        if (fractional) then
           call Parallel_shared_vector_copy(.false.,ne(1)+ne(2),
     >                                      occ2,occ1)
        end if
     

        if (move) then
           call ion_shift()
           if (psp_pawexist()) call nwpw_gintegrals_set(.true.)
        end if

*       ********************************
*       **** generate phaze factors ****
*       ********************************

        
        call nwpw_timing_start(43)
c       these two functions are now threadsafe/multithreaded
        call phafac()
        if (control_version().eq.3) call ewald_phafac()
        call nwpw_timing_end(43)
        call nwpw_timing_start(11)


*       *******************;
*       **** get psi_r ****
*       *******************
        call nwpw_timing_start(44)

!$OMP DO private(n)
        do n=n1(1),n2(ispin)
           call Pack_c_Copy0(1,psi1(1,n),psi_r(1,n))
        end do
!$OMP END DO


        call Grsm_gh_fftb(nfft3d,n2(ispin),psi_r)

        call util_flush(6)
!$OMP BARRIER


        call nwpw_timing_end(44)

!$OMP DO private(n)
        do n=n1(1),n2(ispin)
           call D3dB_r_Zero_Ends0(1,psi_r(1,n))
        end do
!$OMP END DO


        call pspw_hfx_localize2_n(2,psi_r,psi1,psi2,psi2)

c       this function calls non pipelined 3d fft
        call nwpw_meta_gga_gen_tau(ispin,ne,psi1)


*       *********************
*       **** generate dn ****
*       *********************
c        call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
        call Parallel_shared_vector_zero(.true.,ispin*n2ft3d,dn)
        if (fractional) then
          do ms=1,ispin
             do n=n1(ms),n2(ms)
!$OMP DO private(i)
                do i=1,n2ft3d
                   dn(i,ms) = dn(i,ms) 
     >                      + scal2*(psi_r(i,n)**2)
     >                       *occ1(n)
                end do
!$OMP END DO
             end do
             !call D3dB_r_Zero_Ends(1,dn(1,ms))
             !call D1dB_Vector_SumAll(n2ft3d,dn(1,ms))
          end do
        else
          do ms=1,ispin
             do n=n1(ms),n2(ms)
!$OMP DO private(i)
                do i=1,n2ft3d
                   dn(i,ms) = dn(i,ms) + scal2*(psi_r(i,n)**2)
                end do
!$OMP END DO
             end do
             !call D3dB_r_Zero_Ends(1,dn(1,ms))
             !call D1dB_Vector_SumAll(n2ft3d,dn(1,ms))
          end do
        end if
        call D1dB_Vector_SumAll(ispin*n2ft3d,dn)

        call rho_symmetrizer_dn(ispin,n2ft3d,dn)



*       **********************
*       **** generate dng ****
*       **********************
        call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
        call D3dB_r_SMul(1,scal1,dbl_mb(rho(1)),dcpl_mb(tmp1(1)))
        call D3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
c        call mask_C(0,dcpl_mb(dng(1)))
        call Pack_c_pack(0,dcpl_mb(tmp1(1)))
        call Pack_c_Copy(0,dcpl_mb(tmp1(1)),dcpl_mb(dng(1)))
        if (cosmo1_on) then
           call v_local_set_cosmo_BQ(dbl_mb(r_grid(1)),
     >                               dbl_mb(rho(1)),
     >                               dcpl_mb(dng(1)))
        end if



*       ********************************************************
*       **** generate dnall - used for semicore corrections ****
*       ********************************************************
        if (psp_semicore(0)) then
           if (move .or. (it.eq.1)) call semicore_density_update()
           call semicore_density(dcpl_mb(tmp1(1)))
           call D3dB_r_SMul1(1,0.5d0,dcpl_mb(tmp1(1)))
        else
           call D3dB_r_Zero(1,dcpl_mb(tmp1(1)))
        end if
        do ms=1,ispin
          call D3dB_rr_Sum(1,dn(1,ms),
     >                     dcpl_mb(tmp1(1)),
     >                     dbl_mb(dnall(1) +(ms-1)*n2ft3d))
        end do
        
        call nwpw_timing_end(11)


!MATHIAS We might be able to add some parallelism between ions in that function 
        call v_local(dcpl_mb(vl(1)),
     >               move,
     >               dcpl_mb(dng(1)),
     >               dbl_mb(fion(1)))


*       *** long-range psp for charge systems ***
        if (control_version().eq.4) then
          call v_lr_local(dbl_mb(r_grid(1)),
     >                    dbl_mb(vlr_l(1)))
          if (move) then
             call grad_v_lr_local(dbl_mb(r_grid(1)),
     >                            dbl_mb(rho(1)),
     >                            dbl_mb(fion(1)))
          end if
        end if

        if (V_APC_on) then
           call pspw_V_APC(ispin,ne,
     >                           dcpl_mb(dng(1)),
     >                           dcpl_mb(vl(1)),
     >                           Eapc,Papc,
     >                           move,
     >                           dbl_mb(fion(1)))
        end if


*       ************************************
*       **** generate coulomb potential ****
*       ************************************
        if (control_version().eq.3) 
     >     call coulomb_v(dcpl_mb(dng(1)),dcpl_mb(vc(1)))

        if (control_version().eq.4) 
     >     call coulomb2_v(dbl_mb(rho(1)),dbl_mb(vc(1)))


        call nwpw_timing_start(46)
*       *************************************************
*       **** generate exchange-correlation potential ****
*       *************************************************
        call v_bwexc_all_tmp1(gga,n2ft3d,ispin,
     >                        dbl_mb(dnall(1)),
     >                        dbl_mb(xcp(1)),
     >                        dbl_mb(xce(1)),
     >                        dcpl_mb(tmp1(1)))

        call nwpw_timing_end(46)
     
*       *******************************************
*       ***** generate other real-space fields ****
*       *******************************************
        if (field_exist) then
           call D3dB_r_Zero(1,dbl_mb(v_field(1)))

*          **** generate charge potential ****
           call pspw_charge_Generate_V(n2ft3d,
     >                                    dbl_mb(r_grid(1)),
     >                                    dbl_mb(v_field(1)))
*          **** generate Efield potential ****
           call pspw_Efield_Generate_V(n2ft3d,
     >                                    dbl_mb(r_grid(1)),
     >                                    dbl_mb(v_field(1)))
        end if



        call nwpw_timing_start(47)


*       ******************
*       **** get Hpsi ****
*       ******************
        if (control_version().eq.3) 
     >     call psi_H(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,dbl_mb(fion(1)),fractional,occ1)

        if (control_version().eq.4) 
     >     call psi_Hv4(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),dbl_mb(vlr_l(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dbl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,dbl_mb(fion(1)),fractional,occ1)


        call nwpw_timing_end(47)

        
*       ************************************
*       **** do a steepest descent step ****
*       ************************************
        call nwpw_timing_start(45)
        do n=1,n2(ispin)
          call Pack_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
c          call Pack_cc_Sum(1,psi2(1,n),psi1(1,n),psi2(1,n))
          call Pack_cc_Sum2(1,psi1(1,n),psi2(1,n))
        end do


*       *******************************************
*       **** get ewald and semicore forces     ****
*       ****   and do steepest descent on ions ****
*       *******************************************
        if (move) then

*          **** get the ewald force ****
           if (control_version().eq.3) call ewald_f(dbl_mb(fion(1)))


*          **** get the free-space ion force ****
           if (control_version().eq.4) call ion_ion_f(dbl_mb(fion(1)))


*          **** get the semicore force ****
           if (psp_semicore(0)) then
              call semicore_xc_F(ispin,dbl_mb(xcp(1)),dbl_mb(fion(1)))
           end if

*          **** get the paw lagrange multiplier force ****
           if (psp_pawexist()) then
              if (first_iteration) then
                 call Dneall_ffm_sym_Multiply(0,psi1,Hpsi,npack1,hml)
                 call Dneall_m_scal(0,(-1.0d0),hml)
                 call Dneall_mm_copy(0,hml,lmd)
                 call Dneall_mm_copy(0,hml,lmd1)
                 first_iteration = .false.
              end if
              !tmp_L = 2*lmd - lmda1
              call Dneall_mm_copy(0,lmd,dbl_mb(tmp_L(1)))
              call Dneall_m_scal(0,2.0d0,dbl_mb(tmp_L(1)))
              call Dneall_mm_daxpy(0,-1.0d0,lmd1,dbl_mb(tmp_L(1)))
              call psp_paw_overlap_fion(ispin,
     >                             dbl_mb(tmp_L(1)),
     >                             psi1,
     >                             dbl_mb(fion(1)))
              call Dneall_mm_copy(0,lmd,lmd1)
           end if

*          **** get the qmmm force ****
           if (pspw_qmmm_found()) call pspw_qmmm_fion(dbl_mb(fion(1)))

*          **** get dispersion force ****
           if (ion_disp_on()) call ion_disp_force(dbl_mb(fion(1)))

*          **** get forces from external charges ****
           if (pspw_charge_found()) then
              if(pspw_bqext()) then
                 call pspw_charge_charge_Fion(dbl_mb(fion(1)))
              else
                 nion1 = ion_nion()
                 call pspw_charge_Fion_Fcharge(dbl_mb(fion(1)),
     >                                         dbl_mb(fion(1)+3*nion1))
                 call pspw_charge_Fcharge(dbl_mb(fion(1)+3*nion1))
                 call pspw_charge_rho_Fcharge(n2ft3d,dbl_mb(r_grid(1)),
     >                             dbl_mb(rho(1)),
     >                             dv,dbl_mb(fion(1)+3*nion1))
              end if
           end if

*          **** get forces from external Efield ****
           if (pspw_Efield_found()) then
              call pspw_Efield_Fion(dbl_mb(fion(1)))
           end if

*          **** remove ion forces using ion_FixIon ****
           call ion_FixIon(dbl_mb(fion(1)))

           call ion_optimize_step(dbl_mb(fion(1)))
        end if


*       *****************************************        
*       **** lagrange multiplier corrections ****
*       *****************************************        
        if (psp_pawexist()) then
           call psp_overlap_S(ispin,ne,psi1,psi_r)
           if (move) call phafac2()
           if (control_gram_schmidt()) then
              call Dneall_f_Sortho(ms,psi2,psi_r,npack1)
           else
              call psi_lmbda_paw(ispin,ne,nemax,npack1,psi_r,psi2,
     >                           dte,
     >                           lmd,dbl_mb(tmp_L(1)),ierr)
           end if
        else if (fractional) then
           call psi_lmbda2(ispin,ne,nemax,npack1,psi1,psi2,
     >                     dte,occ1,
     >                     lmd,dbl_mb(tmp_L(1)),ierr)
        !else if (sic.or.hfx) then
        else if (sic) then
           call psi_lmbda_sic(ispin,ne,nemax,npack1,psi1,psi2,dte,
     >                    lmd,dbl_mb(tmp_L(1)),ierr)
        else
#ifdef PSI_OMP
!           write(*,*) 'T',tid, ':', dbl_mb(thrlmbda(1))
           call psi_lmbda_omp(ispin,ne,nemax,nida,nidb,psi1,psi2,dte,
     >                  lmd,dbl_mb(tmp_L(1)),ierr,dbl_mb(thrlmbda(1)),
     >                                          taskid,adiff,notgram)
#else
           call psi_lmbda(ispin,ne,nemax,npack1,psi1,psi2,dte,
     >                    lmd,dbl_mb(tmp_L(1)),ierr)

#endif


        end if
        call nwpw_timing_end(45)

      end do

*     *************************************
*     ***** total energy calculation ******
*     *************************************
      call nwpw_timing_start(10)
      call Parallel_np(np)
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)
      if ((move).and.(psp_pawexist())) call phafac() !*** reset phase factors to r1 ***

*     *** get orbital energies ****
      call Dneall_ffm_sym_Multiply(0,psi1,Hpsi,npack1,hml)
      call Dneall_m_scal(0,(-1.0d0),hml)
!$OMP END PARALLEL 



      if (fractional) then
         call Dneall_m_diag_scal(0,occ1,hml)
         eorbit = Dneall_m_trace(0,hml)
         !eorbit = Dneall_m_occtrace(0,hml,occ1)
         call Dneall_m_diag_scal_inv(0,occ1,hml)
      else
         eorbit = Dneall_m_trace(0,hml)
      end if
      if (ispin.eq.1) eorbit = eorbit+eorbit
       

*     **** get ewald energy ****
      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()
    
*     **** get free-space ion-ion energy ****
      if (control_version().eq.4) eion = ion_ion_e()

           
*     **** get coulomb energy ****
      if (control_version().eq.3) ehartr = coulomb_e(dcpl_mb(dng(1)))
      if (control_version().eq.4) then
         call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(vc(1)),ehartr)
         ehartr = 0.5d0*ehartr*dv
      end if
        

*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dbl_mb(dnall(1)),dbl_mb(xce(1)),exc)
      call D3dB_rr_dot(1,dn(1,1),dbl_mb(xcp(1)),pxc)
      if (ispin.eq.1) then
         exc= exc + exc 
         pxc= pxc + pxc 
      else
         call D3dB_rr_dot(1,dbl_mb(dnall(1)+n2ft3d),
     >                      dbl_mb(xce(1)),exc2)
         call D3dB_rr_dot(1,dn(1,2),dbl_mb(xcp(1)+n2ft3d),pxc2)
         exc= exc + exc2
         pxc= pxc + pxc2
      end if
      exc = exc*dv
      pxc = pxc*dv

      if (nwpw_meta_gga_on()) then
         pxc = pxc + nwpw_meta_gga_pxc(ispin,ne,psi1)
      end if

        
  
*     ***** average Kohn-Sham kinetic energy ****
      call ke_ave(ispin,ne,psi1,eke,fractional,occ1)


*     **** average Kohn-Sham v_local energy ****
      call Pack_cc_dot(0,dcpl_mb(dng(1)),dcpl_mb(vl(1)),elocal)


*     *** add in long range part ****
      if (control_version().eq.4) then
      call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(vlr_l(1)),sum)
      elocal = elocal + sum*dv
      end if

*     *** add in other real-space fields ****
      if (field_exist) then
      call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(v_field(1)),sum)
      elocal = elocal + sum*dv
      end if

 
*     ***** average Kohn-Sham v_nonlocal energy ****
c      call dcopy(2*npack1*nemax,0.0d0,0,Hpsi,1)
c      call v_nonlocal(ispin,ne,psi1,Hpsi,
c     >                .false.,dbl_mb(ftest(1)),fractional,occ1)
c      enlocal = 0.0d0
c      do ms=1,ispin
c      do n=n1(ms),n2(ms)
c         call Pack_cc_idot(1,psi1(1,n),Hpsi(1,n),sum)
c         if (fractional) sum=sum*occ1(n)
c         enlocal = enlocal - sum
c      end do
c      end do
c      if (np.gt.1) call Parallel_SumAll(enlocal)
c      if (ispin.eq.1) enlocal = 2.0d0*enlocal

      enlocal = E_vnonlocal(ispin,ne,fractional,occ1)



*     **** total energy ****
!$OMP MASTER
      Eold=E(1)
      E(1) = eorbit + eion + exc - ehartr - pxc
      E(2) = eorbit
      E(3) = ehartr
      E(4) = exc
      E(5) = eion
      E(6) = eke
      E(7) = elocal
      E(8) = enlocal
      E(9) = 2.0d0*ehartr
      E(10) = pxc
!$OMP END MASTER

*     **** QM/MM terms ****
      if (pspw_qmmm_found()) then
         e_lj     = pspw_qmmm_LJ_E()
         e_q      = pspw_qmmm_Q_E()
         e_spring = pspw_qmmm_spring_E()
         E(1)  = E(1) + e_lj + e_q + e_spring

         E(11) = e_lj
         E(12) = e_q
         E(13) = e_spring

         E(14) = pspw_qmmm_LJ_Emix()
         E(14) = E(14) + pspw_qmmm_Q_Emix()
         call v_local_mm(dcpl_mb(vl(1)))
         call Pack_cc_dot(0,dcpl_mb(dng(1)),dcpl_mb(vl(1)),elocal)
         if (control_version().eq.4) then
            call v_lr_local_mm(dbl_mb(r_grid(1)),dbl_mb(vlr_l(1)))
            call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(vlr_l(1)),sum)
            elocal = elocal + sum*dv
         end if
         E(14) = E(14) + elocal
      end if

*     **** get pspw_charge and pspw_Efield energies ****
         !psi_1v_field()
      if (field_exist) then
         if (pspw_Efield_type().eq.0) then
c            call Calculate_Resta_Dipole(.false.,ispin,ne,ne,
c     >                                  npack1,nfft3d,psi1,dipole)
c            write(*,*) "dipole=",dipole  
c            write(*,*) " EP=",pspw_Efield_Energy_dipole(dipole) 

            !elocal = dipole_Efield_e()
            !e_q    = pspw_Efield_Energy_ion()
            E(49) = dipole_Efield_e()
            E(50) = dipole_Efield_p()
            !E(49) = elocal - e_q
            !E(50) = e_q
            E(51) = 0.0d0
            E(1) = E(1) + E(49) - E(50) !elocal - dipole_Efield_p()
         else
            !E(49)  = electron_psi_v_field_ave(psi1,dn)
            call D3dB_rr_dot(1,dn,dbl_mb(v_field(1)),sum)
            E(49) = sum*dv
            call D3dB_rr_dot(1,dn(1,ispin),dbl_mb(v_field(1)),sum)
            E(49) = E(49) + sum*dv

            E(50)  = pspw_charge_Energy_ion()
     >             + pspw_Efield_Energy_ion()
            E(51)  = pspw_charge_Energy_charge()
            E(1)   = E(1) + E(50) + E(51)
         end if
      end if



*     **** PAW ee terms ****
      if (psp_pawexist()) then
         E(35) = psp_kinetic_atom(ispin,ne,psi1)      !*** kinetic atom ***
         E(36) = psp_kinetic_core()                   !*** kinetic core  - independent of psi ***
         E(37) = psp_valence_core_atom(ispin,ne,psi1) !*** valence core energy ***

         E(38) = psp_ncmp_vloc(ispin)         !*** (ncmp)*Vlocal
         E(39) = psp_vloc_atom(ispin,ne,psi1) !*** Vlocal (local basis) -  nv*Vsrloc + nloc*Vloc - ncmp*Vloc ***

         ehartree_atom = psp_hartree_atom(ispin,ne,psi1)
         ecmp_cmp      = psp_hartree_cmp_cmp(ispin)
         ecmp_pw       = psp_hartree_cmp_pw(ispin,dcpl_mb(dng(1)),dn)
         E(40) = ehartree_atom             !*** vcoulomb atom  ***
         E(41) = ecmp_cmp                  !*** ncmp-ncmp coulomb energy ***
         E(42) = ecmp_pw                   !*** ncmp-pw coulomb energy ***

         call psp_xc_atom(ispin,ne,psi1,exc_atom,pxc_atom)
         E(43) = exc_atom                  !*** exc atom  ***
         E(44) = pxc_atom                  !*** pxc atom  ***

         E(45) = psp_ion_core()            !*** ion core energy - independent of psi ***


c         E(1) = E(1) + E(36) + E(45) 
c     >        + exc_atom - pxc_atom
c     >        - ehartree_atom - ecmp_cmp - ecmp_pw
         E(1) = E(1) 
     >        + exc_atom - pxc_atom
     >        - ehartree_atom - ecmp_cmp - ecmp_pw
      end if

*     **** COSMO terms ****
      if (cosmo_on) then

         !*** cosmo1 ****
         if (cosmo1_on) then
            call v_local_cosmo(dcpl_mb(vl(1)))
            call Pack_cc_dot(0,dcpl_mb(dng(1)),dcpl_mb(vl(1)),elocal)
            if (control_version().eq.4) then
               call v_lr_local_cosmo(dbl_mb(r_grid(1)),dbl_mb(vlr_l(1)))
               call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(vlr_l(1)),sum)
               elocal = elocal + sum*dv
            end if
            E(46) = elocal
            E(47) = nwpw_cosmo_EQionq()
            E(48) = nwpw_cosmo_Eqq()

         !*** cosmo2 ****
         else
            E(52) = Eapc
            E(53) = Papc
            E(46) = Eapc
            E(47) = nwpw_cosmo_EQionq() !** E(Qion|q)
            E(48) = nwpw_cosmo_Eqq()    !** E(q|q)

            !E(1)  = E(1) + E(52) - E(53) + E(47) + E(48)
            E(1)  = E(1) + E(52) - E(53) 
         end if
     
*     **** AQC Energy ****
      else if (V_APC_on) then
         E(52) = Eapc
         E(53) = Papc
         E(1)  = E(1) + E(52) - E(53)
      end if


*     **** SIC corrections ****
      if (pspw_SIC()) then
         call pspw_energy_SIC(ispin,psi_r,ehsic,phsic,exsic,pxsic)
         E(1) = E(1) + ehsic + exsic
         E(16) = ehsic
         E(17) = exsic
         if (pspw_SIC_relaxed()) then
            E(1)  = E(1) - phsic - pxsic
            E(18) = phsic
            E(19) = pxsic
         end if
      end if

*     **** HFX corrections ****
      if (pspw_HFX()) then
         call pspw_energy_HFX(ispin,psi_r,ehfx,phfx)
         E(1) = E(1) + ehfx
         E(20) = ehfx
         if (pspw_HFX_relaxed()) then
            E(1)  = E(1) - phfx 
            E(21) = phfx
         end if
      end if

*     **** DFT+U terms ****
      if (psp_U_psputerm()) then
         call psp_U_psputerm_energy(ehfx,phfx)
         E(29) =  ehfx
         E(30) =  phfx
         E(1)  = E(1) + E(29) - E(30)
      end if

*     **** metadynamics energy ****
      if (meta_found()) then
         call meta_energypotential(ispin,ne,psi1,E(31),E(32))
         E(1)  = E(1) + E(31) - E(32)
      end if

*     **** Dispersion energy ****
      if (ion_disp_on()) then
         E(33) = ion_disp_energy()
         E(1)  = E(1) + E(33)
      end if

*     **** output Forces for Fei ***
      if (control_Fei()) call fei_output(E(1),dbl_mb(fion(1)))


*     **** set convergence variables ****
      deltae = (E(1)-Eold)/(dt*dble(it_in))

*     *** deltac ***
      do n=n1(1),n2(ispin)
         do i=1,npack1
            Hpsi(i,n) = psi2(i,n) - psi1(i,n)
         end do
      end do

      do n=n1(1),n2(ispin)
         call Pack_cc_idot(1,Hpsi(1,n),Hpsi(1,n),dbl_mb(sumi(1)+n-1))
      end do
     
      if (np_i.gt.1) 
     >     call D3dB_Vector_SumAll((n2(ispin)-n1(1)+1),
     >                             dbl_mb(sumi(1)))

      deltac = 0.0d0
      do n=n1(1),n2(ispin)
         if (dbl_mb(sumi(1)+n-1).gt.deltac) deltac=dbl_mb(sumi(1)+n-1)
      end do
      if (np_j.gt.1) call D1dB_MaxAll(deltac)
      deltac = deltac/dte

c     do i=1,ion_nion()
c         fion(1,i) = (r2(1,i) - r1(1,i))/ion_dti(i)
c         fion(2,i) = (r2(2,i) - r1(2,i))/ion_dti(i)
c         fion(3,i) = (r2(3,i) - r1(3,i))/ion_dti(i)
c     end do

*     *** deltar  ***
      deltar = deltamm
      if (move) then
        do i=1,ion_nion()
           sum = dsqrt( dbl_mb(fion(1)+(i-1)*3  )**2 
     >                + dbl_mb(fion(1)+(i-1)*3+1)**2 
     >                + dbl_mb(fion(1)+(i-1)*3+2)**2)
           if (sum.gt.deltar) deltar = sum
        end do
      end if
      
       
      call nwpw_timing_end(10)

*     **** dealocate MA local variables ****
      call nwpw_timing_start(12)
#ifdef PSI_OMP
      value = BA_pop_stack(thrlmbda(2))

      call omp_destroy_nest_lock(reduce_lock1)
      call omp_destroy_nest_lock(reduce_lock2)
      call omp_destroy_nest_lock(reduce_lock3)
#endif
      value = BA_pop_stack(sumi(2))
      value = BA_pop_stack(ftest(2))
      value = BA_pop_stack(fion(2))
      value = BA_pop_stack(dnall(2))
      value = BA_pop_stack(xce(2))
      value = BA_pop_stack(xcp(2))
      value = BA_pop_stack(dng(2))
      value = BA_pop_stack(rho(2))
      value = BA_pop_stack(vl(2))
      value = BA_pop_stack(v_field(2))

      if ((control_version().eq.4).or.(field_exist))
     >   value = BA_pop_stack(r_grid(2))

      if (control_version().eq.4) 
     >   value = BA_pop_stack(vlr_l(2))

      value = BA_pop_stack(vc(2))
      value = BA_pop_stack(tmp2(2))
      value = BA_pop_stack(tmp1(2))
      value = Dneall_m_pop_stack(tmp_L)

      call nwpw_timing_end(12)

      call nwpw_timing_end(42)

      return
      end
 

c      subroutine nwpw_baddn_dum(tag,ispin,ne,npack1,nfft3d,psi1)
c      implicit none
c      character*(*) tag
c      integer ispin,ne(2),npack1,nfft3d
c      complex*16 psi1(npack1,*)
c
c      integer MASTER,taskid
c      parameter (MASTER=0)
c      logical value,oprint
c      integer psi_r(2),dn(2)
c      integer ms,n,i,n1(2),n2(2),nx,ny,nz,n2ft3d
c      real*8   exc(2),scal1,dv,scal2
c
c      real*8   lattice_omega
c      external lattice_omega
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c      call Parallel_taskid(taskid)
c      oprint = (taskid.eq.MASTER)
c
c      call D3dB_n2ft3d(1,n2ft3d)
c
c      n1(1) = 1
c      n2(1) = ne(1)
c      n1(2) = ne(1) + 1
c      n2(2) = ne(1) + ne(2)
c
c
c      call D3dB_nx(1,nx)
c      call D3dB_ny(1,ny)
c      call D3dB_nz(1,nz)
c
c      scal1 = 1.0d0/dble(nx*ny*nz)
c      scal2 = 1.0d0/lattice_omega()
c      dv    = scal1*lattice_omega()
c
c      value = BA_push_get(mt_dbl,(ne(1)+ne(2))*n2ft3d,'psi_rrr',
c     >                    psi_r(2),psi_r(1))
c      value = BA_push_get(mt_dbl,2*n2ft3d,'dnnn',dn(2),dn(1))
c      if (.not.value) write(*,*) "Crap"
c
c      call dcopy((ne(1)+ne(2))*n2ft3d,0.0d0,0,dbl_mb(psi_r(1)),1)
c      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)
c
c
c        do n=n1(1),n2(ispin)
c           call Pack_c_Copy(1,psi1(1,n),dbl_mb(psi_r(1)+(n-1)*n2ft3d))
c        end do
c
c        call Grsm_gh_fftb(nfft3d,n2(ispin),dbl_mb(psi_r(1)))
c
c        do n=n1(1),n2(ispin)
c           call D3dB_r_Zero_Ends(1,dbl_mb(psi_r(1)+(n-1)*n2ft3d))
c        end do
c
c         call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)
c         !call D3dB_r_nZero(1,ispin,dbl_mb(dn(1)))
c         do ms=1,ispin
c             do n=n1(ms),n2(ms)
c                do i=1,n2ft3d
c                   dbl_mb(dn(1)+i-1+(ms-1)*n2ft3d) 
c     >           = dbl_mb(dn(1)+i-1+(ms-1)*n2ft3d)
c     >              + scal2*((dbl_mb(psi_r(1)+i-1+(n-1)*n2ft3d))**2)
c                end do
c             end do
c             call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)+(ms-1)*n2ft3d))
c             call D1dB_Vector_SumAll(n2ft3d,dbl_mb(dn(1)+(ms-1)*n2ft3d))
c          end do
c
c        exc(1) = 0.0d0
c        exc(2) = 0.0d0
c        do ms=1,ispin
c           do i=1,n2ft3d
c             exc(ms) = exc(ms) + dbl_mb(dn(1)+i-1+(ms-1)*n2ft3d)
c           end do
c           call D3dB_SumAll(exc(ms))
c        end do
c        if (oprint) write(*,*) tag//" DNALL=",
c     >              exc(1)*dv,exc(2)*dv,exc,dv,
c     >              ispin,ne,n2ft3d,npack1
c
c      value = BA_pop_stack(dn(2))
c      value = BA_pop_stack(psi_r(2))
c
c
c      return
c      end

