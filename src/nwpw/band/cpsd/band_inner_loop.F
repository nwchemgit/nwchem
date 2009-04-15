*
* $Id: inner_loop.F,v 1.33 2009-02-07 03:50:36 bylaska Exp $
*
*  ************************************************************
*  *                                                          *
*  *                   BAND cpsd routine                      *
*  *                                                          *
*  ************************************************************
  
      subroutine band_inner_loop(ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi1_tag,psi2_tag,dn,
     >                      it_in,E,deltae,deltac,deltar,
     >                      hml_tag,lmd_tag,
     >                      psir_tag,Hpsi_tag)
      implicit none
      integer ispin,ne(2)
      integer npack1,nfft3d,nemax
      integer psi1_tag,psi2_tag
      real*8  dn(nfft3d,2)
      integer it_in
      real*8  E(*)
      real*8  deltae,deltac,deltar
      integer hml_tag,lmd_tag

*     **** very big workspace variables ****
      integer psir_tag,Hpsi_tag

#include "mafdecls.fh"
#include "errquit.fh"

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
      real*8  e_lj, e_q, e_spring
      !real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav



*     **** MA local variables ****
      logical value,field_exist
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

*     ***** external functions ****
      logical  control_move,psp_semicore
      logical  pspw_SIC,pspw_SIC_relaxed
      logical  pspw_HFX,pspw_HFX_relaxed
      integer  ion_nion,control_gga
      real*8   control_time_step,control_fake_mass,ion_dti
      real*8   lattice_omega,coulomb_e,ewald_e
      external control_move,psp_semicore
      external pspw_SIC,pspw_SIC_relaxed
      external pspw_HFX,pspw_HFX_relaxed
      external ion_nion,control_gga
      external control_time_step,control_fake_mass,ion_dti
      external lattice_omega,coulomb_e,ewald_e
      integer  control_version
      external control_version
      real*8   ion_ion_e
      external ion_ion_e

      double precision Dneall_m_trace
      external         Dneall_m_trace
      logical  Dneall_m_push_get_block,Dneall_m_pop_stack
      external Dneall_m_push_get_block,Dneall_m_pop_stack

*     ***** QM/MM external functions ****
      logical  pspw_charge_found,pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_charge_found,pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E

c      call Pack_npack(0,npack0)
c      n2ft3d = 2*nfft3d
c      deltamm = 0.0d0
c
c      field_exist = pspw_charge_found()
c      sic        = pspw_SIC()
c      hfx        = pspw_HFX()
c      gga = control_gga()
c
c      call nwpw_timing_start(12)
c*     **** allocate MA local variables ****
c      value = Dneall_m_push_get_block(1,8,tmp_L)
c      value = value.and.
c     >        MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
c      value = value.and.
c     >        MA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
c
c      if (control_version().eq.3) then
c       value = value.and.
c     >        MA_push_get(mt_dcpl,(npack0),'vc',  vc(2),  vc(1))
c      end if
c
c      if (control_version().eq.4) then
c       value = value.and.
c     >        MA_push_get(mt_dbl,(n2ft3d),'vc',vc(2),vc(1))
c       value = value.and.
c     >        MA_push_get(mt_dbl,(n2ft3d),'vlr_l',vlr_l(2),vlr_l(1))
c      end if
c
c      if ((field_exist).or.(control_version().eq.4)) then
c       value = value.and.
c     >    MA_push_get(mt_dbl,(3*n2ft3d),'r_grid',r_grid(2),r_grid(1))
c      end if
c
c       value = value.and.
c     >  MA_push_get(mt_dbl,(n2ft3d),'v_field',v_field(2),v_field(1))
c
c      value = value.and.
c     >        MA_push_get(mt_dcpl,(npack0),'vloc', vl(2), vl(1))
c      value = value.and.
c     >        MA_push_get(mt_dbl,(n2ft3d),'rho',rho(2),rho(1))
c      value = value.and.
c     >        MA_push_get(mt_dcpl,(npack0),'dng',dng(2), dng(1))
c      value = value.and.
c     >        MA_push_get(mt_dbl,(4*nfft3d),'xcp',xcp(2), xcp(1))
c      value = value.and.
c     >        MA_push_get(mt_dbl,(4*nfft3d),'xce',xce(2), xce(1))
c      value = value.and.
c     >        MA_push_get(mt_dbl,(4*nfft3d),'dnall',dnall(2),dnall(1))
c      natmx = ion_nion()
c      value = value.and.
c     >        MA_push_get(mt_dbl,(3*natmx),'fion',fion(2),fion(1))
c      value = value.and.
c     >        MA_push_get(mt_dbl,(3*natmx),'ftest',ftest(2),ftest(1))
c      value = value.and.
c     >        MA_push_get(mt_dbl,(nemax),'sumi',sumi(2),sumi(1))
c      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
c      call dcopy(4*nfft3d,0.0d0,0,dbl_mb(xcp(1)),1)
c      call dcopy(4*nfft3d,0.0d0,0,dbl_mb(xce(1)),1)
c
c      call nwpw_timing_end(12)
c     
c      call D3dB_nx(1,nx)
c      call D3dB_ny(1,ny)
c      call D3dB_nz(1,nz)
c      move = control_move()
c      
c
c      n1(1) = 1
c      n2(1) = ne(1)
c      n1(2) = ne(1) + 1
c      n2(2) = ne(1) + ne(2)
c
c      dt = control_time_step()
c      dte = dt/dsqrt(control_fake_mass())
c      scal1 = 1.0d0/dble(nx*ny*nz)
c      scal2 = 1.0d0/lattice_omega()
c      dv    = scal1*lattice_omega()
c
c      if ((control_version().eq.4).or.(field_exist))
c     >   call lattice_r_grid(dbl_mb(r_grid(1)))
c
c*     ******************************************
c*     ****                                  ****
c*     ****   Start of steepest descent loop ****
c*     ****                                  ****
c*     ******************************************
c      do it=1,it_in
c        call dcopy(2*npack1*nemax,psi2,1,psi1,1)
c        if (move) call ion_shift()
c
c*       ********************************
c*       **** generate phaze factors ****
c*       ********************************
c        call phafac()
c        if (control_version().eq.3) call ewald_phafac()
c        call nwpw_timing_start(11)
c*       *******************;
c*       **** get psi_r ****
c*       *******************
cc        do n=n1(1),n2(ispin)
cc           call Pack_c_Copy(1,psi1(1,n),psi_r(1,n))
cc           call Pack_c_unpack(1,psi_r(1,n))
cc           call D3dB_cr_fft3b(1,psi_r(1,n))
cc           call D3dB_r_Zero_Ends(1,psi_r(1,n))
cc        end do
c        do n=n1(1),n2(ispin)
c           call Pack_c_Copy(1,psi1(1,n),psi_r(1,n))
c        end do
c        call Grsm_gh_fftb(nfft3d,n2(ispin),psi_r)
c        do n=n1(1),n2(ispin)
c           call D3dB_r_Zero_Ends(1,psi_r(1,n))
c        end do
c
c
c
c*       *********************
c*       **** generate dn ****
c*       *********************
c        call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
c        if (fractional) then
c          do ms=1,ispin
c             do n=n1(ms),n2(ms)
c                do i=1,n2ft3d
c                   dn(i,ms) = dn(i,ms) 
c     >                      + scal2*(psi_r(i,n)**2)
c     >                       *occ1(n)
c                end do
c             end do
c             call D3dB_r_Zero_Ends(1,dn(1,ms))
c             call D1dB_Vector_SumAll(n2ft3d,dn(1,ms))
c          end do
c        else
c          do ms=1,ispin
c             do n=n1(ms),n2(ms)
c                do i=1,n2ft3d
c                   dn(i,ms) = dn(i,ms) + scal2*(psi_r(i,n)**2)
c                end do
c             end do
c             call D3dB_r_Zero_Ends(1,dn(1,ms))
c             call D1dB_Vector_SumAll(n2ft3d,dn(1,ms))
c          end do
c        end if
c
c*       **********************
c*       **** generate dng ****
c*       **********************
c        call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
c        call D3dB_r_SMul(1,scal1,dbl_mb(rho(1)),dcpl_mb(tmp1(1)))
c        call D3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
cc        call mask_C(0,dcpl_mb(dng(1)))
c        call Pack_c_pack(0,dcpl_mb(tmp1(1)))
c        call Pack_c_Copy(0,dcpl_mb(tmp1(1)),dcpl_mb(dng(1)))
c
c
c
c*       ********************************************************
c*       **** generate dnall - used for semicore corrections ****
c*       ********************************************************
c        if (psp_semicore(0)) then
c           if (move .or. (it.eq.1)) call semicore_density_update()
c           call semicore_density(dcpl_mb(tmp1(1)))
cc           call D3dB_r_SMul(1,0.5d0,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)))
c           call D3dB_r_SMul1(1,0.5d0,dcpl_mb(tmp1(1)))
c        else
c           call dcopy(n2ft3d,0.0d0,0,dcpl_mb(tmp1(1)),1)
c        end if
c        do ms=1,ispin
c          call D3dB_rr_Sum(1,dn(1,ms),
c     >                     dcpl_mb(tmp1(1)),
c     >                     dbl_mb(dnall(1) +(ms-1)*n2ft3d))
c        end do
c        
c        call nwpw_timing_end(11)
c
c
c*       *****************************************
c*       **** generate local pseudopotential  ****
c*       **** and also get force if move true ****
c*       *****************************************
c 
c        call v_local(dcpl_mb(vl(1)),
c     >               move,
c     >               dcpl_mb(dng(1)),
c     >               dbl_mb(fion(1)))
c
c
c*       *** long-range psp for charge systems ***
c        if (control_version().eq.4) then
c          call v_lr_local(dbl_mb(r_grid(1)),
c     >                    dbl_mb(vlr_l(1)))
c          if (move) then
c             call grad_v_lr_local(dbl_mb(r_grid(1)),
c     >                            dbl_mb(rho(1)),
c     >                            dbl_mb(fion(1)))
c          end if
c        end if
c
c
c*       ************************************
c*       **** generate coulomb potential ****
c*       ************************************
c        if (control_version().eq.3) 
c     >     call coulomb_v(dcpl_mb(dng(1)),dcpl_mb(vc(1)))
c
c        if (control_version().eq.4) 
c     >     call coulomb2_v(dbl_mb(rho(1)),dbl_mb(vc(1)))
c
c*       *************************************************
c*       **** generate exchange-correlation potential ****
c*       *************************************************
c        if (gga.eq.0) then
c          call vxc(n2ft3d,ispin,dbl_mb(dnall(1)),
c     >                      dbl_mb(xcp(1)),
c     >                      dbl_mb(xce(1)),
c     >                      dcpl_mb(tmp1(1)))
c        else if ((gga.ge.10).and.(gga.lt.100)) then
c          call v_bwexc(gga,n2ft3d,ispin,dbl_mb(dnall(1)),
c     >                      1.0d0,1.0d0,
c     >                      dbl_mb(xcp(1)),
c     >                      dbl_mb(xce(1)))
c
c        else if (gga.eq.110) then
c          call v_bwexc(10,n2ft3d,ispin,dbl_mb(dnall(1)),
c     >                      0.75d0,1.0d0,
c     >                      dbl_mb(xcp(1)),
c     >                      dbl_mb(xce(1)))
c        else if (gga.eq.112) then
c          call v_bwexc(12,n2ft3d,ispin,dbl_mb(dnall(1)),
c     >                      0.75d0,1.0d0,
c     >                      dbl_mb(xcp(1)),
c     >                      dbl_mb(xce(1)))
c
c        end if
c     
c*       *******************************************
c*       ***** generate other real-space fields ****
c*       *******************************************
c        if (field_exist) then
c
c           call dcopy(n2ft3d,0.0d0,0,dbl_mb(v_field(1)),1)
c
c
c*          **** generate charge potential ****
c           if (pspw_charge_found()) then
c              call pspw_charge_Generate_V(n2ft3d,
c     >                                    dbl_mb(r_grid(1)),
c     >                                    dbl_mb(v_field(1)))
c           end if
c        end if
c
c        call nwpw_timing_start(13)
c
c
c
c*       ******************
c*       **** get Hpsi ****
c*       ******************
c        if (control_version().eq.3) 
c     >     call psi_H(ispin,ne,psi1,psi_r,
c     >             dcpl_mb(vl(1)),
c     >             dbl_mb(v_field(1)),field_exist,
c     >             dcpl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
c     >             move,dbl_mb(fion(1)),fractional,occ1)
c
c        if (control_version().eq.4) 
c     >     call psi_Hv4(ispin,ne,psi1,psi_r,
c     >             dcpl_mb(vl(1)),dbl_mb(vlr_l(1)),
c     >             dbl_mb(v_field(1)),field_exist,
c     >             dbl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
c     >             move,dbl_mb(fion(1)),fractional,occ1)
c
c        
c*       ************************************
c*       **** do a steepest descent step ****
c*       ************************************
c        do n=1,n2(ispin)
c          call Pack_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
cc          call Pack_cc_Sum(1,psi2(1,n),psi1(1,n),psi2(1,n))
c          call Pack_cc_Sum2(1,psi1(1,n),psi2(1,n))
c        end do
c        call nwpw_timing_end(13)
c
c*       *******************************************
c*       **** get ewald and semicore forces     ****
c*       ****   and do steepest descent on ions ****
c*       *******************************************
c        if (move) then
c
c*          **** get the ewald force ****
c           if (control_version().eq.3) call ewald_f(dbl_mb(fion(1)))
c
c*          **** get the free-space ion force ****
c           if (control_version().eq.4) call ion_ion_f(dbl_mb(fion(1)))
c
c*          **** get the semicore force ****
c           if (psp_semicore(0)) then
c              call semicore_xc_F(ispin,dbl_mb(xcp(1)),dbl_mb(fion(1)))
c           end if
c
c*          **** get the qmmm force ****
c           if (pspw_qmmm_found()) call pspw_qmmm_fion(dbl_mb(fion(1)))
c
c*          **** remove ion forces using ion_FixIon ****
c           call ion_FixIon(dbl_mb(fion(1)))
c
c           call ion_optimize_step(dbl_mb(fion(1)))
c        end if
c
c
c*       *****************************************        
c*       **** lagrange multiplier corrections ****
c*       *****************************************        
c        if (fractional) then
c        call psi_lmbda2(ispin,ne,nemax,npack1,psi1,psi2,
c     >                  dte,occ1,
c     >                  lmd,dbl_mb(tmp_L(1)),ierr)
c        else if (sic.or.hfx) then
c        call psi_lmbda_sic(ispin,ne,nemax,npack1,psi1,psi2,dte,
c     >                 lmd,dbl_mb(tmp_L(1)),ierr)
c        else
c        call psi_lmbda(ispin,ne,nemax,npack1,psi1,psi2,dte,
c     >                 lmd,dbl_mb(tmp_L(1)),ierr)
c        end if
c
c      end do
c
c*     *************************************
c*     ***** total energy calculation ******
c*     *************************************
c      call nwpw_timing_start(10)
c      call Parallel_np(np)
c      call Parallel2d_np_i(np_i)
c      call Parallel2d_np_j(np_j)
c
c*     *** get orbital energies ****
c      call Dneall_ffm_sym_Multiply(0,psi1,Hpsi,npack1,hml)
c      call Dneall_m_scal(0,(-1.0d0),hml)
c      if (fractional) then
c         !call Dneall_m_diag_scal(0,dbl_mb(fweight(1)),hml)
c         eorbit = Dneall_m_trace(0,hml)
c         !call Dneall_m_diag_scal_inv(0,dbl_mb(fweight(1)),hml)
c      else
c         eorbit = Dneall_m_trace(0,hml)
c      end if
c      if (ispin.eq.1) eorbit = eorbit+eorbit
c       
c
c*     **** get ewald energy ****
c      eion = 0.0d0
c      if (control_version().eq.3) eion = ewald_e()
c    
c*     **** get free-space ion-ion energy ****
c      if (control_version().eq.4) eion = ion_ion_e()
c
c           
c*     **** get coulomb energy ****
c      if (control_version().eq.3) ehartr = coulomb_e(dcpl_mb(dng(1)))
c      if (control_version().eq.4) then
c         call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(vc(1)),ehartr)
c         ehartr = 0.5d0*ehartr*dv
c      end if
c        
c
c*     **** get exchange-correlation energy ****
c      call D3dB_rr_dot(1,dbl_mb(dnall(1)),dbl_mb(xce(1)),exc)
c      call D3dB_rr_dot(1,dn(1,1),dbl_mb(xcp(1)),pxc)
c      if (ispin.eq.1) then
c         exc= exc + exc 
c         pxc= pxc + pxc 
c      else
c         call D3dB_rr_dot(1,dbl_mb(dnall(1)+n2ft3d),
c     >                      dbl_mb(xce(1)),exc2)
c         call D3dB_rr_dot(1,dn(1,2),dbl_mb(xcp(1)+n2ft3d),pxc2)
c         exc= exc + exc2
c         pxc= pxc + pxc2
c      end if
c      exc = exc*dv
c      pxc = pxc*dv
c        
c  
c*     ***** average Kohn-Sham kinetic energy ****
c      call ke_ave(ispin,ne,psi1,eke,fractional,occ1)
c
c
c*     **** average Kohn-Sham v_local energy ****
c      call Pack_cc_dot(0,dcpl_mb(dng(1)),dcpl_mb(vl(1)),elocal)
c
c
c*     *** add in long range part ****
c      if (control_version().eq.4) then
c      call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(vlr_l(1)),sum)
c      elocal = elocal + sum*dv
c      end if
c
c*     *** add in other real-space fields ****
c      if (field_exist) then
c      call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(v_field(1)),sum)
c      elocal = elocal + sum*dv
c      end if
c
c 
c*     ***** average Kohn-Sham v_nonlocal energy ****
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
c
c
c
c*     **** total energy ****
c      Eold=E(1)
c      E(1) = eorbit + eion + exc - ehartr - pxc
c      E(2) = eorbit
c      E(3) = ehartr
c      E(4) = exc
c      E(5) = eion
c      E(6) = eke
c      E(7) = elocal
c      E(8) = enlocal
c      E(9) = 2.0d0*ehartr
c      E(10) = pxc
c
c      if (pspw_qmmm_found()) then
c         e_lj     = pspw_qmmm_LJ_E()
c         e_q      = pspw_qmmm_Q_E()
c         e_spring = pspw_qmmm_spring_E()
c         E(1)  = E(1) + e_lj + e_q + e_spring
c
c         E(11) = e_lj
c         E(12) = e_q
c         E(13) = e_spring
c      end if
c
c*     **** SIC corrections ****
c      if (pspw_SIC()) then
c         call pspw_energy_SIC(ispin,psi_r,ehsic,phsic,exsic,pxsic)
c         E(1) = E(1) + ehsic + exsic
c         E(16) = ehsic
c         E(17) = exsic
c         if (pspw_SIC_relaxed()) then
c            E(1)  = E(1) - phsic - pxsic
c            E(18) = phsic
c            E(19) = pxsic
c         end if
c      end if
c
c*     **** HFX corrections ****
c      if (pspw_HFX()) then
c         call pspw_energy_HFX(ispin,psi_r,ehfx,phfx)
c         E(1) = E(1) + ehfx
c         E(20) = ehfx
c         if (pspw_HFX_relaxed()) then
c            E(1)  = E(1) - phfx 
c            E(21) = phfx
c         end if
c      end if
c
c
c*     **** set convergence variables ****
c      deltae = (E(1)-Eold)/(dt*dble(it_in))
c
c*     *** deltac ***
c      do n=n1(1),n2(ispin)
c         do i=1,npack1
c            Hpsi(i,n) = psi2(i,n) - psi1(i,n)
c         end do
c      end do
c
c      do n=n1(1),n2(ispin)
c         call Pack_cc_idot(1,Hpsi(1,n),Hpsi(1,n),dbl_mb(sumi(1)+n-1))
c      end do
c     
c      if (np_i.gt.1) 
c     >     call D3dB_Vector_SumAll((n2(ispin)-n1(1)+1),
c     >                             dbl_mb(sumi(1)))
c
c      deltac = 0.0d0
c      do n=n1(1),n2(ispin)
c         if (dbl_mb(sumi(1)+n-1).gt.deltac) deltac=dbl_mb(sumi(1)+n-1)
c      end do
c      if (np_j.gt.1) call D1dB_MaxAll(deltac)
c      deltac = deltac/dte
c
cc     do i=1,ion_nion()
cc         fion(1,i) = (r2(1,i) - r1(1,i))/ion_dti(i)
cc         fion(2,i) = (r2(2,i) - r1(2,i))/ion_dti(i)
cc         fion(3,i) = (r2(3,i) - r1(3,i))/ion_dti(i)
cc     end do
c
c*     *** deltar ***
c      deltar = deltamm
c      if (move) then
c        do i=1,ion_nion()
c           sum = dsqrt( dbl_mb(fion(1)+(i-1)*3  )**2 
c     >                + dbl_mb(fion(1)+(i-1)*3+1)**2 
c     >                + dbl_mb(fion(1)+(i-1)*3+2)**2)
c           if (sum.gt.deltar) deltar = sum
c        end do
c      end if
c       
c      call nwpw_timing_end(10)
c
c*     **** dealocate MA local variables ****
c      call nwpw_timing_start(12)
c      value = MA_pop_stack(sumi(2))
c      value = MA_pop_stack(ftest(2))
c      value = MA_pop_stack(fion(2))
c      value = MA_pop_stack(dnall(2))
c      value = MA_pop_stack(xce(2))
c      value = MA_pop_stack(xcp(2))
c      value = MA_pop_stack(dng(2))
c      value = MA_pop_stack(rho(2))
c      value = MA_pop_stack(vl(2))
c      value = MA_pop_stack(v_field(2))
c
c      if ((control_version().eq.4).or.(field_exist))
c     >   value = MA_pop_stack(r_grid(2))
c
c      if (control_version().eq.4) 
c     >   value = MA_pop_stack(vlr_l(2))
c
c      value = MA_pop_stack(vc(2))
c      value = MA_pop_stack(tmp2(2))
c      value = MA_pop_stack(tmp1(2))
c      value = Dneall_m_pop_stack(tmp_L)
c
c      call nwpw_timing_end(12)

      return
      end
 
