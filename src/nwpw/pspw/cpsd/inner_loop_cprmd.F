*
* $Id$
*
      subroutine inner_loop_cprmd(verlet,sa_alpha,
     >                      ispin,ne,neq,
     >                      nfft3d,nemaxq,
     >                      psi0,psi1,psi2,dn,
     >                      it_in,it_sum,E,
     >                      hml,lmd,
     >                      Hpsi,psi_r,
     >                      calc_pressure,pressure,p1,p2,
     >                      fractional,occ0,occ1,occ2)
      implicit none
      logical    verlet
      real*8     sa_alpha(2)
      integer    ispin,ne(2),neq(2)
      integer    nfft3d,nemaxq
      complex*16 psi0(nfft3d,nemaxq)
      complex*16 psi1(nfft3d,nemaxq)
      complex*16 psi2(nfft3d,nemaxq)
      real*8     dn(2*nfft3d,2)
      integer    it_in,it_sum
      real*8     E(*)
      real*8     hml(*),lmd(*)

*     **** very big workspace variables ****
      complex*16 Hpsi(nfft3d,nemaxq)
      real*8     psi_r(2*nfft3d,nemaxq)

      logical calc_pressure
      real*8  pressure,p1,p2,stress(3,3)

      logical fractional
      real*8 occ0(*),occ1(*),occ2(*)

#include "bafdecls.fh"
#include "errquit.fh"
ccccccc#include "frac_occ.fh"


*     **** local variables ****
      logical move,fei
      integer n2ft3d,np,np_i,np_j
      integer i,j,ii,jj,n,n1(2),n2(2),it,ms,nn,ierr,gga
      integer nx,ny,nz
      integer index,indext
      real*8  sum,Eold,eorbit,eion,ehartr,eke,eki,sse,ssr,sa1,sa2
      real*8  exc,exc2,pxc,pxc2,dte,dte0,scal1,scal2,dv,dt,fmass,h
      real*8  ehsic,phsic,exsic,pxsic,ehfx,phfx,espring,enlocal
      !real*8  e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      real*8   e_lj,e_q,e_spring
      real*8  s,r


*     **** MA local variables ****
      logical value,nose,field_exist,sic,allow_translation
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
      integer rho(2),vlr_l(2),r_grid(2)
      integer xcp(2),xce(2),dnall(2)
      integer v_field(2)
      integer natmx,fion(2),ftest(2)

*     ***** external functions ****
      integer  ion_nion,control_gga
      real*8   ion_ke,ion_ion_e,E_vnonlocal
      real*8   control_time_step,control_fake_mass,ion_dti
      real*8   lattice_omega,coulomb_e,ewald_e
      external ion_nion,control_gga
      external ion_ke,ion_ion_e,E_vnonlocal
      external control_time_step,control_fake_mass,ion_dti
      external lattice_omega,coulomb_e,ewald_e
      logical  psp_semicore
      external psp_semicore
      integer  control_version
      external control_version
      logical  control_Nose,control_Fei
      external control_Nose,control_Fei
      real*8   Nose_e_energy,Nose_r_energy,Nose_sse,Nose_ssr
      real*8   Nose_dXe,Nose_dXr
      external Nose_e_energy,Nose_r_energy,Nose_sse,Nose_ssr
      external Nose_dXe,Nose_dXr

*     ***** QM/MM external functions ****
      logical  pspw_charge_found,pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_charge_found,pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E

      logical  dplot_iteration_check
      external dplot_iteration_check

      logical  pspw_SIC,pspw_SIC_relaxed,control_allow_translation
      logical  pspw_HFX,pspw_HFX_relaxed
      external pspw_SIC,pspw_SIC_relaxed,control_allow_translation
      external pspw_HFX,pspw_HFX_relaxed

      double precision Dneall_m_trace
      external         Dneall_m_trace
      logical  Dneall_m_push_get_block,Dneall_m_pop_stack
      external Dneall_m_push_get_block,Dneall_m_pop_stack

      logical  meta_found,tamd_found,psp_U_psputerm
      external meta_found,tamd_found,psp_U_psputerm
      logical  nwpw_meta_gga_on,ion_disp_on
      external nwpw_meta_gga_on,ion_disp_on
      real*8   nwpw_meta_gga_pxc,ion_disp_energy
      external nwpw_meta_gga_pxc,ion_disp_energy



      call Parallel_np(np)

      n2ft3d = 2*nfft3d
      field_exist       = pspw_charge_found()
      sic               = pspw_SIC()
      gga               = control_gga()
      allow_translation = control_allow_translation()
      fei = control_Fei()

*     **** allocate MA local variables ****
      call nwpw_timing_start(12)
      value = Dneall_m_push_get_block(1,8,tmp_L)
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))

      if (control_version().eq.3) then
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'vc',  vc(2),  vc(1))
      end if

      if (field_exist) then
       value = value.and.
     >    BA_push_get(mt_dbl,(3*n2ft3d),'r_grid',r_grid(2),r_grid(1))
      end if

      value = value.and.
     >  BA_push_get(mt_dcpl,(n2ft3d),'v_field',v_field(2),v_field(1))

      value = value.and.
     >         BA_push_get(mt_dcpl,(nfft3d),'vl',  vl(2),  vl(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'rho',rho(2), rho(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'dng',dng(2), dng(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xcp',xcp(2), xcp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xce',xce(2), xce(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'dnall',dnall(2),dnall(1))
      natmx = ion_nion()
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'fion',fion(2),fion(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'ftest',ftest(2),ftest(1))
      if (.not.value) 
     > call errquit('inner_loop_md:pushing stack',0, MA_ERR)
      call dcopy(4*nfft3d,0.0d0,0,dbl_mb(dnall(1)),1)
      call dcopy(4*nfft3d,0.0d0,0,dbl_mb(xcp(1)),1)
      call dcopy(4*nfft3d,0.0d0,0,dbl_mb(xce(1)),1)

      call nwpw_timing_end(12)
     
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      move = .true.
      
      nose = control_Nose()
      sse = 1.0d0
      ssr = 1.0d0

      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1) + 1
      n2(2) = neq(1) + neq(2)

      dt    = control_time_step()
      fmass = control_fake_mass()
      dte   = dt*dt/fmass
      if (.not. verlet) dte=0.5d0*dte
      if (.not.nose) then
        sa1 =    1.0d0/(2.0d0-sa_alpha(1))
        sa2 = sa_alpha(1)/(2.0d0-sa_alpha(1))
      end if

      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()

      if (field_exist)
     >   call lattice_r_grid(dbl_mb(r_grid(1)))

      espring = 0.0d0


*     ******************************************
*     ****                                  ****
*     **** Start of molecular dynamics loop ****
*     ****                                  ****
*     ******************************************
      do it=1,it_in
        call dcopy(2*nfft3d*nemaxq,psi1,1,psi0,1)
        call dcopy(2*nfft3d*nemaxq,psi2,1,psi1,1)
*       *** skip ion_shift if newton step ***
        if (verlet) call ion_shift()
        if (nose.and.verlet) call Nose_shift()


*       ********************************
*       **** generate phaze factors ****
*       ********************************
        call phafac()
        if (control_version().eq.3) call ewald_phafac()

        call nwpw_timing_start(11)
*       *******************
*       **** get psi_r ****
*       *******************
        do n=n1(1),n2(ispin)
           call D3dB_c_Copy(1,psi1(1,n),psi_r(1,n))
        end do
        call Grsm_gh_fftb_nopack(nfft3d,n2(ispin),psi_r)
        do n=n1(1),n2(ispin)
           call D3dB_r_Zero_Ends(1,psi_r(1,n))
        end do

        call nwpw_meta_gga_gen_tau(ispin,neq,psi1)

*       *********************
*       **** generate dn ****
*       *********************
        call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
        if (fractional) then
          do ms=1,ispin
             do n=n1(ms),n2(ms)
                do i=1,n2ft3d
                   dn(i,ms) = dn(i,ms)
     >                      + scal2*(psi_r(i,n)**2)
     >                       *occ1(n)
                end do
             end do
             call D3dB_r_Zero_Ends(1,dn(1,ms))
             call D1dB_Vector_SumAll(n2ft3d,dn(1,ms))
          end do
        else
          do ms=1,ispin
             do n=n1(ms),n2(ms)
                do i=1,n2ft3d
                   dn(i,ms) = dn(i,ms) + scal2*(psi_r(i,n)**2)
                end do
             end do
             call D3dB_r_Zero_Ends(1,dn(1,ms))
             call D1dB_Vector_SumAll(n2ft3d,dn(1,ms))
          end do
        end if


*       **********************
*       **** generate dng ****
*       **********************
        call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
        call D3dB_r_SMul(1,scal1,dbl_mb(rho(1)),dcpl_mb(tmp1(1)))
        call D3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
c       ----**** need to to masking here ****
        call D3dB_c_Copy(0,dcpl_mb(tmp1(1)),dcpl_mb(dng(1)))

*       ********************************************************
*       **** generate dnall - used for semicore corrections ****
*       ********************************************************
        if (psp_semicore(0)) then
           call semicore_density_update()
           call semicore_density(dcpl_mb(tmp1(1)))
c           call D3dB_r_SMul(1,0.5d0,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)))
           call D3dB_r_SMul1(1,0.5d0,dcpl_mb(tmp1(1)))
        else
           call dcopy(n2ft3d,0.0d0,0,dcpl_mb(tmp1(1)),1)
        end if
        do ms=1,ispin
          call D3dB_rr_Sum(1,dn(1,ms),
     >                     dcpl_mb(tmp1(1)),
     >                     dbl_mb(dnall(1) +(ms-1)*n2ft3d))
        end do
        call nwpw_timing_end(11)





*       *****************************************
*       **** generate local pseudopotential  ****
*       **** and also get force if move true ****
*       *****************************************
        call v_local(dcpl_mb(vl(1)),
     >               move,
     >               dcpl_mb(dng(1)),
     >               dbl_mb(fion(1)))


*       ************************************
*       **** generate coulomb potential ****
*       ************************************
        if (control_version().eq.3)
     >     call coulomb_v_nopack(dcpl_mb(dng(1)),dcpl_mb(vc(1)))


*       *************************************************
*       **** generate exchange-correlation potential ****
*       *************************************************
        call v_bwexc_all_tmp1(gga,n2ft3d,ispin,
     >                        dbl_mb(dnall(1)),
     >                        dbl_mb(xcp(1)),
     >                        dbl_mb(xce(1)),
     >                        dcpl_mb(tmp1(1)))

     
*       **********************************************
*       **** generate other real-space potentials ****
*       **********************************************
        if (field_exist) then

           call dcopy(n2ft3d,0.0d0,0,dbl_mb(v_field(1)),1)


*          **** generate charge potential ****
           if (pspw_charge_found()) then
            call pspw_charge_Generate_V(n2ft3d,
     >                               dbl_mb(r_grid(1)),
     >                               dbl_mb(v_field(1)))
           end if

        end if


*       ******************
*       **** get Hpsi ****
*       ******************
        if (control_version().eq.3)
     >  call psi_H_nopack(ispin,neq,psi1,psi_r,
     >             dcpl_mb(vl(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,dbl_mb(fion(1)),fractional,occ1)


*       **********************
*       **** get ewald force *
*       **********************
*       **** get the ewald force ****
        if (control_version().eq.3) call ewald_f(dbl_mb(fion(1)))

*       ************************
*       **** get semicoreforce *
*       ************************
        if (psp_semicore_nopack(0)) then
           call semicore_xc_F_nopack(ispin,dbl_mb(xcp(1)),
     >                                     dbl_mb(fion(1)))
        end if

*       **** get the qmmm force ****
        if (pspw_qmmm_found()) call pspw_qmmm_fion(dbl_mb(fion(1)))

*       **** get the dispersion force ****
        if (ion_disp_on()) call ion_disp_force(dbl_mb(fion(1)))


*       *****************************************
*       **** remove ion forces using ion_FixIon *
*       *****************************************
        if (fei) 
     >    call dcopy(3*natmx,dbl_mb(fion(1)),1,dbl_mb(ftest(1)),1)
        call ion_FixIon(dbl_mb(fion(1)))
c        if (meta) call meta_force(ispin,neq,psi1,Hpsi,dbl_mb(fion(1)))

        !**** center of mass constraint ****
c        if (.not.allow_translation) then
c          call remove_center_F_mass(dbl_mb(fion(1)))
c        end if

*       **************************
*       **** do a verlet step ****
*       **************************
        if (verlet) then
*          **** constant temperature ****
           if (nose) then 
             sse = Nose_sse()
             ssr = Nose_ssr()
             do n=1,n2(ispin)
              call D3dB_c_SMul(1,0.5d0*dte,Hpsi(1,n),psi2(1,n))
              call D3dB_cc_daxpy(1,-1.0d0,psi0(1,n),psi2(1,n))
              call D3dB_cc_daxpy(1,1.0d0,psi1(1,n),psi2(1,n))
              call D3dB_c_SMul1(1,2.0d0*sse,psi2(1,n))
              call D3dB_cc_daxpy(1,1.0d0,psi0(1,n),psi2(1,n))
             end do
             call ion_nose_step(ssr,dbl_mb(fion(1)))

*          **** constant energy ****
           else
             do n=1,n2(ispin)
              call D3dB_c_SMul(1,dte*sa1,Hpsi(1,n),psi2(1,n))
              call D3dB_cc_daxpy(1,-1.0d0*sa2,psi0(1,n),psi2(1,n))
              call D3dB_cc_daxpy(1,2.0d0*sa1,psi1(1,n),psi2(1,n))
             end do

*            **** QM/MM Verlet update ****
             call ion_verlet_step(dbl_mb(fion(1)),sa_alpha(2))
           end if

*       **************************
*       **** do a newton step ****
*       **************************
        else
           r = 1.0d0
           s = 1.0d0
           if (nose) then
             r =  (1.0d0-0.5d0*dt*Nose_dXr())
             s =  (1.0d0-0.5d0*dt*Nose_dXe())
           end if
           do n=1,n2(ispin)
              call D3dB_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
              call D3dB_cc_daxpy(1,s*dt*sa_alpha(1),psi0(1,n),psi2(1,n))
              call D3dB_cc_Sum2(1,psi1(1,n),psi2(1,n))
           end do

*          **** QM/MM Newton update ****
           call ion_newton_step(dbl_mb(fion(1)),sa_alpha(2)*r)

        end if


*       *****************************************        
*       **** lagrange multiplier corrections ****
*       *****************************************        

        !**** orthoganality constraint ****
        dte0 = dte
        if (nose.and.verlet) dte0 = dte*sse
        if (fractional) then
        call psi_lmbda2_nopack(ispin,neq,nemaxq,nfft3d,psi1,psi2,
     >                  dte0,occ1,
     >                  lmd,dbl_mb(tmp_L(1)),ierr)
        else if (sic) then
        call psi_lmbda_sic_nopack(ispin,neq,nemaxq,nfft3d,psi1,psi2,dte0,
     >                 lmd,dbl_mb(tmp_L(1)),ierr)
        else
        call psi_lmbda_nopack(ispin,neq,nemaxq,nfft3d,psi1,psi2,dte0,
     >                 lmd,dbl_mb(tmp_L(1)),ierr)
        end if

        !**** center of mass constraint ****

        !**** total angular momentum constraint ****


*       **************************
*       *** update thermostats ***
*       **************************
        if (nose) then
          if (verlet) then
             eke = 0.0d0
             do i=1,n2(ispin)
                call D3dB_cc_idot(1,psi2(1,i),psi0(1,i),sum)
                eke = eke+sum
             end do
             if (np.gt.1) call Parallel_SumAll(eke)
             eke = (ne(1)+ne(2) - eke)
             if (ispin.eq.1) eke = 2.0d0*eke
             eke = 0.5d0*(fmass/(dt*dt))*eke
             eki = ion_ke()
             call Nose_Verlet_Step(eke,eki)
          else
              eke = 0.0d0
              do i=1,n2(ispin)
                call D3dB_cc_idot(1,psi0(1,i),psi0(1,i),sum)
                eke = eke+sum
              end do
              if (np.gt.1) call Parallel_SumAll(eke)
              if (ispin.eq.1) eke = 2.0d0*eke
              eke = eke*fmass
              eki = ion_ke()
              call Nose_Newton_Step(eke,eki)
          end if
        end if


*       ********************
*       **** call dplot ****
*       ********************
        if (dplot_iteration_check(it+it_sum)) then
         call dplot_iteration((it+it_sum),ispin,neq,psi1,dn,psi_r)
        end if


      end do
*     ******************************************************
*     ***** end main loop **********************************
*     ******************************************************

*     **** if newton then skip energy calculations ****
      if (.not. verlet) goto 333


*     *************************************
*     ***** total energy calculation ******
*     *************************************
      call nwpw_timing_start(10)
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

*     *** get orbital energies ****
      call Dneall_ffm_sym_Multiply_nopack(0,psi1,Hpsi,nfft3d,hml)
      call Dneall_m_scal(0,(-1.0d0),hml)
      if (fractional) then
         call Dneall_m_diag_scal(0,occ1,hml)
         eorbit = Dneall_m_trace(0,hml)
         call Dneall_m_diag_scal_inv(0,occ1,hml)
      else
         eorbit = Dneall_m_trace(0,hml)
      end if
      if (ispin.eq.1) eorbit = eorbit+eorbit


*     **** get ewald energy ****
      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()

           
*     **** get coulomb energy ****
      if (control_version().eq.3) ehartr = coulomb_e_nopack(dcpl_mb(dng(1)))


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
         pxc = pxc + nwpw_meta_gga_pxc(ispin,neq,psi1)
      end if



        
*     **** velocity and kinetic energy of psi ****
      h = 1.0d0/(2.0d0*dt)
      eke = 0.0d0
      do i=1,n2(ispin)
         call D3dB_c_SMul1(1,-h,psi0(1,i))
         call D3dB_cc_daxpy(1,h,psi2(1,i),psi0(1,i))
         call D3dB_cc_idot(1,psi0(1,i),psi0(1,i),sum)
         eke = eke+sum
      end do
      if (np.gt.1) call Parallel_SumAll(eke)
      eke = eke*fmass
      if (ispin.eq.1) eke = 2.0d0*eke

 

*     **** total energy ****
      Eold=E(1)
      E(2) = eorbit + eion + exc - ehartr - pxc + espring
      E(3) = eke
      E(4) = ion_ke()
      E(5) = eorbit
      E(6) = ehartr
      E(7) = exc
      E(8) = eion
      E(22) = espring

*     ******** QM/MM energies ******
      if (pspw_qmmm_found()) then
         e_lj     = pspw_qmmm_LJ_E()
         e_q      = pspw_qmmm_Q_E()
         e_spring = pspw_qmmm_spring_E()
         E(2)  = E(2) + e_lj + e_q + e_spring

         E(11) = e_lj
         E(12) = e_q
         E(13) = e_spring
      end if

*     **** SIC corrections ****
      if (pspw_SIC()) then
         call pspw_energy_SIC(ispin,psi_r,ehsic,phsic,exsic,pxsic)
         E(2) = E(2) + ehsic + exsic
         E(16) = ehsic
         E(17) = exsic
         if (pspw_SIC_relaxed()) then
            E(2)  = E(2) - phsic - pxsic
            E(18) = phsic
            E(19) = pxsic
         end if
      end if

*     **** HFX corrections ****
      if (pspw_HFX()) then
         call pspw_energy_HFX(ispin,psi_r,ehfx,phfx)
         E(2) = E(2) + ehfx
         E(20) = ehfx
         if (pspw_HFX_relaxed()) then
            E(2)  = E(2) - phfx
            E(21) = phfx
         end if
      end if

*     **** DFT+U terms ****
      if (psp_U_psputerm()) then
         call psp_U_psputerm_energy(ehfx,phfx)
         E(29) =  ehfx
         E(30) =  phfx
         E(2)  = E(2) + E(29) - E(30)
      end if

*     **** metadynamics energy ****
      if (meta_found()) then
         call meta_energypotential(ispin,neq,psi1,E(31),E(32))
         E(2)  = E(2) + E(31) - E(32)
      end if

*     **** dispersion energy ***
      if (ion_disp_on()) then
         E(33) = ion_disp_energy()
         E(2)  = E(2) + E(33)
      end if

*     **** tamd energy ****
      if (tamd_found()) then
         call tamd_energypotential(ispin,neq,psi1,E(34),E(35))
         E(2)  = E(2) + E(34) - E(35)
      end if


*     **** Energy and Energy**2 sum ***
      E(25) = E(25) + E(2)
      E(26) = E(26) + E(2)*E(2)
      E(27) = E(27) + E(2)+E(3)+E(4)
      E(28) = E(28) + (E(2)+E(3)+E(4))**2

*     **** output Forces for Fei ***
      if (fei) call fei_output(E(2),dbl_mb(ftest(1)))


*     **** Nose thermostat energies ****
      if (nose) then
        E(9)  = Nose_e_energy()
        E(10) = Nose_r_energy()
        E(1)  = E(2)+E(3)+E(4)+E(9)+E(10)
      else
        E(1) = E(2)+E(3)+E(4)
      end if


*     ******** pressure ******
      if (calc_pressure) then

        enlocal = E_vnonlocal(ispin,neq,fractional,occ1)

        call cgsd_pressure_stress(ispin,neq,psi1,
     >                            dbl_mb(dnall(1)),
     >                            dcpl_mb(dng(1)),
     >                            dbl_mb(xcp(1)),
     >                            enlocal,exc,pxc,
     >                            pressure,p1,p2,stress)
      end if


*      **** write ecce data ****
       call ecce_print_module_entry('task car-parrinello')
       call ion_ecce()

       call ecce_print1('total energy', mt_dbl,     E(1), 1)
       call ecce_print1('total kinetic', mt_dbl,    E(3)+E(4), 1)
       call ecce_print1('potential energy', mt_dbl, E(2), 1)
       call ecce_print1('electron kinetic', mt_dbl, E(3), 1)
       call ecce_print1('ion kinetic', mt_dbl,      E(4), 1)
       call ecce_print1('time', mt_dbl,      (it_in+it_sum)*dt, 1)

       call ecce_print2('total gradient', mt_dbl, dbl_mb(fion(1)),
     $        3,3,natmx)
       call ecce_print1('gradient norm', mt_dbl, E(1), 1)
       call ecce_print1('orbital gradient norm', mt_dbl, E(4), 1)
c       call ecce_print1('gradient max', mt_dbl, E(1), 1)
       call ecce_print_module_exit('task car-parrinello', 'ok')


      call nwpw_timing_end(10)

*     **** dealocate MA local variables ****
 333  continue
      call nwpw_timing_start(12)

      value = BA_pop_stack(ftest(2))
      value = value.and.BA_pop_stack(fion(2))
      value = value.and.BA_pop_stack(dnall(2))
      value = value.and.BA_pop_stack(xce(2))
      value = value.and.BA_pop_stack(xcp(2))
      value = value.and.BA_pop_stack(dng(2))
      value = value.and.BA_pop_stack(rho(2))
      value = value.and.BA_pop_stack(vl(2))
      value = value.and.BA_pop_stack(v_field(2))

      if (field_exist) 
     >  value = value.and.BA_pop_stack(r_grid(2))

      value = value.and.BA_pop_stack(vc(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and.Dneall_m_pop_stack(tmp_L)
      if (.not.value) 
     > call errquit('inner_loop_md:popping stack',1, MA_ERR)

      call nwpw_timing_end(12)

      return
      end
 
