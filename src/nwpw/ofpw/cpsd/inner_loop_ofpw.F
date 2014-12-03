*,eke
* $Id$
*

*  ************************************************************
*  *                cpsd routine                              *
*  *                                                          *
*  ************************************************************
  
      subroutine inner_loop_ofpw(ispin,ne,
     >                      npack1,nfft3d,
     >                      psi1,psi2,dn,
     >                      it_in,E,deltae,deltac,deltar,
     >                      hml,lmd,
     >                      psi_r,Hpsi)
      implicit none
      integer    ispin,ne(2)
      integer    npack1,nfft3d
      complex*16 psi1(npack1,ispin)
      complex*16 psi2(npack1,ispin)
      real*8     dn(2*nfft3d,2)
      integer    it_in
      real*8     E(*)
      real*8     deltae,deltac,deltar
      real*8     hml(*),lmd(*)

*     **** very big workspace variables ****
      real*8     psi_r(2*nfft3d,ispin)
      complex*16 Hpsi(npack1,ispin)


#include "bafdecls.fh"
#include "errquit.fh"


*     **** local variables ****
      logical move,sic,hfx
      integer n2ft3d,np_i,np_j
      integer i,n,n1(2),n2(2),it,ms,nn,ierr,nne(2)
      integer nx,ny,nz
      integer gga
      real*8  sum,Eold,eorbit,eion,ehartr,eke,elocal,enlocal
      real*8  A,B,C,x
      real*8  exc,exc2,pxc,pxc2,dte,scal1,scal2,dv,dt
      real*8  etf,etf2,ptf,ptf2
      real*8  evw,ewgc,pwgc,pwgc2
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
      integer v_field(2),v_tmp(2)
      integer xcp(2),xce(2),dnall(2),v_tf(2),e_tf(2),v_wgc(2)
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

      double precision Dneall_m_trace,wgc_e
      external         Dneall_m_trace,wgc_e
      logical  Dneall_m_push_get_block,Dneall_m_pop_stack
      external Dneall_m_push_get_block,Dneall_m_pop_stack

*     ***** QM/MM external functions ****
      logical  pspw_charge_found,pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_charge_found,pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E

      call Pack_npack(0,npack0)
      call D3dB_n2ft3d(1,n2ft3d)
      deltamm = 0.0d0

      field_exist = pspw_charge_found()
      gga = control_gga()

      call nwpw_timing_start(12)


*     **** allocate MA local variables ****
      value = BA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
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
       call dcopy(2*npack0,0.0d0,0,dcpl_mb(vl(1)),1)
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'rho',rho(2),rho(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'dng',dng(2), dng(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xcp',xcp(2), xcp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xce',xce(2), xce(1))

      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'v_tmp',v_tmp(2), v_tmp(1))

      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'v_tf',v_tf(2), v_tf(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nfft3d),'e_tf',e_tf(2), e_tf(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'v_wgc',v_wgc(2), v_wgc(1))
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(v_wgc(1)),1)

      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'dnall',dnall(2),dnall(1))
      natmx = ion_nion()
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'fion',fion(2),fion(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'ftest',ftest(2),ftest(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(ispin),'sumi',sumi(2),sumi(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      call nwpw_timing_end(12)
     
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      move = control_move()
      

      dt = control_time_step()
      dte = dt/dsqrt(control_fake_mass())
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()

      if ((control_version().eq.4).or.(field_exist))
     >   call lattice_r_grid(dbl_mb(r_grid(1)))

*     ******************************************
*     ****                                  ****
*     ****   Start of steepest descent loop ****
*     ****                                  ****
*     ******************************************
      do it=1,it_in
        call dcopy(2*npack1*ispin,psi2,1,psi1,1)
        if (move) call ion_shift()

*       ********************************
*       **** generate phaze factors ****
*       ********************************
        call phafac()
        if (control_version().eq.3) call ewald_phafac()


        call nwpw_timing_start(11)
*       *******************
*       **** get psi_r ****
*       *******************
        do n=1,ispin
           call Pack_c_Copy(1,psi1(1,n),psi_r(1,n))
           call Pack_c_unpack(1,psi_r(1,n))
           call D3dB_cr_pfft3b(1,1,psi_r(1,n))
        end do




*       *********************
*       **** generate dn ****
*       *********************
        call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
        do ms=1,ispin
           do i=1,n2ft3d
              dn(i,ms) = dn(i,ms) + dble(ne(ms))*scal2*(psi_r(i,ms)**2)
           end do
           call D3dB_r_Zero_Ends(1,dn(1,ms))
        end do

*       **********************
*       **** generate dng ****
*       **********************
        call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
        call D3dB_r_SMul(1,scal1,dbl_mb(rho(1)),dcpl_mb(tmp1(1)))
        !call D3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
        call D3dB_rc_pfft3f(1,0,dcpl_mb(tmp1(1)))
        call Pack_c_pack(0,dcpl_mb(tmp1(1)))
        call Pack_c_Copy(0,dcpl_mb(tmp1(1)),dcpl_mb(dng(1)))




*       ********************************************************
*       **** generate dnall - used for semicore corrections ****
*       ********************************************************
        if (psp_semicore(0)) then
           if (move .or. (it.eq.1)) call semicore_density_update()
           call semicore_density(dcpl_mb(tmp1(1)))
           call D3dB_r_SMul(1,0.5d0,dcpl_mb(tmp1(1)),dcpl_mb(tmp1(1)))
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


*       ************************************
*       **** generate coulomb potential ****
*       ************************************
        if (control_version().eq.3) 
     >     call coulomb_v(dcpl_mb(dng(1)),dcpl_mb(vc(1)))

        if (control_version().eq.4) 
     >     call coulomb2_v(dbl_mb(rho(1)),dbl_mb(vc(1)))



*       *************************************************
*       **** generate exchange-correlation potential ****
*       *************************************************
        if (gga.eq.0) then
          call vxc(n2ft3d,ispin,dbl_mb(dnall(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)),
     >                      dcpl_mb(tmp1(1)))
        else if ((gga.ge.10).and.(gga.lt.100)) then
          call v_bwexc(gga,n2ft3d,ispin,dbl_mb(dnall(1)),
     >                      1.0d0,1.0d0,
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)))

        else if (gga.eq.110) then
          call v_bwexc(10,n2ft3d,ispin,dbl_mb(dnall(1)),
     >                      0.75d0,1.0d0,
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)))
        else if (gga.eq.112) then
          call v_bwexc(12,n2ft3d,ispin,dbl_mb(dnall(1)),
     >                      0.75d0,1.0d0,
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)))

        end if

*       *************************************************
*       **** generate Kinetic  potentials            ****
*       *************************************************
       call v_thomasfermi(n2ft3d,ispin,
     >                    dn,
     >                    dbl_mb(v_tf(1)),
     >                    dbl_mb(e_tf(1)))
       call wgc_v(ispin,dn,dbl_mb(v_wgc(1)))

     
*       *******************************************
*       ***** generate other real-space fields ****
*       *******************************************
        if (field_exist) then

           call dcopy(n2ft3d,0.0d0,0,dbl_mb(v_field(1)),1)


*          **** generate charge potential ****
           if (pspw_charge_found()) then
              call pspw_charge_Generate_V(n2ft3d,
     >                                    dbl_mb(r_grid(1)),
     >                                    dbl_mb(v_field(1)))
           end if
        end if

        call nwpw_timing_start(13)


*       ******************
*       **** get Hchi ****
*       ******************
        if (control_version().eq.3) 
     >     call chi_H(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(vc(1)),dbl_mb(xcp(1)),
     >             dbl_mb(v_tf(1)),dbl_mb(v_wgc(1)),
     >             Hpsi,
     >             move,dbl_mb(fion(1)))

        if (control_version().eq.4) 
     >     call chi_Hv4(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),dbl_mb(vlr_l(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dbl_mb(vc(1)),dbl_mb(xcp(1)),
     >             dbl_mb(v_tf(1)),dbl_mb(v_wgc(1)),
     >             Hpsi,
     >             move,dbl_mb(fion(1)))

        
*       ************************************
*       **** do a steepest descent step ****
*       ************************************
        do n=1,ispin
          call Pack_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
          call Pack_cc_Sum2(1,psi1(1,n),psi2(1,n))
        end do
        call nwpw_timing_end(13)

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

*          **** get the qmmm force ****
           if (pspw_qmmm_found()) call pspw_qmmm_fion(dbl_mb(fion(1)))

*          **** remove ion forces using ion_FixIon ****
           call ion_FixIon(dbl_mb(fion(1)))

           call ion_optimize_step(dbl_mb(fion(1)))
        end if


*       *****************************************        
*       **** lagrange multiplier corrections ****
*       *****************************************        
        call chi_lmbda(ispin,npack1,psi1,psi2,dte,lmd,ierr)


      end do

*     *************************************
*     ***** total energy calculation ******
*     *************************************
      call nwpw_timing_start(10)
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

*     *** get orbital energies ****
      call Pack_cc_dot(1,psi1(1,1),Hpsi(1,1),eorbit)
      hml(1) = -eorbit
      if (ispin.eq.2) then
         call Pack_cc_dot(1,psi1(1,2),Hpsi(1,2),sum)
         hml(2) = -sum
         eorbit = eorbit + sum
      else
         eorbit = eorbit + eorbit
      end if
      eorbit = -eorbit


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
      !ehartr = 0.0d0
        

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
        
  
*     ***** Thomas-Fermi kinetic energy ****
      call D3dB_rr_dot(1,dbl_mb(dnall(1)),dbl_mb(e_tf(1)),etf)
      call D3dB_rr_dot(1,dn(1,1),dbl_mb(v_tf(1)),ptf)
      if (ispin.eq.1) then
         etf= etf + etf 
         ptf= ptf + ptf 
      else
         call D3dB_rr_dot(1,dbl_mb(dnall(1)+n2ft3d),
     >                      dbl_mb(e_tf(1)),etf2)
         call D3dB_rr_dot(1,dn(1,2),dbl_mb(v_tf(1)+n2ft3d),ptf2)
         etf= etf + etf2
         ptf= ptf + ptf2
      end if
      etf = etf*dv
      ptf = ptf*dv


*     ***** von-Weizsacker kinetic energy ****
      nne(1) = 1
      nne(2) = 0
      if (ispin.eq.1) then
         call ke_ave(ispin,nne,psi1,evw,.false.,n1)
         evw = evw*ne(1)
      else
         call ke_ave(ispin,nne,psi1,evw,.false.,n1)
         evw = evw*ne(1)
         call ke_ave(ispin,nne,psi1(1,2),ewgc,.false.,n1)
         evw = evw + ewgc*ne(2)
      end if
      


*     ***** WGC kinetic energy ****
      call D3dB_rr_dot(1,dn(1,1),dbl_mb(v_wgc(1)),pwgc)
      if (ispin.eq.1) then
         pwgc = pwgc + pwgc
      else
         call D3dB_rr_dot(1,dn(1,2),dbl_mb(v_wgc(1)),pwgc2)
         pwgc = pwgc + pwgc2
      end if
      pwgc = pwgc*dv
      ewgc = wgc_e(ispin,dn)
      !ewgc = 0.0d0
      !pwgc = 0.0d0



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
      enlocal = 0.0d0


*     **** total energy ****
      Eold=E(1)
      !E(1) = elocal + eion + exc + etf + ehartr 
      E(1) = eorbit + eion +exc+etf+ewgc-ehartr-pxc-ptf-pwgc
      E(2) = eorbit
      E(3) = ehartr
      E(4) = exc
      E(5) = eion
      E(6) = evw
      E(7) = elocal
      E(8) = enlocal
      E(9) = 2.0d0*ehartr
      E(10) = pxc
      E(14) = etf
      E(15) = ptf
      E(16) = ewgc
      E(17) = pwgc

      if (pspw_qmmm_found()) then
         e_lj     = pspw_qmmm_LJ_E()
         e_q      = pspw_qmmm_Q_E()
         e_spring = pspw_qmmm_spring_E()
         E(1)  = E(1) + e_lj + e_q + e_spring

         E(11) = e_lj
         E(12) = e_q
         E(13) = e_spring
      end if

*     **** set convergence variables ****
      deltae = (E(1)-Eold)/(dt*dble(it_in))

*     *** deltac ***
      do n=1,ispin
         do i=1,npack1
            Hpsi(i,n) = psi2(i,n) - psi1(i,n)
         end do
      end do

      do n=1,ispin
         call Pack_cc_idot(1,Hpsi(1,n),Hpsi(1,n),dbl_mb(sumi(1)+n-1))
      end do
      if (np_i.gt.1) 
     >     call D3dB_Vector_SumAll(ispin,dbl_mb(sumi(1)))

      deltac = 0.0d0
      do n=1,ispin
         if (dbl_mb(sumi(1)+n-1).gt.deltac) deltac=dbl_mb(sumi(1)+n-1)
      end do
      deltac = deltac/dte *(ne(1)+ne(2))


*     *** deltar ***
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
      value = BA_pop_stack(sumi(2))
      value = value.and.BA_pop_stack(ftest(2))
      value = value.and.BA_pop_stack(fion(2))
      value = value.and.BA_pop_stack(dnall(2))
      value = value.and.BA_pop_stack(v_wgc(2))
      value = value.and.BA_pop_stack(e_tf(2))
      value = value.and.BA_pop_stack(v_tf(2))
      value = value.and.BA_pop_stack(v_tmp(2))
      value = value.and.BA_pop_stack(xce(2))
      value = value.and.BA_pop_stack(xcp(2))
      value = value.and.BA_pop_stack(dng(2))
      value = value.and.BA_pop_stack(rho(2))
      value = value.and.BA_pop_stack(vl(2))
      value = value.and.BA_pop_stack(v_field(2))

      if ((control_version().eq.4).or.(field_exist))
     >   value = value.and.BA_pop_stack(r_grid(2))

      if (control_version().eq.4) 
     >   value = value.and.BA_pop_stack(vlr_l(2))

      value = value.and.BA_pop_stack(vc(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))

      call nwpw_timing_end(12)


      return
      end
 
