*
* $Id$
*
  
      subroutine inner_loop_qmmm_step(verlet,
     >                      ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi0,psi1,psi2,dn,
     >                      it_in,it_sum,E,
     >                      hml,lmd,
     >                      Hpsi,psi_r,
     >                      rion,nion,rcharge,qcharge,ncharge,
     >                      fion,fcharge)
      implicit none
      logical    verlet
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax
      complex*16 psi0(npack1,nemax)
      complex*16 psi1(npack1,nemax)
      complex*16 psi2(npack1,nemax)
      real*8     dn(2*nfft3d,2)
      integer    it_in,it_sum
      real*8     E(*)
      real*8     hml(2*nemax*nemax),lmd(2*nemax*nemax)
      complex*16 Hpsi(npack1,nemax)
      real*8     psi_r(2*nfft3d,nemax)

      real*8  rion(3,*)
      integer nion
      real*8  rcharge(3,*),qcharge(*)
      integer ncharge
      real*8  fion(3,*)
      real*8  fcharge(3,*)


#include "bafdecls.fh"
#include "errquit.fh"
ccc#include "frac_occ.fh"


*     **** local variables ****
      logical move
      integer n2ft3d,np
      integer i,j,ii,jj,n,n1(2),n2(2),it,ms,nn,ierr,gga
      integer nx,ny,nz
      integer index,indext
      real*8  sum,Eold,eorbit,eion,ehartr,eke,eki,sse,ssr
      real*8  exc,exc2,pxc,pxc2,dte,dte0,scal1,scal2,dv,dt,fmass,h
      real*8  ehsic,phsic,exsic,pxsic,ehfx,phfx
c      real*8  e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav


*     **** MA local variables ****
      logical value,nose,field_exist,sic,allow_translation
      integer tmp_L(2)
      integer tmp1(2),tmp2(2)
      integer vl(2),vc(2),dng(2)
      integer rho(2),vlr_l(2),r_grid(2)
      integer xcp(2),xce(2),dnall(2)
      integer v_field(2)
      integer npack0

*     ***** external functions ****
      integer  ion_nion,control_gga
      real*8   ion_ke,ion_ion_e
      real*8   control_time_step,control_fake_mass,ion_dti
      real*8   lattice_omega,coulomb_e,ewald_e
      external ion_nion,control_gga
      external ion_ke,ion_ion_e
      external control_time_step,control_fake_mass,ion_dti
      external lattice_omega,coulomb_e,ewald_e
      logical  psp_semicore
      external psp_semicore
      integer  control_version
      external control_version
      logical  control_Nose
      external control_Nose
      real*8   Nose_e_energy,Nose_r_energy,Nose_sse,Nose_ssr
      external Nose_e_energy,Nose_r_energy,Nose_sse,Nose_ssr

*     ***** QM/MM external functions ****
c      real*8   pspw_charge_Energy_ion_qmmm
c      external pspw_charge_Energy_ion_qmmm
      logical  dplot_iteration_check
      external dplot_iteration_check

      logical  pspw_SIC,pspw_SIC_relaxed,control_allow_translation
      logical  pspw_HFX,pspw_HFX_relaxed
      external pspw_SIC,pspw_SIC_relaxed,control_allow_translation
      external pspw_HFX,pspw_HFX_relaxed

      call Parallel_np(np)
      call Pack_npack(0,npack0)

      n2ft3d = 2*nfft3d
      field_exist       = .true.
      sic               = pspw_SIC()
      gga               = control_gga()
      allow_translation = control_allow_translation()

*     **** allocate MA local variables ****
      call nwpw_timing_start(12)
      value = BA_push_get(mt_dbl,(8*nemax*nemax),
     >                     'tmp_L',tmp_L(2),tmp_L(1))
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

      if ((control_version().eq.4).or.(field_exist)) then
       value = value.and.
     >    BA_push_get(mt_dbl,(3*n2ft3d),'r_grid',r_grid(2),r_grid(1))
      end if

      value = value.and.
     >  BA_push_get(mt_dcpl,(n2ft3d),'v_field',v_field(2),v_field(1))

      value = value.and.
     >         BA_push_get(mt_dcpl,(npack0),'vl',  vl(2),  vl(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'rho',rho(2), rho(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'dng',dng(2), dng(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xcp',xcp(2), xcp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'xce',xce(2), xce(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(4*nfft3d),'dnall',dnall(2),dnall(1))
      if (.not.value) 
     > call errquit('inner_loop_qmmm_step:pushing stack',0, MA_ERR)
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
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      dt    = control_time_step()
      fmass = control_fake_mass()
      dte   = dt*dt/fmass
      if (.not. verlet) dte=0.5d0*dte

      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()

      if ((control_version().eq.4).or.(field_exist))
     >   call lattice_r_grid(dbl_mb(r_grid(1)))


*     ******************************************
*     ****                                  ****
*     **** Start of molecular dynamics loop ****
*     ****                                  ****
*     ******************************************
      do it=1,it_in
        call dcopy(2*npack1*nemax,psi1,1,psi0,1)
        call dcopy(2*npack1*nemax,psi2,1,psi1,1)
        call dcopy(3*ncharge,0.0d0,0,fcharge,1)

*       *** skip ion_shift if newton step ***
        if (nose)   call Nose_Shift()


*       ********************************
*       **** generate phaze factors ****
*       ********************************
        call ion_rion_reset(rion)
        call phafac()
        if (control_version().eq.3) call ewald_phafac()

        call nwpw_timing_start(11)
*       *******************
*       **** get psi_r ****
*       *******************
        do n=n1(1),n2(ispin)
           call Pack_c_Copy(1,psi1(1,n),psi_r(1,n))
           call Pack_c_unpack(1,psi_r(1,n))
           call D3dB_cr_fft3b(1,psi_r(1,n))
           call D3dB_r_Zero_Ends(1,psi_r(1,n))
        end do


*       *********************
*       **** generate dn ****
*       *********************
        call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
c        if (fractional) then
c          do ms=1,ispin
c             do n=n1(ms),n2(ms)
c                do i=1,n2ft3d
c                   dn(i,ms) = dn(i,ms)
c     >                      + scal2*(psi_r(i,n)**2)
c     >                       *dbl_mb(fweight(1)+n-1)
c                end do
c             end do
c             call D3dB_r_Zero_Ends(1,dn(1,ms))
c          end do
c        else
          do ms=1,ispin
             do n=n1(ms),n2(ms)
                do i=1,n2ft3d
                   dn(i,ms) = dn(i,ms) + scal2*(psi_r(i,n)**2)
                end do
             end do
             call D3dB_r_Zero_Ends(1,dn(1,ms))
          end do
c        end if


*       **********************
*       **** generate dng ****
*       **********************
        call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho(1)))
        call D3dB_r_SMul(1,scal1,dbl_mb(rho(1)),dcpl_mb(tmp1(1)))
        call D3dB_rc_fft3f(1,dcpl_mb(tmp1(1)))
        call Pack_c_pack(0,dcpl_mb(tmp1(1)))
        call Pack_c_Copy(0,dcpl_mb(tmp1(1)),dcpl_mb(dng(1)))

*       ********************************************************
*       **** generate dnall - used for semicore corrections ****
*       ********************************************************
        if (psp_semicore(0)) then
           call semicore_density_update()
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
     >               fion)


*       *** long-range psp for charge systems ***
        if (control_version().eq.4) then
          call v_lr_local(dbl_mb(r_grid(1)),
     >                    dbl_mb(vlr_l(1)))
          if (move) then
             call grad_v_lr_local(dbl_mb(r_grid(1)),
     >                            dbl_mb(rho(1)),
     >                            fion)
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
        else if ((gga.eq.1).or.(gga.eq.2)) then
          call v_bwexc(gga,n2ft3d,ispin,dbl_mb(dnall(1)),
     >                      1.0d0,1.0d0,
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)))
        else if (gga.eq.3) then
          call v_bwexc(1,n2ft3d,ispin,dbl_mb(dnall(1)),
     >                      0.75d0,1.0d0,
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)))
        end if

     
*       **********************************************
*       **** generate other real-space potentials ****
*       **********************************************
        if (field_exist) then

           call dcopy(n2ft3d,0.0d0,0,dbl_mb(v_field(1)),1)


*          **** generate charge potential ****
c          call pspw_charge_Generate_V_qmmm(rcharge,qcharge,ncharge,
c     >                               n2ft3d,
c     >                               dbl_mb(r_grid(1)),
c     >                               dbl_mb(v_field(1)))
c
        end if


*       ******************
*       **** get Hpsi ****
*       ******************
        if (control_version().eq.3)
     >  call psi_H(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dcpl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,fion)

        if (control_version().eq.4)
     >     call psi_Hv4(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),dbl_mb(vlr_l(1)),
     >             dbl_mb(v_field(1)),field_exist,
     >             dbl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,fion)


*       **********************
*       **** get ewald force *
*       **********************
*       **** get the ewald force ****
        if (control_version().eq.3) call ewald_f(fion)

*       **** get the free-space ion force ****
        if (control_version().eq.4) call ion_ion_f(fion)
        
*       ************************
*       **** get semicoreforce *
*       ************************
        if (psp_semicore(0)) then
           call semicore_xc_F(ispin,dbl_mb(xcp(1)),fion)
        end if

*       **************************
*       **** get qmmm  forces ****
*       **************************

c        call pspw_charge_rho_Fcharge_qmmm(rcharge,qcharge,ncharge,
c     >                                    n2ft3d,
c     >                                    dbl_mb(r_grid(1)),
c     >                                    dbl_mb(rho(1)),
c     >                                    dv,
c     >                                    fcharge)
c        call pspw_charge_Fion_Fcharge_qmmm(rcharge,qcharge,ncharge,
c     >                                         fion,fcharge)


*       *****************************************
*       **** remove ion forces using ion_FixIon *
*       *****************************************
        call ion_FixIon(fion)


*       **************************
*       **** do a verlet step ****
*       **************************
        if (verlet) then
*          **** constant temperature ****
           if (nose) then 
             sse = Nose_sse()
             ssr = Nose_ssr()
             do n=1,n2(ispin)
              call Pack_c_SMul(1,0.5d0*dte,Hpsi(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,-1.0d0,psi0(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,1.0d0,psi1(1,n),psi2(1,n))
              call Pack_c_SMul(1,2.0d0*sse,psi2(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,1.0d0,psi0(1,n),psi2(1,n))
             end do

*          **** constant energy ****
           else
             do n=1,n2(ispin)
              call Pack_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,-1.0d0,psi0(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,2.0d0,psi1(1,n),psi2(1,n))
             end do

           end if

*       **************************
*       **** do a newton step ****
*       **************************
        else
           do n=1,n2(ispin)
              call Pack_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,dt,psi0(1,n),psi2(1,n))
              call Pack_cc_Sum(1,psi2(1,n),psi1(1,n),psi2(1,n))
           end do

        end if


*       *****************************************        
*       **** lagrange multiplier corrections ****
*       *****************************************        

        !**** orthoganality constraint ****
        dte0 = dte
        if (nose) dte0 = dte*sse
c        if (fractional) then
c        call psi_lmbda2(ispin,ne,nemax,npack1,psi1,psi2,
c     >                  dte0,dbl_mb(fweight(1)),
c     >                  lmd,dbl_mb(tmp_L(1)),ierr)
        if (sic) then
        call psi_lmbda_sic(ispin,ne,nemax,npack1,psi1,psi2,dte0,
     >                 lmd,dbl_mb(tmp_L(1)),ierr)
        else
        call psi_lmbda(ispin,ne,nemax,npack1,psi1,psi2,dte0,
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
                call Pack_cc_idot(1,psi2(1,i),psi0(1,i),sum)
                eke = eke+sum
             end do
             if (np.gt.1) call D3dB_SumAll(eke)
             eke = (ne(1)+ne(2) - eke)
             if (ispin.eq.1) eke = 2.0d0*eke
             eke = 0.5d0*(fmass/(dt*dt))*eke
             eki = ion_ke()
             call Nose_Verlet_Step(eke,eki)
          else
              eke = 0.0d0
              do i=1,n2(ispin)
                call Pack_cc_idot(1,psi0(1,i),psi0(1,i),sum)
                eke = eke+sum
              end do
              if (np.gt.1) call D3dB_SumAll(eke)
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
         call dplot_iteration((it+it_sum),ispin,ne,psi1,dn,psi_r)
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

*     *** get orbital energies ****
      n = ne(1)
      nn = n*n
      do ms=1,ispin
         do ii=n1(ms),n2(ms)
           i = ii-n1(ms)
           index = (i+1) + i*n + (ms-1)*nn
           call Pack_cc_idot(1,psi1(1,ii),Hpsi(1,ii),sum)
           hml(index) = -sum
           do jj=ii+1,n2(ms)
              j = jj-n1(ms)
              index  = (i+1) + j*n + (ms-1)*nn
              indext = (j+1) + i*n + (ms-1)*nn
              call Pack_cc_idot(1,psi1(1,ii),Hpsi(1,jj),sum)
              hml(index)  = -sum
              hml(indext) = -sum
           end do 
         end do
      end do
      if (np.gt.1) call D3dB_Vector_SumAll(ispin*nn,hml)

      eorbit = 0.0d0
      do ms=1,ispin
         do ii=1,ne(ms)
            index = (ii) + (ii-1)*n + (ms-1)*nn
            sum = hml(index)
c            if (fractional) sum=sum*dbl_mb(fweight(1)+ii-1+(ms-1)*ne(1))
            eorbit = eorbit + sum
         end do
      end do
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



        
*     **** velocity an kinetic energy of psi ****
      h = 1.0d0/(2.0d0*dt)
      eke = 0.0d0
      do i=1,n2(ispin)
         call Pack_c_SMul(1,-h,psi0(1,i),psi0(1,i))
         call Pack_cc_daxpy(1,h,psi2(1,i),psi0(1,i))
         call Pack_cc_idot(1,psi0(1,i),psi0(1,i),sum)
         eke = eke+sum
      end do
      if (np.gt.1) call D3dB_SumAll(eke)
      eke = eke*fmass
      if (ispin.eq.1) eke = 2.0d0*eke

 
c      e_ionmm = pspw_charge_Energy_ion_qmmm(rcharge,qcharge,ncharge)


*     **** total energy ****
      Eold=E(1)
      E(2) = eorbit + eion + exc - ehartr - pxc 
      E(3) = eke
      E(4) = 0.0d0
      E(5) = eorbit
      E(6) = ehartr
      E(7) = exc
      E(8) = eion


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



*     **** Energy and Energy**2 sum ***
      E(25) = E(25) + E(2)
      E(26) = E(26) + E(2)*E(2)

*     **** output Forces for Fei ***
      call fei_output(E(2),fion)



      if (nose) then
        E(9) = Nose_e_energy()
        E(10) = Nose_r_energy()
        E(1) = E(2)+E(3)+E(4)+E(9)+E(10)
      else
        E(1) = E(2)+E(3)+E(4)
      end if




      call nwpw_timing_end(10)

*     **** dealocate MA local variables ****
 333  continue
      call nwpw_timing_start(12)

      value =           BA_pop_stack(dnall(2))
      value = value.and.BA_pop_stack(xce(2))
      value = value.and.BA_pop_stack(xcp(2))
      value = value.and.BA_pop_stack(dng(2))
      value = value.and.BA_pop_stack(rho(2))
      value = value.and.BA_pop_stack(vl(2))
      value = value.and.BA_pop_stack(v_field(2))

      if ((control_version().eq.4).or.(field_exist)) 
     >  value = value.and.BA_pop_stack(r_grid(2))

      if (control_version().eq.4) 
     >  value = value.and.BA_pop_stack(vlr_l(2))

      value = value.and.BA_pop_stack(vc(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and.BA_pop_stack(tmp_L(2))
      if (.not.value) 
     > call errquit('inner_loop_md:popping stack',1, MA_ERR)

      call nwpw_timing_end(12)

      return
      end
 
