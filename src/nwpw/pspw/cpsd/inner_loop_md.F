*
* $Id: inner_loop_md.F,v 1.4 2001-11-27 00:05:03 bylaska Exp $
*
*  ************************************************************
*  *                MPI cpmd routine                          *
*  *                                                          *
*  *  This is a developing cpmdv3 parallel code wrtten in     *
*  *  Fortran and MPI.                                        *
*  *                                                          *
*  *    + mpl message passing library used                    *
*  *                                                          *
*  *    + ngp is used instead of nfft in this proceudure      *
*  *                                                          *
*  *    + error checking is based on aimd.h parameters        *
*  *      then control file                                   *
*  ************************************************************
  
      subroutine inner_loop_md(verlet,
     >                      ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi0,psi1,psi2,dn,
     >                      it_in,E,
     >                      hml,lmd,
     >                      Hpsi,psi_r)
      implicit none
      logical    verlet
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax
      complex*16 psi0(npack1,nemax)
      complex*16 psi1(npack1,nemax)
      complex*16 psi2(npack1,nemax)
      real*8     dn(2*nfft3d,2)
      integer    it_in
      real*8     E(*)
      real*8     hml(2*nemax*nemax),lmd(2*nemax*nemax)

*     **** very big workspace variables ****
      complex*16 Hpsi(npack1,nemax)
      real*8     psi_r(2*nfft3d,nemax)


#include "mafdecls.fh"


*     **** local variables ****
      logical move
      integer n2ft3d,np
      integer i,j,ii,jj,n,n1(2),n2(2),it,ms,nn,ierr
      integer nx,ny,nz
      integer index,indext
      real*8  sum,Eold,eorbit,eion,ehartr,eke,eki,sse,ssr
      real*8  exc,exc2,pxc,pxc2,dte,dte0,scal1,scal2,dv,dt,fmass,h
      real*8  e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav


*     **** MA local variables ****
      logical value,nose
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
      integer natmx,fion(2),ftest(2)
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
      logical  pspw_qmmm_found
      real*8   pspw_qmmm_ke
      real*8   pspw_qmmm_Energy_ionMM
      real*8   pspw_qmmm_Energy_QMMM
      real*8   pspw_qmmm_Energy_MMMM
      real*8   pspw_qmmm_Energy_pol
      real*8   pspw_qmmm_Energy_vib
      real*8   pspw_qmmm_Energy_cav
      external pspw_qmmm_found
      external pspw_qmmm_ke
      external pspw_qmmm_Energy_ionMM
      external pspw_qmmm_Energy_QMMM
      external pspw_qmmm_Energy_MMMM
      external pspw_qmmm_Energy_pol
      external pspw_qmmm_Energy_vib
      external pspw_qmmm_Energy_cav



*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
*    10 - energy_time
*    11 - density_time
*    12 - allocate/deallocate_time
      real*8 times(50),tim1,tim2
      common / timing / times



      call Parallel_np(np)
      call Pack_npack(0,npack0)

      n2ft3d = 2*nfft3d

*     **** allocate MA local variables ****
      call current_second(tim1)
      value = MA_push_get(mt_dbl,(8*nemax*nemax),
     >                     'tmp_L',tmp_L(2),tmp_L(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))

      if (control_version().eq.3) then
      value = value.and.
     >        MA_push_get(mt_dcpl,(npack0),'vc',  vc(2),  vc(1))
      end if

      if (control_version().eq.4) then
       value = value.and.
     >        MA_push_get(mt_dbl,(n2ft3d),'vc',vc(2),vc(1))

       value = value.and.
     >        MA_push_get(mt_dbl,(n2ft3d),'vlr_l',vlr_l(2),vlr_l(1))
       value = value.and.
     >    MA_push_get(mt_dbl,(3*n2ft3d),'r_grid',r_grid(2),r_grid(1))
      end if

      value = value.and.
     >         MA_push_get(mt_dcpl,(npack0),'vl',  vl(2),  vl(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(n2ft3d),'rho',rho(2), rho(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(npack0),'dng',dng(2), dng(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(4*nfft3d),'xcp',xcp(2), xcp(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(4*nfft3d),'xce',xce(2), xce(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(4*nfft3d),'dnall',dnall(2),dnall(1))
      natmx = ion_nion()
      value = value.and.
     >        MA_push_get(mt_dbl,(3*natmx),'fion',fion(2),fion(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(3*natmx),'ftest',ftest(2),ftest(1))
      if (.not.value) call errquit('out of stack memory',0)

      call current_second(tim2)
      times(12) = times(12) + (tim2-tim1)
     
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

      if (control_version().eq.4)
     >   call lattice_r_grid(dbl_mb(r_grid(1)))


*     ******************************************
*     ****                                  ****
*     **** Start of molecular dynamics loop ****
*     ****                                  ****
*     ******************************************
      do it=1,it_in
        call dcopy(2*npack1*nemax,psi1,1,psi0,1)
        call dcopy(2*npack1*nemax,psi2,1,psi1,1)
*       *** skip ion_shift if newton step ***
        if (verlet) call ion_shift()
        if (nose)   call Nose_Shift()
        if ((verlet).and.(pspw_qmmm_found())) call pspw_qmmm_shift()

*       ********************************
*       **** generate phaze factors ****
*       ********************************
        call phafac()


        call current_second(tim1)
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
        do ms=1,ispin
           do n=n1(ms),n2(ms)
              do i=1,n2ft3d
                 dn(i,ms) = dn(i,ms) + scal2*(psi_r(i,n)**2)
              end do
           end do
           call D3dB_r_Zero_Ends(1,dn(1,ms))
        end do
        

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
        call current_second(tim2)
        times(11) = times(11) + (tim2-tim1)





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
        if (control_gga().eq.0) then
          call vxc(n2ft3d,ispin,dbl_mb(dnall(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)),
     >                      dcpl_mb(tmp1(1)))
        else
          call v_bwexc(n2ft3d,ispin,dbl_mb(dnall(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)))
        end if
     
*       **** generate QM/MM potential ****
        if (pspw_qmmm_found()) then
            call pspw_qmmm_Generate_V(n2ft3d,dbl_mb(r_grid(1)),
     >                                      dbl_mb(vlr_l(1)))
        end if


*       ******************
*       **** get Hpsi ****
*       ******************
        if (control_version().eq.3)
     >  call psi_H(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),dcpl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,dbl_mb(fion(1)))

        if (control_version().eq.4)
     >     call psi_Hv4(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),dbl_mb(vlr_l(1)),
     >             dbl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,dbl_mb(fion(1)))



*       **********************
*       **** get ewald force *
*       **********************
*       **** get the ewald force ****
        if (control_version().eq.3) call ewald_f(dbl_mb(fion(1)))

*       **** get the free-space ion force ****
        if (control_version().eq.4) call ion_ion_f(dbl_mb(fion(1)))
        
*       ************************
*       **** get semicoreforce *
*       ************************
        if (psp_semicore(0)) then
           call semicore_xc_F(ispin,dbl_mb(xcp(1)),dbl_mb(fion(1)))
        end if

*       *****************************************
*       **** remove ion forces using ion_FixIon *
*       *****************************************
        call ion_FixIon(dbl_mb(fion(1)))

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
             call ion_nose_step(ssr,dbl_mb(fion(1)))

*          **** constant energy ****
           else
             do n=1,n2(ispin)
              call Pack_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,-1.0d0,psi0(1,n),psi2(1,n))
              call Pack_cc_daxpy(1,2.0d0,psi1(1,n),psi2(1,n))
             end do

*            **** QM/MM Verlet update ****
             if (pspw_qmmm_found()) then
                call pspw_qmmm_Update(2,n2ft3d,dbl_mb(r_grid(1)),
     >                                       dbl_mb(rho(1)),
     >                                       dv,dt,
     >                                       dbl_mb(fion(1)),
     >                                       sum)
             end if
             call ion_verlet_step(dbl_mb(fion(1)))
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

*          **** QM/MM Newton update ****
           if (pspw_qmmm_found()) then
              call pspw_qmmm_Update(1,n2ft3d,dbl_mb(r_grid(1)),
     >                                       dbl_mb(rho(1)),
     >                                       dv,dt,
     >                                       dbl_mb(fion(1)),
     >                                       sum)
           end if
           call ion_newton_step(dbl_mb(fion(1)))

        end if


*       *****************************************        
*       **** lagrange multiplier corrections ****
*       *****************************************        
        dte0 = dte
        if (nose) dte0 = dte*sse
        call psi_lmbda(ispin,ne,nemax,npack1,psi1,psi2,dte0,
     >                 lmd,dbl_mb(tmp_L(1)),ierr)

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
             eke = 0.5d0*(fmass/(dt*dt))*(ne(1)+ne(2) - eke)
             eki = ion_ke()
             call Nose_Verlet_Step(eke,eki)
          else
              eke = 0.0d0
              do i=1,n2(ispin)
                call Pack_cc_idot(1,psi0(1,i),psi0(1,i),sum)
                eke = eke+sum
              end do
              if (np.gt.1) call D3dB_SumAll(eke)
              eke = eke*fmass
              eki = ion_ke()
              call Nose_Newton_Step(eke,eki)
          end if
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
      call current_second(tim1)

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
            eorbit = eorbit + hml(index)
         end do
      end do
      if (ispin.eq.1) eorbit = eorbit+eorbit
        

*     **** get ewald energy ****
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

 

*     **** total energy ****
      Eold=E(1)
      E(2) = eorbit + eion + exc - ehartr - pxc
      E(3) = eke
      E(4) = ion_ke()
      E(5) = eorbit
      E(6) = ehartr
      E(7) = exc
      E(8) = eion

*     ******** QM/MM energies ******
      if (pspw_qmmm_found()) then
         e_ionmm = pspw_qmmm_Energy_ionMM()
         e_qmmm  = pspw_qmmm_Energy_QMMM()
         e_mmmm  = pspw_qmmm_Energy_MMMM()
         e_pol   = pspw_qmmm_Energy_pol()
         e_vib   = pspw_qmmm_Energy_vib()
         e_cav   = pspw_qmmm_Energy_cav()

         E(2) = E(2) + e_ionmm
     >               + e_qmmm
     >               + e_mmmm
     >               + e_pol
     >               + e_vib
     >               + e_cav
         E(4) = E(4) + pspw_qmmm_ke()
         E(11) = e_qmmm
         E(12) = e_mmmm
         E(13) = e_pol
         E(14) = e_vib
         E(15) = e_cav
      end if



      if (nose) then
        E(9) = Nose_e_energy()
        E(10) = Nose_r_energy()
        E(1) = E(2)+E(3)+E(4)+E(9)+E(10)
      else
        E(1) = E(2)+E(3)+E(4)
      end if

      call current_second(tim2)
      times(10) = times(10) + (tim2-tim1)

*     **** dealocate MA local variables ****
 333  continue
      call current_second(tim1)

      value = MA_pop_stack(ftest(2))
      value = MA_pop_stack(fion(2))
      value = MA_pop_stack(dnall(2))
      value = MA_pop_stack(xce(2))
      value = MA_pop_stack(xcp(2))
      value = MA_pop_stack(dng(2))
      value = MA_pop_stack(rho(2))
      value = MA_pop_stack(vl(2))
      if (control_version().eq.4) value = MA_pop_stack(r_grid(2))
      if (control_version().eq.4) value = MA_pop_stack(vlr_l(2))
      value = MA_pop_stack(vc(2))
      value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp1(2))
      value = MA_pop_stack(tmp_L(2))

      call current_second(tim2)
      times(12) = times(12) + (tim2-tim1)

      return
      end
 
