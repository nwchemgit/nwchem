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
  
      subroutine inner_loop(ispin,ne,
     >                      npack1,nfft3d,nemax,
     >                      psi1,psi2,dn,
     >                      it_in,E,deltae,deltac,deltar,
     >                      hml,lmd,
     >                      psi_r,Hpsi)
*
* $Id: inner_loop.F,v 1.4 2001-10-17 23:48:14 bylaska Exp $
*
      implicit none
      integer    ispin,ne(2)
      integer    npack1,nfft3d,nemax
      complex*16 psi1(npack1,nemax)
      complex*16 psi2(npack1,nemax)
      real*8     dn(2*nfft3d,2)
      integer    it_in
      real*8     E(*)
      real*8     deltae,deltac,deltar
      real*8     hml(2*nemax*nemax),lmd(2*nemax*nemax)

*     **** very big workspace variables ****
      real*8     psi_r(2*nfft3d,nemax)
      complex*16 Hpsi(npack1,nemax)


#include "mafdecls.fh"


*     **** local variables ****
      logical move
      integer n2ft3d,np
      integer i,j,ii,jj,n,n1(2),n2(2),it,ms,nn,ierr
      integer nx,ny,nz
      integer index,indext
      real*8  sum,Eold,eorbit,eion,ehartr,eke,elocal,enlocal
      real*8  exc,exc2,pxc,pxc2,dte,scal1,scal2,dv,dt
      real*8  e_elc_water,e_ion_water,e_water_water,deltawater


*     **** MA local variables ****
      logical value
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
      integer xcp(2),xce(2),dnall(2)
      integer natmx,fion(2),ftest(2)
      integer sumi(2)
      integer npack0

*     ***** external functions ****
      logical  control_move,psp_semicore,Waterpsp_found
      logical  pspw_SIC,pspw_SIC_relaxed
      integer  ion_nion,control_gga
      real*8   control_time_step,control_fake_mass,ion_dti
      real*8   lattice_omega,coulomb_e,ewald_e
      external control_move,psp_semicore,Waterpsp_found
      external pspw_SIC,pspw_SIC_relaxed
      external ion_nion,control_gga
      external control_time_step,control_fake_mass,ion_dti
      external lattice_omega,coulomb_e,ewald_e
      integer  control_version
      external control_version
      real*8   ion_ion_e
      external ion_ion_e

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
*    12 - allocate/deallocate time
*    13 - Hpsi and update
      real*8 times(50),tim1,tim2
      common / timing / times





      call Pack_npack(0,npack0)
      n2ft3d = 2*nfft3d
      deltawater = 0.0d0


      call current_second(tim1)
*     **** allocate MA local variables ****
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
     >        MA_push_get(mt_dcpl,(npack0),'vloc', vl(2), vl(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(n2ft3d),'rho',rho(2),rho(1))
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
      value = value.and.
     >        MA_push_get(mt_dbl,(nemax),'sumi',sumi(2),sumi(1))
      if (.not. value) call errquit('out of stack memory',0)

      call current_second(tim2)
      times(12) = times(12) + (tim2-tim1)
     
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

      if (control_version().eq.4) 
     >   call lattice_r_grid(dbl_mb(r_grid(1)))

*     ******************************************
*     ****                                  ****
*     ****   Start of steepest descent loop ****
*     ****                                  ****
*     ******************************************
      do it=1,it_in
        call dcopy(2*npack1*nemax,psi2,1,psi1,1)
        if (move) call ion_shift()
        if ((move).and.(Waterpsp_found())) call Waterpsp_shift()

*       ********************************
*       **** generate phaze factors ****
*       ********************************
        call phafac()


        call current_second(tim1)
*       *******************
*       **** get psi_r ****
*       *******************
        do n=n1(1),n2(ispin)
c           call mask_C(1,psi1(1,n))
c           call D3dB_c_Copy(1,psi1(1,n),psi_r(1,n))
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
c        call mask_C(0,dcpl_mb(dng(1)))
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
     
*       **** generate Waterpsp potential ****
        if (Waterpsp_found()) then
            call Waterpsp_Generate_V(n2ft3d,dbl_mb(r_grid(1)),
     >                                      dbl_mb(vlr_l(1)))
        end if

        call current_second(tim1)
*       ******************
*       **** get Hpsi ****
*       ******************
        if (control_version().eq.3) 
     >     call psi_H(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),dcpl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,dbl_mb(fion(1)))

        if (control_version().eq.4) 
     >     call psi_Hv4(ispin,ne,psi1,psi_r,
     >             dcpl_mb(vl(1)),dbl_mb(vlr_l(1)),
     >             dbl_mb(vc(1)),dbl_mb(xcp(1)),Hpsi,
     >             move,dbl_mb(fion(1)))


        
*       ************************************
*       **** do a steepest descent step ****
*       ************************************
        do n=1,n2(ispin)
c          call D3dB_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
c          call D3dB_cc_Sum(1,psi2(1,n),psi1(1,n),psi2(1,n))
c          call mask_C(1,psi2(1,n))
          call Pack_c_SMul(1,dte,Hpsi(1,n),psi2(1,n))
          call Pack_cc_Sum(1,psi2(1,n),psi1(1,n),psi2(1,n))
        end do
        call current_second(tim2)
        times(13) = times(13) + (tim2-tim1)

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

*          **** Waterpsp update ****
           if (Waterpsp_found()) then
              call Waterpsp_Update(0,n2ft3d,dbl_mb(r_grid(1)),
     >                                       dbl_mb(rho(1)),
     >                                       dv,dt,
     >                                       dbl_mb(fion(1)),
     >                                       deltawater)
           end if

*          **** remove ion forces using ion_FixIon ****
           call ion_FixIon(dbl_mb(fion(1)))

           call ion_optimize_step(dbl_mb(fion(1)))
        end if


*       *****************************************        
*       **** lagrange multiplier corrections ****
*       *****************************************        
        call psi_lmbda(ispin,ne,nemax,npack1,psi1,psi2,dte,
     >                 lmd,dbl_mb(tmp_L(1)),ierr)

      end do

*     *************************************
*     ***** total energy calculation ******
*     *************************************
      call current_second(tim1)
      call Parallel_np(np)

*     *** get orbital energies ****
      n = ne(1)
      nn = n*n
      do ms=1,ispin
         do ii=n1(ms),n2(ms)
           i = ii-n1(ms)
           index = (i+1) + i*n + (ms-1)*nn
c           call D3dB_cc_dot(1,psi1(1,ii),Hpsi(1,ii),sum)
           call Pack_cc_idot(1,psi1(1,ii),Hpsi(1,ii),sum)

           hml(index) =  -sum
           do jj=ii+1,n2(ms)
              j = jj-n1(ms)
              index  = (i+1) + j*n + (ms-1)*nn
              indext = (j+1) + i*n + (ms-1)*nn
c              call D3dB_cc_dot(1,psi1(1,ii),Hpsi(1,jj),sum)
              call Pack_cc_idot(1,psi1(1,ii),Hpsi(1,jj),sum)

              hml(index)  =  -sum
              hml(indext) =  -sum
           end do 
         end do
      end do
      if (np.gt.1)  call D3dB_Vector_SumAll((ispin*nn),hml)
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
        
  
*     ***** average Kohn-Sham kinetic energy ****
      call ke_ave(ispin,ne,psi1,eke)


*     **** average Kohn-Sham v_local energy ****
      call Pack_cc_dot(0,dcpl_mb(dng(1)),dcpl_mb(vl(1)),elocal)

c     call dcopy(2*npack0,dcpl_mb(vl(1)),1,dcpl_mb(tmp1(1)),1)
c     call Pack_c_unpack(0,dcpl_mb(tmp1(1)))
c     call D3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
c     elocal = 0.0d0
c     do ms=1,ispin
c        do n=n1(ms),n2(ms)
c           call D3dB_rr_Mul(1,dcpl_mb(tmp1(1)),psi_r(1,n),
c    >                       dcpl_mb(tmp2(1)))
c
c           call D3dB_rr_idot(1,psi_r(1,n),dcpl_mb(tmp2(1)),sum)
c           elocal = elocal + sum*scal1*scal2
c        end do
c     end do
c     if (np.gt.1) call D3dB_SumAll(elocal)
c     if (ispin.eq.1) elocal = 2.0d0*elocal


*     *** add in long range part ****
      if (control_version().eq.4) then
      call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(vlr_l(1)),sum)
      elocal = elocal + sum*dv
      end if

 
*     ***** average Kohn-Sham v_nonlocal energy ****
      call dcopy(2*npack1*nemax,0.0d0,0,Hpsi,1)
      call v_nonlocal(ispin,ne,psi1,Hpsi,
     >                move,dbl_mb(ftest(1)))
      enlocal = 0.0d0
      do ms=1,ispin
      do n=n1(ms),n2(ms)
         call Pack_cc_idot(1,psi1(1,n),Hpsi(1,n),sum)
         enlocal = enlocal - sum
      end do
      end do
      if (np.gt.1) call D3dB_SumAll(enlocal)
      if (ispin.eq.1) enlocal = 2.0d0*enlocal


*     *********************************
*     ******** Waterpsp energies ******
*     *********************************
      if (Waterpsp_found()) then
*        **** get free-space electron-water energy ****
         call dcopy(n2ft3d,0.0d0,0,dbl_mb(vlr_l(1)),1)
         call Waterpsp_Generate_V(n2ft3d,dbl_mb(r_grid(1)),
     >                                   dbl_mb(vlr_l(1)))
         call D3dB_rr_dot(1,dbl_mb(rho(1)),dbl_mb(vlr_l(1)),e_elc_water)
         e_elc_water = e_elc_water*dv

*        **** subtract e_elc_water from elocal ****
         elocal = elocal - e_elc_water

*        **** get the ion-water energy ****
         call Waterpsp_Energy_ion(e_ion_water)

*        **** get the water-water energy ****
         call Waterpsp_Energy_Water(e_water_water)
      end if


*     **** total energy ****
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

      if (Waterpsp_found()) then
         E(1) = E(1) + e_ion_water + e_water_water
         E(11) = e_elc_water
         E(12) = e_ion_water
         E(13) = e_water_water
      end if

*     **** SIC corrections ****
      if (pspw_SIC()) then
         call pspw_energy_SIC(ispin,psi_r,exc,pxc)
         E(1) = E(1) + exc
         E(14) = exc
         E(15) = 0.0d0
         if (pspw_SIC_relaxed()) then
            E(1)  = E(1) - pxc
            E(15) = pxc
         end if
      end if


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
      if (np.gt.1) 
     >     call D3dB_Vector_SumAll((n2(ispin)-n1(1)),
     >                             dbl_mb(sumi(1)))
      deltac = 0.0d0
      do n=n1(1),n2(ispin)
         if (dbl_mb(sumi(1)+n-1).gt.deltac) deltac=dbl_mb(sumi(1)+n-1)
      end do
      deltac = deltac/dte

c     do i=1,ion_nion()
c         fion(1,i) = (r2(1,i) - r1(1,i))/ion_dti(i)
c         fion(2,i) = (r2(2,i) - r1(2,i))/ion_dti(i)
c         fion(3,i) = (r2(3,i) - r1(3,i))/ion_dti(i)
c     end do

*     *** deltar ***
      deltar = deltawater
      if (move) then
        do i=1,ion_nion()
           sum = dsqrt( dbl_mb(fion(1)+(i-1)*3  )**2 
     >                + dbl_mb(fion(1)+(i-1)*3+1)**2 
     >                + dbl_mb(fion(1)+(i-1)*3+2)**2)
           if (sum.gt.deltar) deltar = sum
        end do
      end if
       
      call current_second(tim2)
      times(10) = times(10) + (tim2-tim1)

*     **** dealocate MA local variables ****
      call current_second(tim1)
      value = MA_pop_stack(sumi(2))
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
 
