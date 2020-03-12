*
* $Id$
*
*  ************************************************************
*  *                                                          *
*  *                   BAND cpmd routine                      *
*  *                                                          *
*  ************************************************************
  
      subroutine band_inner_loop_md(verlet,sa_alpha,
     >                      ispin,ispinq,ne,neq,nbrill,nbrillq,
     >                      nfft3d,
     >                      psi0_tag,psi1_tag,psi2_tag,dn,
     >                      it_in,it_sum,
     >                      E,
     >                      hml_tag,
     >                      psir_tag,Hpsi_tag,
     >                      calc_pressure,pressure,p1,p2)

      implicit none
      logical verlet
      real*8  sa_alpha(2)
      integer ispin,ispinq,ne(2),neq(2),nbrill,nbrillq,nfft3d
      integer psi0_tag,psi1_tag,psi2_tag
      real*8  dn(nfft3d,*)
      integer it_in,it_sum
      real*8  E(*)
      integer hml_tag

*     **** very big workspace variables ****
      integer psir_tag,Hpsi_tag

      logical calc_pressure
      real*8  pressure,p1,p2,stress(3,3)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical move,sic,hfx,done,nose
      integer n2ft3d,np_i,np_j,np_k,np,nbq,nemaxq
      integer i,n,n1(2),n2(2),it,ms,ierr,shifto
      integer shiftA,shiftB,nshiftA,nshiftB
      integer psi0_shift,psi1_shift,psi2_shift
      integer nx,ny,nz,nbi,nbj,ki,kj,occ_tag,nb,lmd_shift
      integer gga,psi_shift,hpsi_shift,hml_shift,nshift1,nshift2
      real*8  sum,Eold,eorbit,eion,ehartr,eke,eki,elocal,enlocal
      real*8  exc,exc2,pxc,pxc2,dte,scal1,scal2,dv,dt
      real*8  ehsic,phsic,exsic,pxsic,ehfx,phfx
      real*8  deltamm,weight,sa1,sa2,fmass,ssr,sse,r,s,h
      real*8  e_lj, e_q, e_spring
      !real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav



*     **** MA local variables ****
      logical value,field_exist,fei
      integer tmp_L(2),trho(2),txcp(2),txce(2)
      integer tmp1(2),tmp2(2)
      integer vl(2),vc(2),dng(2)
      integer rho(2)
      integer xcp(2),xce(2),dnall(2)
      integer natmx,fion(2),ftest(2)
      integer sumi(2)
      integer npack0
      integer npack1

*     ***** external functions ****
      logical  control_move,cpsp_semicore,control_Nose
      logical  pspw_SIC,pspw_SIC_relaxed,control_Fei
      logical  pspw_HFX,pspw_HFX_relaxed
      integer  ion_nion,control_gga,cpsi_data_nsize
      real*8   control_time_step,control_fake_mass,ion_dti
      real*8   lattice_omega,c_coulomb_e,ewald_e
      external control_move,cpsp_semicore,control_Nose
      external pspw_SIC,pspw_SIC_relaxed,control_Fei
      external pspw_HFX,pspw_HFX_relaxed
      external ion_nion,control_gga,cpsi_data_nsize
      external control_time_step,control_fake_mass,ion_dti
      external lattice_omega,c_coulomb_e,ewald_e
      integer  control_version,cpsi_data_get_next,cpsi_data_get_chnk
      external control_version,cpsi_data_get_next,cpsi_data_get_chnk
      real*8   brillioun_weight,Nose_sse,Nose_ssr,Nose_dXe,Nose_dXr
      external brillioun_weight,Nose_sse,Nose_ssr,Nose_dXe,Nose_dXr
      real*8   Nose_e_energy,Nose_r_energy,ion_ke
      external Nose_e_energy,Nose_r_energy,ion_ke

      double precision Dneall_m_trace
      external         Dneall_m_trace
      logical  C3dB_cr_pfft3_queue_filled
      external C3dB_cr_pfft3_queue_filled
      logical  Pneb_w_push_get_block,Pneb_w_pop_stack
      external Pneb_w_push_get_block,Pneb_w_pop_stack

*     ***** QM/MM external functions ****
      logical  pspw_charge_found,pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_charge_found,pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      logical  ion_disp_on,control_only_lda
      external ion_disp_on,control_only_lda
      real*8   ion_disp_energy
      external ion_disp_energy

      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)
      deltamm = 0.0d0
      gga     = control_gga()
      nshiftA = 2*npack1
      nshiftB = 2*nfft3d
      fei = control_Fei()

      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1)+1
      n2(2) = neq(1)+neq(2)
      nemaxq = neq(1)+neq(2)

      call nwpw_timing_start(12)

*     **** allocate MA local variables ****
      value = Pneb_w_push_get_block(0,1,8,tmp_L)
      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vc',  vc(2),  vc(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'vloc', vl(2), vl(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'rho',rho(2),rho(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(npack0),'dng',dng(2), dng(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nfft3d),'xcp',xcp(2), xcp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nfft3d),'xce',xce(2), xce(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nfft3d),'dnall',dnall(2),dnall(1))
      natmx = ion_nion()
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'fion',fion(2),fion(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'ftest',ftest(2),ftest(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemaxq),'sumi',sumi(2),sumi(1))
      if (.not.control_only_lda()) then
         call D3dB_n2ft3d(1,n2ft3d)
         value = BA_push_get(mt_dbl,(2*n2ft3d),'trho',trho(2),trho(1))
         value = value.and.
     >           BA_push_get(mt_dbl,(2*n2ft3d),'txcp',txcp(2),txcp(1))
         value = value.and.
     >           BA_push_get(mt_dbl,(2*n2ft3d),'txce',txce(2),txce(1))
      end if
      if (.not. value) 
     > call errquit('band_inner_loop:out of stack memory',0, MA_ERR)
      call dcopy(2*nfft3d,0.0d0,0,dbl_mb(xcp(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,dbl_mb(xce(1)),1)

      call nwpw_timing_end(12)
     
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      move = .true.

      nose = control_Nose()
      sse = 1.0d0
      ssr = 1.0d0

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

*     ******************************************
*     ****                                  ****
*     **** Start of molecular dynamics loop ****
*     ****                                  ****
*     ******************************************

      do it=1,it_in
        call cpsi_data_copyall(psi1_tag,psi0_tag)
        call cpsi_data_copyall(psi2_tag,psi1_tag)
        if (verlet) call ion_shift()
        if (nose.and.verlet) call Nose_shift()


*       ********************************
*       **** generate phaze factors ****
*       ********************************
        call cphafac()
        call cphafac_k()
        call ewald_phafac()
        call nwpw_timing_start(11)

*       *******************
*       **** get psi_r ****
*       *******************
        nbi = 1
        ki  = 1
        shiftA =  cpsi_data_get_chnk(psi1_tag,nbi)

        nbj = 1
        kj  = 1
        shiftB =  cpsi_data_get_chnk(psir_tag,nbj)

        done = .false.
        do while (.not.done)

           if (nbi.le.nbrillq) then
              call C3dB_cr_pfft3b_queuein(nbi,dbl_mb(shiftA))
              shiftA = shiftA + nshiftA
              ki = ki+1
              if (ki.gt.nemaxq) then
                 nbi = nbi + 1
                 ki  = 1
                 if (nbi.le.nbrillq) then
                    shiftA =  cpsi_data_get_chnk(psi1_tag,nbi)
                 endif
              end if
           end if

           if ((C3dB_cr_pfft3_queue_filled()).or.(nbi.gt.nbrillq)) then
              call C3dB_cr_pfft3b_queueout(nbj,dbl_mb(shiftB))
              shiftB = shiftB + nshiftB
              kj = kj+1
              if (kj.gt.nemaxq) then
                 nbj = nbj + 1
                 kj  = 1
                 if (nbj.le.nbrillq) then
                    shiftB =  cpsi_data_get_chnk(psir_tag,nbj)
                 endif
              end if
           endif
           done = ((nbi.gt.nbrillq).and.(nbj.gt.nbrillq))
        end do

*       **** set the occupations ****
        call cpsi_data_set_next(psir_tag,cpsi_data_get_next(psi1_tag))


*       *********************
*       **** generate dn ****
*       *********************
        call dcopy(ispin*nfft3d,0.0d0,0,dn,1)


        occ_tag = cpsi_data_get_next(psir_tag)
        if (occ_tag.gt.0) then
         do nbq=1,nbrillq
          weight = brillioun_weight(nbq)
          shiftB = cpsi_data_get_chnk(psir_tag,nbq)
          shifto = cpsi_data_get_chnk(occ_tag,nbq)
          do ms=1,ispin
          do n=n1(ms),n2(ms)
             call C3dB_cr_aSqrpy(1,dbl_mb(shifto)*weight,
     >                             dbl_mb(shiftB),
     >                             dn(1,ms))
             shiftB = shiftB + nshiftB
             shifto = shifto + 1
          end do
          end do
         end do
        else
         do nbq=1,nbrillq
          weight = brillioun_weight(nbq)
          shiftB = cpsi_data_get_chnk(psir_tag,nbq)
          do ms=1,ispin
          do n=n1(ms),n2(ms)
             call C3dB_cr_aSqrpy(1,weight,dbl_mb(shiftB),dn(1,ms))
             shiftB = shiftB + nshiftB
          end do
          end do
         end do
        end if

*       *** reduce over brillioun zone ***
        call K1dB_Vector_SumAll(ispin*nfft3d,dn)
        call dscal(ispin*nfft3d,scal2,dn,1)


*       **********************
*       **** generate dng ****
*       **********************
        call C3dB_rrc_Sum(1,dn,dn(1,ispin),dcpl_mb(rho(1)))
        call C3dB_rc_pfft3f(1,0,dcpl_mb(rho(1)))
        call C3dB_c_SMul1(1,scal1,dcpl_mb(rho(1)))
        call Cram_c_pack(0,dcpl_mb(rho(1)))
        call Cram_c_Copy(0,dcpl_mb(rho(1)),dcpl_mb(dng(1)))


*       ********************************************************
*       **** generate dnall - used for semicore corrections ****
*       ********************************************************
        if (cpsp_semicore(0)) then
           if (move .or. (it.eq.1)) call c_semicore_density_update()
           call c_semicore_density(dcpl_mb(rho(1)))
           call C3dB_r_SMul1(1,0.5d0,dcpl_mb(rho(1)))
           do ms=1,ispin
             call C3dB_rr_Sum(1,dn(1,ms),
     >                        dcpl_mb(rho(1)),
     >                        dbl_mb(dnall(1)+(ms-1)*nfft3d))
           end do
        else
           call dcopy(ispin*nfft3d,dn,1,dbl_mb(dnall(1)),1)
        end if
        call nwpw_timing_end(11)


*       *****************************************
*       **** generate local pseudopotential  ****
*       **** and also get force if move true ****
*       *****************************************
        call cpsp_v_local(dcpl_mb(vl(1)),
     >               move,
     >               dcpl_mb(dng(1)),
     >               dbl_mb(fion(1)))


*       ************************************
*       **** generate coulomb potential ****
*       ************************************
        call c_coulomb_v(dcpl_mb(dng(1)),dcpl_mb(vc(1)))


*       *************************************************
*       **** generate exchange-correlation potential ****
*       *************************************************
        if (control_only_lda()) then
           call vxc(nfft3d,ispin,
     >              dbl_mb(dnall(1)),
     >              dbl_mb(xcp(1)),
     >              dbl_mb(xce(1)),
     >              dcpl_mb(rho(1)))

c        !*** this code needs to be optimized!!! ****
         !*** GGA calculation ****
         else
            do ms=1,ispin
              call C3dB_D3dB_r_Copy(1,dbl_mb(dnall(1)+(ms-1)*nfft3d),
     >                                dbl_mb(trho(1) +(ms-1)*n2ft3d))
            end do

            call v_bwexc_all(gga,n2ft3d,ispin,
     >                       dbl_mb(trho(1)),
     >                       dbl_mb(txcp(1)),
     >                       dbl_mb(txce(1)))

            do ms=1,ispin
              call D3dB_C3dB_r_Copy(1,dbl_mb(txcp(1)+(ms-1)*n2ft3d),
     >                                dbl_mb(xcp(1) +(ms-1)*nfft3d))
            end do
            call D3dB_C3dB_r_Copy(1,dbl_mb(txce(1)),
     >                              dbl_mb(xce(1)))
         end if



*       ******************
*       **** get Hpsi ****
*       ******************
        call nwpw_timing_start(13)
        call cpsi_H(ispin,neq,
     >              psi1_tag,
     >              psir_tag,
     >              dcpl_mb(vl(1)),
     >              dcpl_mb(vc(1)),
     >              dbl_mb(xcp(1)),
     >              Hpsi_tag,
     >              move,dbl_mb(fion(1)))


*       *******************************************
*       **** get ewald and semicore forces     ****
*       *******************************************
*       **** get the ewald force ****
        call ewald_f(dbl_mb(fion(1)))

*       **** get the semicore force ****
        if (cpsp_semicore(0)) then
           call c_semicore_xc_F(ispin,dbl_mb(xcp(1)),dbl_mb(fion(1)))
        end if

c*          **** get the qmmm force ****
c           if (pspw_qmmm_found()) call pspw_qmmm_fion(dbl_mb(fion(1)))

*       **** get the dispersion force ****
        if (ion_disp_on()) call ion_disp_force(dbl_mb(fion(1)))

*       **** remove ion forces using ion_FixIon ****
        if (fei) 
     >    call dcopy(3*natmx,dbl_mb(fion(1)),1,dbl_mb(ftest(1)),1)
        call ion_FixIon(dbl_mb(fion(1)))

        
*       **********************************
*       **** do an verlet/newton step ****
*       **********************************

*       **************************
*       **** do a verlet step ****
*       **************************
        call cpsi_data_update(psi2_tag)
        if (verlet) then
*          **** constant temperature ****
           if (nose) then
             sse = Nose_sse()
             ssr = Nose_ssr()
             do nb=1,nbrillq
                psi0_shift = cpsi_data_get_chnk(psi0_tag,nb)
                psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
                psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
                Hpsi_shift = cpsi_data_get_chnk(Hpsi_tag,nb)
                call band_cpmd_subupdate1(nb,npack1,(neq(1)+neq(2)),
     >                             dbl_mb(psi0_shift),
     >                             dbl_mb(psi1_shift),
     >                             dbl_mb(psi2_shift),
     >                             dbl_mb(Hpsi_shift),
     >                             dte,sse)
             end do
             call ion_nose_step(ssr,dbl_mb(fion(1)))

*          **** constant energy ****
           else
             do nb=1,nbrillq
                psi0_shift = cpsi_data_get_chnk(psi0_tag,nb)
                psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
                psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
                Hpsi_shift = cpsi_data_get_chnk(Hpsi_tag,nb)
                call band_cpmd_subupdate2(nb,npack1,(neq(1)+neq(2)),
     >                             dbl_mb(psi0_shift),
     >                             dbl_mb(psi1_shift),
     >                             dbl_mb(psi2_shift),
     >                             dbl_mb(Hpsi_shift),
     >                             dte,sa1,sa2)
             end do
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
           do nb=1,nbrillq
              psi0_shift = cpsi_data_get_chnk(psi0_tag,nb)
              psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
              psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
              Hpsi_shift = cpsi_data_get_chnk(Hpsi_tag,nb)
              call band_cpmd_subupdate3(nb,npack1,(neq(1)+neq(2)),
     >                             dbl_mb(psi0_shift),
     >                             dbl_mb(psi1_shift),
     >                             dbl_mb(psi2_shift),
     >                             dbl_mb(Hpsi_shift),
     >                             dte,s*dt*sa_alpha(1))
           end do
           call ion_newton_step(dbl_mb(fion(1)),sa_alpha(2)*r)
        end if
        call cpsi_data_noupdate(psi2_tag)


        call nwpw_timing_end(13)

*       *****************************************        
*       **** lagrange multiplier corrections ****
*       *****************************************        
        call cpsi_data_update(psi2_tag)
        do nb=1,nbrillq
          psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
          psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
c          lmd_shift  = cpsi_data_get_chnk(lmd_tag,nb)
          do ms=1,ispinq
          if (neq(ms).gt.0) then
             call cpsi_lmbda(ms,nb,npack1,
     >                   dbl_mb(psi1_shift),
     >                   dbl_mb(psi2_shift),
     >                   dte,
     >                   dcpl_mb(tmp_L(1)),ierr)

          end if
          end do
        end do
        call cpsi_data_noupdate(psi2_tag)


*       **************************
*       *** update thermostats ***
*       **************************
        if (nose) then
          if (verlet) then
             eke = 0.0d0
             do nb=1,nbrillq
                psi0_shift = cpsi_data_get_chnk(psi0_tag,nb)
                psi2_shift = cpsi_data_get_chnk(psi2_tag,nb)
                do i=1,n2(ispin)
                   call Cram_cc_idot(nb,
     >                               dbl_mb(psi2_shift),
     >                               dbl_mb(psi0_shift),
     >                               sum)
                   eke = eke+sum*brillioun_weight(nb)
                   psi0_shift = psi0_shift + (2*npack1)
                   psi2_shift = psi2_shift + (2*npack1)
                end do
             end do
             if (np.gt.1) call Parallel_SumAll(eke)
             eke = (ne(1)+ne(2) - eke)
             if (ispin.eq.1) eke = 2.0d0*eke
             eke = 0.5d0*(fmass/(dt*dt))*eke
             eki = ion_ke()
             call Nose_Verlet_Step(eke,eki)
          else
              eke = 0.0d0
              do nb=1,nbrillq
                psi0_shift = cpsi_data_get_chnk(psi0_tag,nb)
                do i=1,n2(ispin)
                   call Cram_cc_idot(nb,
     >                               dbl_mb(psi0_shift),
     >                               dbl_mb(psi0_shift),
     >                               sum)
                   eke = eke+sum*brillioun_weight(nb)
                   psi0_shift = psi0_shift + (2*npack1)
                end do
              end do
              if (np.gt.1) call Parallel_SumAll(eke)
              if (ispin.eq.1) eke = 2.0d0*eke
              eke = eke*fmass
              eki = ion_ke()
              call Nose_Newton_Step(eke,eki)
          end if
        end if


      end do

*     **** if newton then skip energy calculations ****
      if (.not. verlet) goto 333


*     *************************************
*     ***** total energy calculation ******
*     *************************************
      call nwpw_timing_start(10)
      call Parallel_np(np)
      call Parallel3d_np_i(np_i)
      call Parallel3d_np_j(np_j)
      call Parallel3d_np_k(np_k)

*     *** get orbital energies ****
      eorbit = 0.0d0
      call cpsi_data_update(hml_tag)
      do nb=1,nbrillq
         psi_shift  = cpsi_data_get_chnk(psi1_tag,nb)
         hpsi_shift = cpsi_data_get_chnk(Hpsi_tag,nb)
         hml_shift  = cpsi_data_get_chnk(hml_tag,nb)
         call Pneb_ffw_Multiply(0,nb,
     >                  dbl_mb(psi_shift),
     >                  dbl_mb(hpsi_shift),npack1,
     >                  dbl_mb(hml_shift))
         call Pneb_w_scal(0,nb,-1.0d0,dbl_mb(hml_shift))
         call Pneb_w_trace(0,nb,dbl_mb(hml_shift),sum)
         eorbit = eorbit + sum*brillioun_weight(nb)
      end do
      call cpsi_data_noupdate(hml_tag)
      call K1dB_SumAll(eorbit)
      if (ispin.eq.1) eorbit = eorbit + eorbit

*     **** get ewald energy ****
      eion = ewald_e()

*     **** get coulomb energy ****
      ehartr = c_coulomb_e(dcpl_mb(dng(1)))

*     **** get exchange-correlation energy ****
      call C3dB_rr_dot(1,dbl_mb(dnall(1)),dbl_mb(xce(1)),exc)
      call C3dB_rr_dot(1,dn,dbl_mb(xcp(1)),pxc)
      if (ispin.eq.1) then
         exc= exc + exc 
         pxc= pxc + pxc 
      else
         call C3dB_rr_dot(1,dbl_mb(dnall(1)+nfft3d),
     >                      dbl_mb(xce(1)),exc2)
         call C3dB_rr_dot(1,dn(1,2),dbl_mb(xcp(1)+nfft3d),pxc2)
         exc= exc + exc2
         pxc= pxc + pxc2
      end if
      exc = exc*dv
      pxc = pxc*dv


*     **** velocity and kinetic energy of psi ****
      h   = 1.0d0/(2.0d0*dt)
      eke = 0.0d0
      call cpsi_data_update(psi0_tag)
      do nb=1,nbrillq
         psi0_shift  = cpsi_data_get_chnk(psi0_tag,nb)
         psi2_shift  = cpsi_data_get_chnk(psi2_tag,nb)
         do i=1,nemaxq
            call Cram_c_SMul1(nb,-h,dbl_mb(psi0_shift))
            call Cram_cc_daxpy(nb,h,
     >                         dbl_mb(psi2_shift),
     >                         dbl_mb(psi0_shift))
            call Cram_cc_idot(nb,dbl_mb(psi0_shift),
     >                           dbl_mb(psi0_shift),sum)
            eke = eke+sum*brillioun_weight(nb)
            psi0_shift = psi0_shift + (2*npack1)
            psi2_shift = psi2_shift + (2*npack1)
         end do
      end do
      call cpsi_data_noupdate(psi0_tag)
      if (np.gt.1) call Parallel_SumAll(eke)
      eke = eke*fmass
      if (ispin.eq.1) eke = eke + eke


*     **** total energy ****
      Eold=E(1)
      E(2) = eorbit + eion + exc - ehartr - pxc
      E(3) = eke
      E(4) = ion_ke()
      E(5) = eorbit
      E(6) = ehartr
      E(7) = exc
      E(8) = eion
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

*     **** Dispersion energy ****
      if (ion_disp_on()) then
         E(33) = ion_disp_energy()
         E(2)  = E(2) + E(33)
      end if


*     **** Energy and Energy**2 sum ***
      E(25) = E(25) + E(2)
      E(26) = E(26) + E(2)*E(2)

*     **** output Forces for Fei ***
      if (fei) call fei_output(E(2),dbl_mb(ftest(1)))

      if (nose) then
        E(9)  = Nose_e_energy()
        E(10) = Nose_r_energy()
        E(1)  = E(2)+E(3)+E(4)+E(9)+E(10)
      else
        E(1) = E(2)+E(3)+E(4)
      end if


c*     ******** pressure ******
c      if (calc_pressure) then
c
c*        ***** average Kohn-Sham v_nonlocal energy ****
c        call dcopy(2*npack1*nemaxq,0.0d0,0,Hpsi,1)
c        call v_nonlocal(ispin,neq,psi1,Hpsi,
c     >                .false.,dbl_mb(ftest(1)),fractional,occ1)
c        enlocal = 0.0d0
c        do ms=1,ispin
c        do n=n1(ms),n2(ms)
c         call Pack_cc_idot(1,psi1(1,n),Hpsi(1,n),sum)
c         if (fractional) sum=sum*occ1(n)
c         enlocal = enlocal - sum
c        end do
c        end do
c        if (np.gt.1) call Parallel_SumAll(enlocal)
c        if (ispin.eq.1) enlocal = 2.0d0*enlocal
c
c
c        call cgsd_pressure_stress(ispin,neq,psi1,
c     >                            dbl_mb(dnall(1)),
c     >                            dcpl_mb(dng(1)),
c     >                            enlocal,exc,pxc,
c     >                            pressure,p1,p2,stress)
c      end if


*      **** write ecce data ****
       call ecce_print_module_entry('task gradient')

       call ion_ecce()
       call ecce_print1('total energy', mt_dbl, E(2), 1)
       call ecce_print2('total gradient', mt_dbl, dbl_mb(fion(1)),
     $        3,3,natmx)
       call ecce_print1('gradient norm', mt_dbl, E(1), 1)
       call ecce_print1('orbital gradient norm', mt_dbl, E(4), 1)
c       call ecce_print1('gradient max', mt_dbl, E(1), 1)
       call ecce_print_module_exit('task gradient', 'ok')
       
      call nwpw_timing_end(10)

*     **** dealocate MA local variables ****
 333  continue
      call nwpw_timing_start(12)
      value = .true.
      if (.not.control_only_lda()) then
         value = value.and.BA_pop_stack(txce(2))
         value = value.and.BA_pop_stack(txcp(2))
         value = value.and.BA_pop_stack(trho(2))
      end if
      value = value.and.BA_pop_stack(sumi(2))
      value = value.and.BA_pop_stack(ftest(2))
      value = value.and.BA_pop_stack(fion(2))
      value = value.and.BA_pop_stack(dnall(2))
      value = value.and.BA_pop_stack(xce(2))
      value = value.and.BA_pop_stack(xcp(2))
      value = value.and.BA_pop_stack(dng(2))
      value = value.and.BA_pop_stack(rho(2))
      value = value.and.BA_pop_stack(vl(2))
      value = value.and.BA_pop_stack(vc(2))
      value = value.and.Pneb_w_pop_stack(tmp_L)
      if (.not. value) 
     > call errquit('band_inner_loop:popping stack',0,MA_ERR)

      call nwpw_timing_end(12)

      return
      end
 
      subroutine band_cpmd_subupdate1(nb,nfft3d,nn,
     >                             psi0,psi1,psi2,Hpsi,dte,sse)
      implicit none
      integer    nb,nfft3d,nn
      complex*16 psi0(nfft3d,nn)
      complex*16 psi1(nfft3d,nn)
      complex*16 psi2(nfft3d,nn)
      complex*16 Hpsi(nfft3d,nn)
      real*8     dte,sse

      integer n
*     ************************************
*     **** do a nose step ****
*     ************************************
      do n=1,nn
        call Cram_c_SMul(nb,0.5d0*dte,Hpsi(1,n),psi2(1,n))
        call Cram_cc_daxpy(nb,-1.0d0,psi0(1,n),psi2(1,n))
        call Cram_cc_Sum2(nb,psi1(1,n),psi2(1,n))
        call Cram_c_SMul1(nb,2.0d0*sse,psi2(1,n))
        call Cram_cc_Sum2(nb,psi0(1,n),psi2(1,n))
      end do

      return
      end


      subroutine band_cpmd_subupdate2(nb,nfft3d,nn,
     >                                psi0,psi1,psi2,Hpsi,
     >                                dte,sa1,sa2)
      implicit none
      integer    nb,nfft3d,nn
      complex*16 psi0(nfft3d,nn)
      complex*16 psi1(nfft3d,nn)
      complex*16 psi2(nfft3d,nn)
      complex*16 Hpsi(nfft3d,nn)
      real*8     dte,sa1,sa2

      integer n
*     **********************
*     **** verlet  step ****
*     **********************
      do n=1,nn
        call Cram_c_SMul(nb,sa1*dte,Hpsi(1,n),psi2(1,n))
        call Cram_cc_daxpy(nb,-1.0d0*sa2,psi0(1,n),psi2(1,n))
        call Cram_cc_daxpy(nb,2.0d0*sa1,psi1(1,n),psi2(1,n))
      end do

      return
      end


      subroutine band_cpmd_subupdate3(nb,nfft3d,nn,
     >                                psi0,psi1,psi2,Hpsi,
     >                                dte,dt)
      implicit none
      integer    nb,nfft3d,nn
      complex*16 psi0(nfft3d,nn)
      complex*16 psi1(nfft3d,nn)
      complex*16 psi2(nfft3d,nn)
      complex*16 Hpsi(nfft3d,nn)
      real*8     dte,dt

      integer n
*     **********************
*     **** newton  step ****
*     **********************
      do n=1,nn
        call Cram_c_SMul(nb,dte,Hpsi(1,n),psi2(1,n))
        call Cram_cc_daxpy(nb,dt,psi0(1,n),psi2(1,n))
        call Cram_cc_Sum2(nb,psi1(1,n),psi2(1,n))
      end do

      return
      end



      
