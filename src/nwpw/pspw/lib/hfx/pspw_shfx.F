
*     ****************************
*     *                    	 *
*     *   pspw_potential_sHFX0   *
*     *                          *
*     ****************************
      subroutine pspw_potential_sHFX0(ispin0,psi_r,Hpsi_r)
      implicit none
      integer ispin0
      real  psi_r(*)
      real Hpsi_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

      integer istart,iend,jstart,jend,imodn,imodtask
      integer ms,l,q,n,indx1,indx2,Levels,neq(2)
      integer requests(5),reqcnt

      integer  Butter_Levels,Dneall_na_ptr
      external Butter_Levels,Dneall_na_ptr

*     ***** now do exchange as normal ****
!$OMP MASTER
      ehfx = 0.0d0
      phfx = 0.0d0
!$OMP END MASTER
      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then

         if (replicated) then

*           **** butterfly algorithm ****
            if (butterfly) then
               call Dneall_neq(neq)
               Levels = Butter_Levels(npj)
               do ms=1,ispin
                  call Parallel_shared_vector_szero(.false.,nrsize,
     >                          real_mb(Hpsi_r_replicated(1)))
                  call Parallel_shared_vector_scopy(.true.,
     >                       neq(ms)*n2ft3d,
     >                       psi_r(1+(ms-1)*neq(1)*n2ft3d),
     >                       real_mb(psi_r_replicated(1)))

                  do l=0,Levels-1
                     call D1dBs_Brdcst_step(l,
     >                       int_mb(Dneall_na_ptr(ms)),-1,
     >                       n2ft3d,
     >                       real_mb(psi_r_replicated(1)),
     >                       requests,reqcnt)

                     call Butter_indexes(l,taskid_j,npj,
     >                       int_mb(Dneall_na_ptr(ms)),
     >                       istart,iend,jstart,jend,
     >                       imodn,imodtask)
                     call pspw_potential_sHFX_sub2(solver_type,
     >                                  istart,iend,
     >                                  jstart,jend,
     >                                  imodn,imodtask,
     >                                  n2ft3d,
     >                                  real_mb(psi_r_replicated(1)),
     >                                  real_mb(Hpsi_r_replicated(1)),
     >                                  ehfx)

                     call D1dB_WaitAll(requests,reqcnt)
                  end do

                  call Butter_indexes_L1(taskid_j,npj,
     >                       int_mb(Dneall_na_ptr(ms)),
     >                       istart,iend,jstart,jend,
     >                       imodn,imodtask)
                  if (jend.ge.jstart)
     >               call pspw_potential_sHFX_sub2(solver_type,
     >                               istart,iend,
     >                               jstart,jend,
     >                               imodn,imodtask,
     >                               n2ft3d,
     >                               real_mb(psi_r_replicated(1)),
     >                               real_mb(Hpsi_r_replicated(1)),
     >                               ehfx)
                  call Butter_indexes_L2(taskid_j,npj,
     >                       int_mb(Dneall_na_ptr(ms)),
     >                       istart,iend,jstart,jend,
     >                       imodn,imodtask)
                  call pspw_potential_sHFX_sub2(solver_type,
     >                               istart,iend,
     >                               jstart,jend,
     >                               imodn,imodtask,
     >                               n2ft3d,
     >                               real_mb(psi_r_replicated(1)),
     >                               real_mb(Hpsi_r_replicated(1)),
     >                               ehfx)

                  do l=Levels-1,0,-1
                     call D1dBs_Reduce_step(l,
     >                       int_mb(Dneall_na_ptr(ms)),-1,
     >                       n2ft3d,
     >                       real_mb(Hpsi_r_replicated(1)),
     >                       real_mb(psi_r_replicated(1)))
                  end do
                  call SAXPY_OMP(neq(ms)*n2ft3d,hfx_parameter,
     >                       real_mb(Hpsi_r_replicated(1)),1,
     >                       Hpsi_r(1+(ms-1)*neq(1)*n2ft3d),1)
               end do


*              *** apply hfx_parameter ****
!$OMP MASTER
               ehfx = ehfx*hfx_parameter

               if (ispin.eq.1) ehfx = ehfx + ehfx
!$OMP END MASTER
               call Parallel_SumAll(ehfx)
!$OMP MASTER
               phfx = 2.0d0*ehfx
!$OMP END MASTER

*           **** reduceall algorithm ****
            else
            call Parallel_shared_vector_szero(.false.,nrsize,
     >                                 real_mb(psi_r_replicated(1)))
            call Parallel_shared_vector_szero(.true.,nrsize,
     >                                 real_mb(Hpsi_r_replicated(1)),1)
            do q=1,neqall
               call Dneall_qton(q,n)
               indx1 = (q-1)*n2ft3d + 1
               indx2 = psi_r_replicated(1)+(n-1)*n2ft3d
               call Parallel_shared_vector_scopy(.true.,n2ft3d,
     >                                      psi_r(indx1),real_mb(indx2))
            end do
            call D1dBs_Vector_SumAll(nrsize,
     >                       real_mb(psi_r_replicated(1)))
            call pspw_potential_sHFX_sub(ispin0,
     >                                  real_mb(psi_r_replicated(1)),
     >                                  real_mb(Hpsi_r_replicated(1)))
            call D1dBs_Vector_SumAll(nrsize,
     >                       real_mb(Hpsi_r_replicated(1)))
            do q=1,neqall
               call Dneall_qton(q,n)
               indx1 = Hpsi_r_replicated(1)+(n-1)*n2ft3d
               indx2 = (q-1)*n2ft3d + 1
               call SAXPY_OMP(n2ft3d,1.0d0,real_mb(indx1),1,
     >                                     Hpsi_r(indx2),1)
            end do
            end if

         else
            call pspw_potential_sHFX_sub(ispin0,psi_r,Hpsi_r)
         end if

      end if

      return
      end



*     *************************
*     *                       *
*     *    pspw_energy_sHFX0  *
*     *                       *
*     *************************
      subroutine pspw_energy_sHFX0(ispin0,psi_r,ehfx_out,phfx_out)
      implicit none
      integer ispin0
      real  psi_r(*)
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

      integer q,n,indx1,indx2

      if (((norbs(1)+norbs(2)).ne.0).and.(.not.relaxed)) then

c     **** calculate HFX energy  ****
         if (replicated) then

            call Parallel_shared_vector_szero(.true.,nrsize,
     >                                    real_mb(psi_r_replicated(1)))
            do q=1,neqall
               call Dneall_qton(q,n)
               indx1 = (q-1)*n2ft3d + 1
               indx2 = psi_r_replicated(1)+(n-1)*n2ft3d
               call Parallel_shared_vector_scopy(.true.,n2ft3d,
     >                                     psi_r(indx1),real_mb(indx2))
            end do
            call D1dBs_Vector_SumAll(nrsize,
     >                       real_mb(psi_r_replicated(1)))
            call pspw_energy_sHFX_sub(ispin0,
     >                               real_mb(psi_r_replicated(1)),
     >                               ehfx_out,phfx_out)

         else

            call pspw_energy_sHFX_sub(ispin0,psi_r,ehfx_out,phfx_out)

         end if

c     **** nothing to do ****
      else
         ehfx_out = ehfx
         phfx_out = phfx
      end if

      return
      end




*     ********************************
*     *                    	     *
*     *    pspw_potential_sHFX_orb   *
*     *                              *
*     ********************************
      subroutine pspw_potential_sHFX_orb(ms,psi_r,
     >                                  orb_r,Horb_r)
      implicit none
      integer    ms
      real     psi_r(*)
      real     orb_r(*)
      real     Horb_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n1,n2,n3
      integer dn(2),vij(2),tmp1(2),index2
      real  scal1,scal2,dv,seh
      real*8  eh,ph

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      real     coulomb_screened_s_e
      external coulomb_screened_s_e


      call nwpw_timing_start(33)
      if ((norbs(ms).ne.0).and.relaxed) then
        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_real,(n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'vij_hfx',vij(2),vij(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call Parallel_shared_vector_szero(.false.,n2ft3d,
     >                                    real_mb(dn(1)))
        call Parallel_shared_vector_szero(.false.,n2ft3d,
     >                                   real_mb(vij(1)))
        call Parallel_shared_vector_szero(.true.,n2ft3d,
     >                                    real_mb(tmp1(1)))

        scal1 = real(1.0d0/dble(n1*n2*n3))
        scal2 = real(1.0d0/lattice_omega())
        dv = scal1/scal2

        do j=1,norbs(ms)
           index2 = (int_mb(orbital_list(1,ms)+j-1)-1)*n2ft3d + 1

*          **** generate dnij for Vij  ****
           call D3dBs_rr_Mul(1,psi_r(index2),orb_r,real_mb(dn(1)))
           call D3dBs_r_SMul1(1,scal2,real_mb(dn(1)))
           call D3dBs_r_Zero_Ends(1,real_mb(dn(1)))

*          ***** screened coulomb solver ****
           if (solver_type.eq.1) then
             call D3dBs_r_SMul1(1,scal1,real_mb(dn(1)))
             call D3dBs_rc_pfft3f(1,0,real_mb(dn(1)))
             call Packs_c_pack(0,real_mb(dn(1)))

*            **** get Ecoul energy ****
             eh = dble(coulomb_screened_s_e(real_mb(dn(1))))

*            **** generate Vcoul ****
             call coulomb_screened_s_v(real_mb(dn(1)),real_mb(vij(1)))
             call Packs_c_unpack(0,real_mb(vij(1)))
             call D3dBs_cr_pfft3b(1,0,real_mb(vij(1)))

*          ***** free-space coulomb solver ****
           else
              call coulomb2_s_v(real_mb(dn(1)),real_mb(vij(1)))
              call D3dBs_rr_dot(1,real_mb(dn(1)),real_mb(vij(1)),seh)
              eh = dble(0.5*seh*dv)
           end if

*          **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
           eh = eh*hfx_parameter
           call D3dBs_r_SMul1(1,real(hfx_parameter),real_mb(vij(1)))
           if (ispin.eq.1) eh = eh + eh
           ph = 2.0d0*eh

*          **** generate (Vij)*psi_r ***
           call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                        psi_r(index2),
     >                        real_mb(tmp1(1)))
           call D3dBs_r_Zero_Ends(1,real_mb(tmp1(1)))

*          **** add -(Vij)*psi_r to Hpsi_r ***
           call D3dBs_rr_Sub(1,Horb_r,
     >                        real_mb(tmp1(1)),
     >                        Horb_r)
        end do

        value = value.and.BA_pop_stack(tmp1(2))
        value = value.and.BA_pop_stack(vij(2))
        value = value.and.BA_pop_stack(dn(2))
        if (.not. value) 
     >    call errquit('pspw_potential_sHFX_orb:popping stack memory',0,
     &       MA_ERR)
      end if
      call nwpw_timing_end(33)
      return
      end


c***************** sub/replicated routines *****************************

*     ********************************
*     *                    	     *
*     *    pspw_potential_sHFX_sub   *
*     *                              *
*     ********************************
      subroutine pspw_potential_sHFX_sub(ispin0,psi_r,Hpsi_r)
      implicit none
      integer ispin0
      real    psi_r(*)
      real    Hpsi_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value,done
      integer i,j,n1,n2,n3,ms
      integer dn(2),vij(2),tmp1(2),tmp2(2),index1,index2
      integer i1,j1,k1,NN
      integer i2,j2,k2
      integer i3,j3,k3
      real    scal1,scal2,dv,seh
      real*8  eh,ph,ss,teh
      integer center(3)
      real*8  rcenter(3)
      integer taskid,icount
      real*8 cpu0,cpu1
      integer ktaskjid,kcompute(2)

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      logical  D3dBs_rc_pfft3_queue_filled,D3dBs_cr_pfft3_queue_filled
      external D3dBs_rc_pfft3_queue_filled,D3dBs_cr_pfft3_queue_filled
      logical  pspw_hfx_localize_closeenough
      external pspw_hfx_localize_closeenough

      real*8   pspw_hfx_localize_switchr
      external pspw_hfx_localize_switchr

      real     icoulomb_screened_s_e
      external icoulomb_screened_s_e
      real     icoulomb_screened_small_s_e
      external icoulomb_screened_small_s_e

      call Parallel2d_taskid_i(taskid)
      icount = 0
      call current_second(cpu0)

!$OMP MASTER
      ehfx = 0.0d0
      phfx = 0.0d0
!$OMP END MASTER
      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then
        !call D3dB_nx(1,n1)
        !call D3dB_ny(1,n2)
        !call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_real,(n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'vij_hfx',vij(2),vij(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'tmp2_hfx',tmp2(2),tmp2(1))
        NN = norbs(1)*(norbs(1)+1)/2
        value = value.and.
     >          BA_push_get(mt_int,NN,'kcmpute',kcompute(2),kcompute(1))
        if (.not.value) 
     >     call errquit('pspw_potential_HFX_sub:out of stack memory',
     >                  0,MA_ERR)
        call Parallel_shared_vector_szero(.false.,
     >                                    n2ft3d,real_mb(dn(1)))
        call Parallel_shared_vector_szero(.false.,
     >                                    n2ft3d,real_mb(vij(1)))
        call Parallel_shared_vector_szero(.false.,
     >                                    n2ft3d,real_mb(tmp1(1)))
        call Parallel_shared_vector_szero(.true.,
     >                                    n2ft3d,real_mb(tmp2(1)))

        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        scal1 = real(1.0d0/dble(n1*n2*n3))
        scal2 = real(1.0d0/lattice_omega())
        dv = scal1/scal2

        if (localize_on.and.has_smallgrid) then
           call D3dB_nx(3,n1)
           call D3dB_ny(3,n2)
           call D3dB_nz(3,n3)
           scal1 = real(1.0d0/dble(n1*n2*n3))
        end if

*       **** compute kcompute ****
        ktaskjid = 0

*       ***** screened coulomb solver ****
        if (solver_type.eq.1) then
        do ms=1,ispin0
           if (norbs(ms).eq.0) go to 898
           call Parallel_shared_vector_zero(.true.,norbs(ms),
     >                                dbl_mb(ehfx_orb(1,ms)))
           NN = norbs(ms)*(norbs(ms)+1)/2

*          **** compute kcompute ****
           i1 = 1
           j1 = 1
           do k1=1,NN
              if (pspw_hfx_localize_closeenough(
     >                   int_mb(orbital_list(1,ms)+i1-1),
     >                   int_mb(orbital_list(1,ms)+j1-1))) then
                 int_mb(kcompute(1)+k1-1) = ktaskjid
                 ktaskjid = mod(ktaskjid+1,npj)
              else
                 int_mb(kcompute(1)+k1-1) = npj+1
              end if
              j1 = j1 + 1
              if (j1.gt.i1) then
                 j1 = 1
                 i1 = i1 + 1
              end if
           end do

           i1 = 1
           j1 = 1
           k1 = 1
           i2 = 1
           j2 = 1
           k2 = 1
           i3 = 1
           j3 = 1
           k3 = 1
           done = .false.
           do while (.not.done)

              if ((k1.le.NN).and.
     >            (.not.D3dBs_rc_pfft3_queue_filled())) then

                 !if (mod(k1,npj).eq.taskid_j) then
                 if (int_mb(kcompute(1)+k1-1).eq.taskid_j) then
                 if (pspw_hfx_localize_closeenough(
     >                   int_mb(orbital_list(1,ms)+i1-1),
     >                   int_mb(orbital_list(1,ms)+j1-1))) then

                    index1 =(int_mb(orbital_list(1,ms)+i1-1)-1)*n2ft3d+1
                    index2 =(int_mb(orbital_list(1,ms)+j1-1)-1)*n2ft3d+1

*                   **** generate dnij for Vij  ****
                    call D3dBs_rr_Mul(1,psi_r(index2),
     >                                  psi_r(index1),real_mb(dn(1)))

                    call D3dBs_r_SMul1(1,scal2*scal1,real_mb(dn(1)))
                    call D3dBs_r_Zero_Ends(1,real_mb(dn(1)))

                    if (localize_on.and.has_smallgrid) then
                       call pspw_hfx_localize_center_ovlp(
     >                          int_mb(orbital_list(1,ms)+i1-1),
     >                          int_mb(orbital_list(1,ms)+j1-1),center)
                       call D3dBs_rc_pfft3f_queuein_center(2,center,
     >                                             real_mb(dn(1)))
                    else
                       call D3dBs_rc_pfft3f_queuein(0,real_mb(dn(1)))
                    end if
c                   
                   

                 end if
                 end if

                 k1 = k1 + 1
                 j1 = j1 + 1
                 if (j1.gt.i1) then
                    j1 = 1
                    i1 = i1 + 1
                 end if
              end if
 
              if ((     ((D3dBs_rc_pfft3_queue_filled()).or.(k1.gt.NN))
     >            .and.(k2.le.NN)).and.
     >            (.not.D3dBs_cr_pfft3_queue_filled())) then

                 !if (mod(k2,npj).eq.taskid_j) then
                 if (int_mb(kcompute(1)+k2-1).eq.taskid_j) then
                 if (pspw_hfx_localize_closeenough(
     >                   int_mb(orbital_list(1,ms)+i2-1),
     >                   int_mb(orbital_list(1,ms)+j2-1))) then

                    ss = pspw_hfx_localize_switchr(
     >                   int_mb(orbital_list(1,ms)+i2-1),
     >                   int_mb(orbital_list(1,ms)+j2-1))

                    if (localize_on.and.has_smallgrid) then
                       call D3dBs_rc_pfft3f_queueout_center(2,
     >                                              real_mb(dn(1)))
                       eh = dble(icoulomb_screened_small_s_e(
     >                                                real_mb(dn(1))))
                       call coulomb_screened_small_s_v(real_mb(dn(1)),
     >                                                 real_mb(vij(1)))
                    else
                       call D3dBs_rc_pfft3f_queueout(0,real_mb(dn(1)))
                       eh = dble(icoulomb_screened_s_e(real_mb(dn(1))))
                       call coulomb_screened_s_v(real_mb(dn(1)),
     >                                           real_mb(vij(1)))
                    end if

*                   **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
                    eh = eh*hfx_parameter*ss
                    if (ispin0.eq.1) eh = eh + eh
                    ph = 2.0d0*eh
!$OMP MASTER
                    ehfx = ehfx - eh
                    phfx = phfx - ph
                    dbl_mb(ehfx_orb(1,ms)+i2-1) 
     >               = dbl_mb(ehfx_orb(1,ms)+i2-1) - eh
!$OMP END MASTER
                    if (i2.ne.j2) then
!$OMP MASTER
                       ehfx = ehfx - eh
                       phfx = phfx - ph
                       dbl_mb(ehfx_orb(1,ms)+i2-1) 
     >                  = dbl_mb(ehfx_orb(1,ms)+i2-1) - eh
!$OMP END MASTER
                    end if

                    if (localize_on.and.has_smallgrid) then
                       call pspw_hfx_localize_center_ovlp(
     >                          int_mb(orbital_list(1,ms)+i2-1),
     >                          int_mb(orbital_list(1,ms)+j2-1),center)
                       call D3dBs_cr_pfft3b_queuein_center(2,center,
     >                                             real_mb(vij(1)))
                    else
                       call D3dBs_cr_pfft3b_queuein(0,real_mb(vij(1)))
                    end if
                    
                 end if
                 end if

                 k2 = k2 + 1
                 j2 = j2 + 1
                 if (j2.gt.i2) then
                    j2 = 1
                    i2 = i2 + 1
                 end if
              end if

              if ((D3dBs_cr_pfft3_queue_filled()).or.(k2.gt.NN)) then

                 !if (mod(k3,npj).eq.taskid_j) then
                 if (int_mb(kcompute(1)+k3-1).eq.taskid_j) then
                 if (pspw_hfx_localize_closeenough(
     >                   int_mb(orbital_list(1,ms)+i3-1),
     >                   int_mb(orbital_list(1,ms)+j3-1))) then

                    index1 =(int_mb(orbital_list(1,ms)+i3-1)-1)*n2ft3d+1
                    index2 =(int_mb(orbital_list(1,ms)+j3-1)-1)*n2ft3d+1

                    ss = pspw_hfx_localize_switchr(
     >                   int_mb(orbital_list(1,ms)+i3-1),
     >                   int_mb(orbital_list(1,ms)+j3-1))

                    if (localize_on.and.has_smallgrid) then
                       call D3dBs_cr_pfft3b_queueout_center(2,
     >                                              real_mb(vij(1)))
                    else
                       call D3dBs_cr_pfft3b_queueout(0,real_mb(vij(1)))
                    end if


*                   **** apply hfx_parameter ****
                    call D3dBs_r_SMul1(1,real(hfx_parameter*ss),
     >                                real_mb(vij(1)))

*                   **** generate (Vij)*psi_r ***
                    call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                                 psi_r(index2),
     >                                 real_mb(tmp1(1)))
                    call D3dBs_r_Zero_Ends(1,real_mb(tmp1(1)))

*                   **** add -(Vij)*psi_r to Hpsi_r ***
                    call D3dBs_rr_Sub2(1,real_mb(tmp1(1)),
     >                                    Hpsi_r(index1))

                    !**** include off diagonal terms ****
                    if (i3.ne.j3) then
*                      **** generate (Vij)*psi_r ***
                       call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                                 psi_r(index1),
     >                                 real_mb(tmp2(1)))
                       call D3dBs_r_Zero_Ends(1,real_mb(tmp2(1)))

*                      **** add -(Vij)*psi_r to Hpsi_r ***
                       call D3dBs_rr_Sub2(1,real_mb(tmp2(1)),
     >                                       Hpsi_r(index2))
                    end if
                 end if
                 end if

                 k3 = k3 + 1
                 j3 = j3 + 1
                 if (j3.gt.i3) then
                    j3 = 1
                    i3 = i3 + 1
                 end if
              end if
              done = (k1.gt.NN).and.(k2.gt.NN).and.(k3.gt.NN)
           end do !**** while ****
           call Parallel_Vector_SumAll(norbs(ms),
     >                                 dbl_mb(ehfx_orb(1,ms)))

  898      continue
       end do !**** ms *****

*       ***** free-space coulomb solver ****
        else
        k1 = 1
        do ms=1,ispin0
        do i=1,norbs(ms)
!$OMP MASTER
         dbl_mb(ehfx_orb(1,ms)+i-1) = 0.0d0
!$OMP END MASTER
         do j=1,i
           if (mod(k1,npj).eq.taskid_j) then
           if (pspw_hfx_localize_closeenough(
     >                   int_mb(orbital_list(1,ms)+i-1),
     >                   int_mb(orbital_list(1,ms)+j-1))) then

              index1 = (int_mb(orbital_list(1,ms)+i-1)-1)*n2ft3d + 1
              index2 = (int_mb(orbital_list(1,ms)+j-1)-1)*n2ft3d + 1

              ss = pspw_hfx_localize_switchr(
     >             int_mb(orbital_list(1,ms)+i-1),
     >             int_mb(orbital_list(1,ms)+j-1))

*             **** generate dnij for Vij  ****
              call D3dBs_rr_Mul(1,psi_r(index2),psi_r(index1),
     >                           real_mb(dn(1)))
              call D3dBs_r_SMul1(1,scal2,real_mb(dn(1)))
              call D3dBs_r_Zero_Ends(1,real_mb(dn(1)))
   
              call coulomb2_s_v(real_mb(dn(1)),real_mb(vij(1)))
              call D3dBs_rr_idot(1,real_mb(dn(1)),real_mb(vij(1)),seh)
              eh = dble(0.5*seh*dv)

*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
              eh = eh*hfx_parameter*ss
              call D3dBs_r_SMul1(1,real(hfx_parameter*ss),
     >                              real_mb(vij(1)))
              if (ispin0.eq.1) eh = eh + eh
              ph = 2.0d0*eh


!$OMP MASTER
              ehfx = ehfx - eh
              phfx = phfx - ph
              dbl_mb(ehfx_orb(1,ms)+i-1) =dbl_mb(ehfx_orb(1,ms)+i-1)-eh
!$OMP END MASTER

*             **** generate (Vij)*psi_r ***
              call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                             psi_r(index2),
     >                             real_mb(tmp1(1)))
              call D3dBs_r_Zero_Ends(1,real_mb(tmp1(1)))


*             **** add -(Vij)*psi_r to Hpsi_r ***
              call D3dBs_rr_Sub2(1,real_mb(tmp1(1)),Hpsi_r(index1))

              !**** include off diagonal terms ****
              if (i.ne.j) then
!$OMP MASTER
                 ehfx = ehfx - eh
                 phfx = phfx - ph
                 dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)
     >                                      - eh
!$OMP END MASTER
*                **** generate (Vij)*psi_r ***
                 call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                           psi_r(index1),
     >                           real_mb(tmp2(1)))
                 call D3dBs_r_Zero_Ends(1,real_mb(tmp2(1)))

*                **** add -(Vij)*psi_r to Hpsi_r ***
                 call D3dBs_rr_Sub2(1,real_mb(tmp2(1)),
     >                               Hpsi_r(index2))
              end if
           end if
           end if
           k1 = k1 + 1
         end do
        end do
        call Parallel_Vector_SumAll(norbs(ms),
     >                        dbl_mb(ehfx_orb(1,ms)))
        end do

        end if

        value =           BA_pop_stack(kcompute(2))
        value = value.and.BA_pop_stack(tmp2(2))
        value = value.and.BA_pop_stack(tmp1(2))
        value = value.and.BA_pop_stack(vij(2))
        value = value.and.BA_pop_stack(dn(2))
        if (.not. value) 
     >    call errquit('pspw_potential_HFX:popping stack memory',0,
     &       MA_ERR)

         call Parallel_SumAll(ehfx)
         call Parallel_SumAll(phfx)

      end if

      return
      end


*     *****************************
*     *                           *
*     *    pspw_energy_sHFX_sub   *
*     *                           *
*     *****************************
      subroutine pspw_energy_sHFX_sub(ispin0,psi_r,ehfx_out,phfx_out)
      implicit none
      integer ispin0
      real  psi_r(*)
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n1,n2,n3,ms,k1
      integer dn(2),tmp1(2),index1,index2
      real    scal1,scal2,dv,seh
      real*8  eh,ph,ss

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      logical  pspw_hfx_localize_closeenough
      external pspw_hfx_localize_closeenough
      real*8   pspw_hfx_localize_switchr
      external pspw_hfx_localize_switchr
      real     coulomb_screened_s_e
      external coulomb_screened_s_e

      
      if (((norbs(1)+norbs(2)).ne.0).and.(.not.relaxed)) then
!$OMP MASTER
        ehfx = 0.0d0
        phfx = 0.0d0
!$OMP END MASTER

        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_real,(2*n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call Parallel_shared_vector_szero(.false.,
     >                                     2*n2ft3d,real_mb(dn(1)))
        call Parallel_shared_vector_szero(.true.,
     >                                    n2ft3d,real_mb(tmp1(1)))

        scal1 = real(1.0d0/dble(n1*n2*n3))
        scal2 = real(1.0d0/lattice_omega())
        dv = scal1/scal2

        k1 = 1
        do ms=1,ispin
        do i=1,norbs(ms)
!$OMP MASTER
         dbl_mb(ehfx_orb(1,ms)+i-1) = 0.0d0
!$OMP END MASTER
         do j=1,i

            if (mod(k1,npj).eq.taskid_j) then
            if (pspw_hfx_localize_closeenough(
     >                   int_mb(orbital_list(1,ms)+i-1),
     >                   int_mb(orbital_list(1,ms)+j-1))) then

              index1 = (int_mb(orbital_list(1,ms)+i-1)-1)*n2ft3d + 1
              index2 = (int_mb(orbital_list(1,ms)+j-1)-1)*n2ft3d + 1

              ss = pspw_hfx_localize_switchr(
     >             int_mb(orbital_list(1,ms)+i-1),
     >             int_mb(orbital_list(1,ms)+j-1))


*             **** generate dnij ****
              call D3dBs_rr_Mul(1,psi_r(index1),psi_r(index2),
     >                         real_mb(dn(1)))
              call D3dBs_r_SMul1(1,scal2,real_mb(dn(1)))
              call D3dBs_r_Zero_Ends(1,real_mb(dn(1)))

*             ***** screened coulomb solver ****
              if (solver_type.eq.1) then

*               **** generate dng ****
                call D3dBs_r_SMul1(1,scal1,real_mb(dn(1)))
                call D3dBs_rc_pfft3f(1,0,real_mb(dn(1)))
                call Packs_c_pack(0,real_mb(dn(1)))

*               **** get Ecoul energy ****
                eh = dble(coulomb_screened_s_e(real_mb(dn(1))))
       
*             ***** free-space coulomb solver ****
              else
                 call coulomb2_s_v(real_mb(dn(1)),real_mb(tmp1(1)))
                 call D3dBs_rr_dot(1,real_mb(dn(1)),
     >                               real_mb(tmp1(1)),seh)
                 eh = dble(0.5*seh*dv)
              end if

*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
              eh = eh*hfx_parameter*ss
              if (ispin0.eq.1) eh = eh + eh
              ph = 2.0d0*eh

!$OMP MASTER
              ehfx = ehfx - eh
              phfx = phfx - ph
              dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)-eh
!$OMP END MASTER

              !**** include off diagonal terms ****
              if (i.ne.j) then
!$OMP MASTER
                 ehfx = ehfx - eh
                 phfx = phfx - ph
                 dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)
     >                                      - eh
!$OMP END MASTER
              end if

            end if
            end if
            k1 = k1 + 1

         end do
        end do
        call D1dB_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
        end do

        value =           BA_pop_stack(tmp1(2))
        value = value.and.BA_pop_stack(dn(2))
        if (.not. value) 
     >     call errquit('pspw_energy_HFX_sub:popping stack memory',0,
     &       MA_ERR)

        call D1dB_SumAll(ehfx)
        call D1dB_SumAll(phfx)
      end if
      ehfx_out = ehfx
      phfx_out = phfx

      return
      end




*     ************************************
*     *                    	         *
*     *     pspw_potential_HFX_orb_sub   *
*     *                                  *
*     ************************************
*
*    Note that orb_r and Horb_r are assumed to be replicated rather than psi_r
*    orb_r is not replicated in this routine
*    Horb_r is not reduced in this routine
*
      subroutine pspw_potential_sHFX_orb_sub(ms,psi_r,
     >                                      orb_r,Horb_r)
      implicit none
      integer    ms
      real  psi_r(*)
      real  orb_r(*)
      real  Horb_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer j,n1,n2,n3,q2,p2
      integer dn(2),vij(2),tmp1(2),tmp2(2),index2
      real    scal1,scal2,dv,seh
      real*8  eh,ph

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      real     coulomb_screened_s_e
      external coulomb_screened_s_e

      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then
        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_real,(n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'vij_hfx',vij(2),vij(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'tmp2_hfx',tmp2(2),tmp2(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call Parallel_shared_vector_szero(.false.,n2ft3d,
     >                                    real_mb(dn(1)))
        call Parallel_shared_vector_szero(.false.,n2ft3d,
     >                                    real_mb(vij(1)))
        call Parallel_shared_vector_szero(.false.,n2ft3d,
     >                                    real_mb(tmp1(1)))
        call Parallel_shared_vector_szero(.true.,n2ft3d,
     >                                    real_mb(tmp2(1)))

        scal1 = real(1.0d0/dble(n1*n2*n3))
        scal2 = real(1.0d0/lattice_omega())
        dv = scal1/scal2

        do j=1,norbs(ms)
           call Dneall_ntoqp(int_mb(orbital_list(1,ms)+j-1),q2,p2)
           index2 = (q2-1)*n2ft3d + 1

           if (p2.eq.taskid_j) then
*             **** generate dnij for Vij  ****
              call D3dBs_rr_Mul(1,psi_r(index2),orb_r,real_mb(dn(1)))
              call D3dBs_r_SMul1(1,scal2,real_mb(dn(1)))
              call D3dBs_r_Zero_Ends(1,real_mb(dn(1)))
   
*             ***** screened coulomb solver ****
              if (solver_type.eq.1) then
                call D3dBs_r_SMul1(1,scal1,real_mb(dn(1)))
                call D3dBs_rc_pfft3f(1,0,real_mb(dn(1)))
                call Packs_c_pack(0,real_mb(dn(1)))

*               **** get Ecoul energy ****
                eh = dble(coulomb_screened_s_e(real_mb(dn(1))))

*               **** generate Vcoul ****
                call coulomb_screened_s_v(real_mb(dn(1)),
     >                                    real_mb(vij(1)))
                call Packs_c_unpack(0,real_mb(vij(1)))
                call D3dBs_cr_pfft3b(1,0,real_mb(vij(1)))

*             ***** free-space coulomb solver ****
              else
                 call coulomb2_s_v(real_mb(dn(1)),real_mb(vij(1)))
                 call D3dBs_rr_dot(1,real_mb(dn(1)),
     >                               real_mb(vij(1)),seh)
                 eh = dble(0.5*seh*dv)
              end if

*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
              eh = eh*hfx_parameter
              call D3dBs_r_SMul1(1,hfx_parameter,real_mb(vij(1)))
              if (ispin.eq.1) eh = eh + eh
              ph = 2.0d0*eh


*             **** generate (Vij)*psi_r ***
              call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                           psi_r(index2),
     >                           real_mb(tmp1(1)))
              call D3dBs_r_Zero_Ends(1,real_mb(tmp1(1)))

*             **** add -(Vij)*psi_r to Hpsi_r ***
              call D3dBs_rr_Sub2(1,real_mb(tmp1(1)),Horb_r)
           end if
        end do
        value =           BA_pop_stack(tmp2(2))
        value = value.and.BA_pop_stack(tmp1(2))
        value = value.and.BA_pop_stack(vij(2))
        value = value.and.BA_pop_stack(dn(2))
        if (.not. value) 
     >  call errquit('pspw_potential_sHFX_orb_sub:popping stack memory',
     >               0,MA_ERR)

c        **** eh and ph not used yet ****
c        call D1dB_SumAll(eh)
c        call D1dB_SumAll(ph)
      end if
      return
      end



*     ***********************************
*     *                    	        *
*     *    pspw_potential_HFX_sub2      *
*     *                                 *
*     ***********************************
*
*   This routine is a kernel for computing exact exchange.
*
*    for i=istart:iend
*    for j=jstart:jend
*       dnij(*) = psi_r(*,j) .* psi_r(*,i)
*       Vij(*)  = Coulomb operator(dnij(*))
*       Hpsi_r(*,i) = Vij(*) .* psi_r(*,j)
*       Hpsi_r(*,j) = Vij(*) .* psi_r(*,i)
*       ehfx += 0.5*<psi_r(*,i)|Hpsi(*,i)> 
*             + 0.5*<psi_r(*,j)|Hpsi(*,j)>
*
*   Entry - solver_type: if solver_type==1 then periodic solver, else aperiodic solver
*           istart,iend: indexes
*           jstart,jend: indexes
*           imodn,imodtask: used to define which (i,j) combinations are computed.
*           n2ft3d: size of realspace grid
*           psi_r: wavenfucntions in realspace.
*           ehfx: running sum of exchange energy, not initialized in this routine
*
*   Exit - Hpsi_r: wavefunction gradients in realspace.
*          ehfx: running sum of exchange energy.

      subroutine pspw_potential_sHFX_sub2(solver_type,
     >                                   istart,iend,
     >                                   jstart,jend,
     >                                   imodn,imodtask,
     >                                   n2ft3d,psi_r,Hpsi_r,
     >                                   ehfx)
      implicit none
      integer solver_type
      integer istart,iend,jstart,jend
      integer imodn,imodtask
      integer n2ft3d
      real  psi_r(n2ft3d,*)
      real  Hpsi_r(n2ft3d,*)
      real*8  ehfx

#include "bafdecls.fh"
#include "errquit.fh"

      integer taskid_j

*     **** local variables ****
      logical value,done,special
      integer n1,n2,n3
      integer dn(2),vij(2),tmp1(2)
      integer i1,j1,k1,NN
      integer i2,j2,k2
      integer i3,j3,k3
      real  scal1,scal2,dv,seh
      real*8  eh,ph

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      logical  D3dBs_rc_pfft3_queue_filled,D3dBs_cr_pfft3_queue_filled
      external D3dBs_rc_pfft3_queue_filled,D3dBs_cr_pfft3_queue_filled
      real     icoulomb_screened_s_e
      external icoulomb_screened_s_e

      call Parallel2d_taskid_j(taskid_j)

      special = ((istart.eq.jstart).and.(iend.eq.jend))

      call D3dB_nx(1,n1)
      call D3dB_ny(1,n2)
      call D3dB_nz(1,n3)
      value = BA_push_get(mt_real,(n2ft3d),'dn_hfx',dn(2),dn(1))
      value = value.and.
     >        BA_push_get(mt_real,(n2ft3d),'vij_hfx',vij(2),vij(1))
      value = value.and.
     >        BA_push_get(mt_real,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('pspw_potential_sHFX_sub:out of stack',0,MA_ERR)
      call Parallel_shared_vector_szero(.false.,
     >                                  n2ft3d,real_mb(dn(1)))
      call Parallel_shared_vector_szero(.false.,
     >                                  n2ft3d,real_mb(vij(1)))
      call Parallel_shared_vector_szero(.true., 
     >                                  n2ft3d,real_mb(tmp1(1)))

      scal1 = real(1.0d0/dble(n1*n2*n3))
      scal2 = real(1.0d0/lattice_omega())
      dv = scal1/scal2

*     *** special if i and j span the same indexes ***
      if (special) then
         NN = (iend-istart+1)*(jend-jstart+2)/2
      else
         NN = (iend-istart+1)*(jend-jstart+1)
      end if

*     ***** screened coulomb solver ****
      if (solver_type.eq.1) then
        i1 = istart
        j1 = jstart
        k1 = 1

        i2 = istart
        j2 = jstart
        k2 = 1

        i3 = istart
        j3 = jstart
        k3 = 1
        done = .false.
        do while (.not.done)

*          *** pipeline step 1 ***
           if (k1.le.NN) then

              if (mod(k1,imodn).eq.imodtask) then

*                **** generate dnij for Vij  ****
                 call D3dBs_rr_Mul(1,psi_r(1,j1),
     >                              psi_r(1,i1),real_mb(dn(1)))
                 call D3dBs_r_SMul1(1,scal2*scal1,real_mb(dn(1)))
                 call D3dBs_r_Zero_Ends(1,real_mb(dn(1)))
                 call D3dBs_rc_pfft3f_queuein(0,real_mb(dn(1)))

              end if

              k1 = k1 + 1
              j1 = j1 + 1
              if (special) then
                 if (j1.gt.i1) then
                    j1 = jstart
                    i1 = i1 + 1
                 end if
              else
                 if (j1.gt.jend) then
                    j1 = jstart
                    i1 = i1 + 1
                 end if
              end if
           end if

*          *** pipeline step 2 ***
           if (     ((D3dBs_rc_pfft3_queue_filled()).or.(k1.gt.NN))
     >         .and.(k2.le.NN)) then

              if (mod(k2,imodn).eq.imodtask) then
                 call D3dBs_rc_pfft3f_queueout(0,real_mb(dn(1)))

*                **** generate Vcoul ****
                 eh = dble(icoulomb_screened_s_e(real_mb(dn(1))))
                 call coulomb_screened_s_v(real_mb(dn(1)),
     >                                     real_mb(vij(1)))


*                **** calculcate ph ****
!OMP MASTER
                 ehfx = ehfx - eh
!OMP END MASTER

*                **** include transpose ***
                 if ((i2.ne.j2).or.(.not.special)) then
!OMP MASTER
                    ehfx = ehfx - eh
!OMP END MASTER
                 end if

                 call D3dBs_cr_pfft3b_queuein(0,real_mb(vij(1)))
              end if

              k2 = k2 + 1
              j2 = j2 + 1
              if (special) then
                 if (j2.gt.i2) then
                    j2 = jstart
                    i2 = i2 + 1
                 end if
              else
                 if (j2.gt.jend) then
                    j2 = jstart
                    i2 = i2 + 1
                 end if
              end if

           end if

*          *** pipeline step 3 ***
           if ((D3dBs_cr_pfft3_queue_filled()).or.(k2.gt.NN)) then

              if (mod(k3,imodn).eq.imodtask) then
                 call D3dBs_cr_pfft3b_queueout(0,real_mb(vij(1)))

*                **** generate (Vij)*psi_r ***
                 call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                              psi_r(1,j3),
     >                              real_mb(tmp1(1)))
                 call D3dBs_r_Zero_Ends(1,real_mb(tmp1(1)))

*                **** add -(Vij)*psi_r to Hpsi_r ***
                 call D3dBs_rr_Sub2(1,real_mb(tmp1(1)),Hpsi_r(1,i3))

                 !**** include transpose ****
                 if ((i3.ne.j3).or.(.not.special)) then

*                   **** generate (Vij)*psi_r ***
                    call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                              psi_r(1,i3),
     >                              real_mb(tmp1(1)))
                    call D3dBs_r_Zero_Ends(1,real_mb(tmp1(1)))

*                   **** add -(Vij)*psi_r to Hpsi_r ***
                    call D3dBs_rr_Sub2(1,real_mb(tmp1(1)),
     >                                Hpsi_r(1,j3))
                 end if
              endif

              k3 = k3 + 1
              j3 = j3 + 1
              if (special) then
                 if (j3.gt.i3) then
                    j3 = jstart
                    i3 = i3 + 1
                 end if
              else
                 if (j3.gt.jend) then
                    j3 = jstart
                    i3 = i3 + 1
                 end if
              end if
                
           end if
           done = (k1.gt.NN).and.(k2.gt.NN).and.(k3.gt.NN)
        end do !**** while ****
       

*     ***** free-space coulomb solver -- not pipelined ****
      else
         k1 = 1
         i1 = istart
         j1 = jstart
         done = .false.
         do while (.not.done)

            if (mod(k1,imodn).eq.imodtask) then

*              **** generate dnij for Vij  ****
               call D3dBs_rr_Mul(1,psi_r(1,j1),psi_r(1,i1),
     >                          real_mb(dn(1)))
               call D3dBs_r_SMul1(1,scal2,real_mb(dn(1)))
               call D3dBs_r_Zero_Ends(1,real_mb(dn(1)))
   
               call coulomb2_s_v(real_mb(dn(1)),real_mb(vij(1)))
               call D3dBs_rr_idot(1,real_mb(dn(1)),
     >                              real_mb(vij(1)),seh)
               eh = dble(0.5*seh*dv)

*              **** calculcate ph ****
!OMP MASTER
               ehfx = ehfx - eh
!OMP END MASTER

*              **** generate (Vij)*psi_r ***
               call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                            psi_r(1,j1),
     >                            real_mb(tmp1(1)))
               call D3dBs_r_Zero_Ends(1,real_mb(tmp1(1)))

*              **** add -(Vij)*psi_r to Hpsi_r ***
               call D3dBs_rr_Sub2(1,real_mb(tmp1(1)),Hpsi_r(1,i1))
 
               !**** include transpose terms ****
               if ((i1.ne.j1).or.(.not.special)) then
!OMP MASTER
                  ehfx = ehfx - eh
!OMP END MASTER
 
*                 **** generate (Vij)*psi_r ***
                  call D3dBs_rr_Mul(1,real_mb(vij(1)),
     >                            psi_r(1,i1),
     >                            real_mb(tmp1(1)))
                  call D3dBs_r_Zero_Ends(1,real_mb(tmp1(1)))
 
*                 **** add -(Vij)*psi_r to Hpsi_r ***
                  call D3dBs_rr_Sub2(1,real_mb(tmp1(1)),
     >                                Hpsi_r(1,j1))
               end if

            end if
            
            k1 = k1 + 1
            j1 = j1 + 1
            if (special) then
               if (j1.gt.i1) then
                  j1 = jstart
                  i1 = i1 + 1
               end if
            else
               if (j1.gt.jend) then
                  j1 = jstart
                  i1 = i1 + 1
               end if
            end if
            done = (k1.gt.NN)
         end do

      end if

      !**** deallocate memory ****
      value =           BA_pop_stack(tmp1(2))
      value = value.and.BA_pop_stack(vij(2))
      value = value.and.BA_pop_stack(dn(2))
      if (.not. value) call errquit(
     >   'pspw_potential_sHFX_sub2:popping stack memory',0,MA_ERR)

      return
      end



*     *******************************************
*     *                    	                *
*     *     pspw_potential_sHFX_orb_replicated   *
*     *                                         *
*     *******************************************
      subroutine pspw_potential_sHFX_orb_replicated(ms,psi_r,
     >                                  orb_r,Horb_r)
      implicit none
      integer    ms
      real*8     psi_r(*)
      real*8     orb_r(*)
      real*8     Horb_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"


      call nwpw_timing_start(33)
      if ((norbs(ms).ne.0).and.relaxed) then
         if (replicated) then
            call Parallel_shared_vector_scopy(.true.,n2ft3d,
     >                                  real_mb(Hpsi_r_replicated(1)))
            call pspw_potential_HFX_orb_sub(ms,psi_r,orb_r,
     >                                  real_mb(Hpsi_r_replicated(1)))
            call D1dBs_Vector_SumAll(n2ft3d,
     >                               real_mb(Hpsi_r_replicated(1)))
            call SAXPY_OMP(n2ft3d,1.0d0,real_mb(Hpsi_r_replicated(1)),1,
     >                 Horb_r,1)

         else
            call pspw_potential_sHFX_orb_sub(ms,psi_r,orb_r,Horb_r)
         end if

      end if
      call nwpw_timing_end(33)
      return
      end






*    ***** routines below used for pspw_et calculations ****


*     *************************
*     *                       *
*     *    pspw_energy_sHFX2  *
*     *                       *
*     *************************
      subroutine pspw_energy_sHFX2(ispin0,psi1_r,psi2_r,
     >                            ehfx_out,phfx_out)
      implicit none
      integer ispin0
      real   psi1_r(*),psi2_r(*)
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

      integer q,n,indx1,indx2

      call nwpw_timing_start(33)

c     **** calculate HFX energy  ****
      if ((norbs(1)+norbs(2)).ne.0) then

         if (replicated) then
            call Parallel_shared_vector_szero(.false.,nrsize,
     >                                 real_mb(psi_r_replicated(1)))
            call Parallel_shared_vector_szero(.true.,nrsize,
     >                                 real_mb(Hpsi_r_replicated(1)))
            do q=1,neqall
               call Dneall_qton(q,n)
               indx1 = (q-1)*n2ft3d + 1
               indx2 = psi_r_replicated(1)+(n-1)*n2ft3d
               call Parallel_shared_vector_scopy(.false.,n2ft3d,
     >                                psi1_r(indx1),real_mb(indx2))

               indx2 = Hpsi_r_replicated(1)+(n-1)*n2ft3d
               call Parallel_shared_vector_scopy(.true.,n2ft3d,
     >                            psi2_r(indx1),real_mb(indx2))
            end do
            call D1dBs_Vector_SumAll(nrsize,
     >                               real_mb(psi_r_replicated(1)))
            call D1dBs_Vector_SumAll(nrsize,
     >                              real_mb(Hpsi_r_replicated(1)))
            call pspw_energy_sHFX2_sub(ispin0,
     >                               real_mb(psi_r_replicated(1)),
     >                               real_mb(Hpsi_r_replicated(1)),
     >                               ehfx_out,phfx_out)

         else
            call pspw_energy_sHFX2_sub(ispin0,psi1_r,psi2_r,
     >                                ehfx_out,phfx_out)
         end if

c     **** nothing to do ****
      else
         ehfx_out = ehfx
         phfx_out = phfx
      end if
      call nwpw_timing_end(33)

      return
      end



*     *****************************
*     *                           *
*     *    pspw_energy_sHFX2_sub  *
*     *                           *
*     *****************************
      subroutine pspw_energy_sHFX2_sub(ispin0,psi1_r,psi2_r,
     >                               ehfx_out,phfx_out)
      implicit none
      integer ispin0
      real   psi1_r(*)
      real   psi2_r(*)
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n1,n2,n3,ms,k1
      integer dn(2),tmp1(2),index1,index2
      real    scal1,scal2,dv,seh
      real*8  eh,ph

*     **** external functions ****
      real*8   lattice_omega,coulomb_screened_s_e
      external lattice_omega,coulomb_screened_s_e

      integer taskid

      call Parallel_taskid(taskid)

      
      if ((norbs(1)+norbs(2)).ne.0) then
!OMP MASTER
        ehfx = 0.0d0
        phfx = 0.0d0
!OMP END MASTER

        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_real,(2*n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_real,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call Parallel_shared_vector_szero(.false.,
     >                                    2*n2ft3d,real_mb(dn(1)),1)
        call Parallel_shared_vector_szero(.true.,
     >                                    n2ft3d,real_mb(tmp1(1)),1)

        scal1 = real(1.0d0/dble(n1*n2*n3))
        scal2 = real(1.0d0/lattice_omega())
        dv = scal1/scal2

        k1 = 1
        do ms=1,ispin
        do i=1,norbs(ms)
!OMP MASTER
         dbl_mb(ehfx_orb(1,ms)+i-1) = 0.0d0
!OMP END MASTER
         do j=1,i

            if (mod(k1,npj).eq.taskid_j) then
              index1 = (int_mb(orbital_list(1,ms)+i-1)-1)*n2ft3d + 1
              index2 = (int_mb(orbital_list(1,ms)+j-1)-1)*n2ft3d + 1

*             **** generate dnij ****
              call D3dBs_rr_Mul(1,psi1_r(index1),psi2_r(index2),
     >                         real_mb(dn(1)))
              call D3dBs_r_SMul1(1,scal2,real_mb(dn(1)))
              call D3dBs_r_Zero_Ends(1,real_mb(dn(1)))

*             **** generate  dnji****
              call D3dBs_rr_Mul(1,psi1_r(index2),psi2_r(index1),
     >                         real_mb(dn(1)+n2ft3d))
              call D3dBs_r_SMul1(1,scal2,real_mb(dn(1)+n2ft3d))
              call D3dBs_r_Zero_Ends(1,real_mb(dn(1)+n2ft3d))

*             ***** screened coulomb solver ****
              if (solver_type.eq.1) then

*               **** generate dng ****
                call D3dBs_r_SMul1(1,scal1,real_mb(dn(1)))
                call D3dBs_rc_pfft3f(1,0,real_mb(dn(1)))
                call Packs_c_pack(0,real_mb(dn(1)))

                call D3dB_r_SMul1(1,scal1,real_mb(dn(1)+n2ft3d))
                call D3dB_rc_pfft3f(1,0,real_mb(dn(1)+n2ft3d))
                call Packs_c_pack(0,real_mb(dn(1)+n2ft3d))

*               **** get Ecoul energy ****
                call coulomb_screened_s_v(real_mb(dn(1)),
     >                                    real_mb(tmp1(1)))
                call Packs_cc_dot(0,real_mb(dn(1)+n2ft3d),
     >                              real_mb(tmp1(1)),seh)
                eh = dble(0.5*seh*lattice_omega())
       
*             ***** free-space coulomb solver ****
              else
                 call coulomb2_s_v(real_mb(dn(1)),real_mb(tmp1(1)))
                 call D3dBs_rr_dot(1,real_mb(dn(1)+n2ft3d),
     >                            real_mb(tmp1(1)),seh)
                 eh = dble(0.5*seh*dv)

                 !if (taskid.eq.0) write(*,*) " - ms,i,j,ehfx=",ms,i,j,eh

              end if

*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
              eh = eh*hfx_parameter
              if (ispin0.eq.1) eh = eh + eh
              ph = 2.0d0*eh

!OMP MASTER
              ehfx = ehfx - eh
              phfx = phfx - ph
              dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)-eh
!OMP END MASTER

              !**** include off diagonal terms ****
              if (i.ne.j) then
!OMP MASTER
                 ehfx = ehfx - eh
                 phfx = phfx - ph
                 dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)
     >                                      - eh
!OMP END MASTER
              end if

            end if
            k1 = k1 + 1

         end do
        end do
        call D1dB_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
        end do

        value =           BA_pop_stack(tmp1(2))
        value = value.and.BA_pop_stack(dn(2))
        if (.not. value) 
     >     call errquit('pspw_energy_sHFX2_sub:popping stack memory',0,
     >       MA_ERR)

        call D1dB_SumAll(ehfx)
        call D1dB_SumAll(phfx)
      end if
      ehfx_out = ehfx
      phfx_out = phfx

      return
      end




