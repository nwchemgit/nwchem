*
* $Id$
*

*     *************************
*     *                       *
*     *     pspw_init_HFX     *
*     *                       *
*     *************************
      subroutine pspw_init_HFX(rtdb,ispin0,ne)
      implicit none
      integer rtdb
      integer ispin0
      integer ne(2)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "pspw_hfx.fh"

*     **** local variables ****
      logical value
      integer ma_type
      integer n1,n2,n3,mapping,ms,neq(2)

*     **** external functions ****
      logical  control_has_ngrid_small
      external control_has_ngrid_small
      integer  control_version,control_mapping,Butter_size
      integer  Dneall_na_ptr
      external control_version,control_mapping,Butter_size
      external Dneall_na_ptr

      ispin = ispin0
      norbs(1) = 0
      norbs(2) = 0
      ehfx = 0.0d0
      hfx_on = .false.
      hfx_virtual_on = .true.
      localize_on  = .false.
      localize2_on = .false.
      has_smallgrid = control_has_ngrid_small()
      call D3dB_n2ft3d(1,n2ft3d)

      if (.not.btdb_get(rtdb,'pspw:HFX',mt_log,1,hfx_on))
     >   hfx_on = .false.

      if (.not.btdb_get(rtdb,'pspw:HFX_Lin',mt_log,1,lin_on)) 
     >   lin_on = .false.
      lin_run = .false.
      lin_off = .not.lin_on
      lin_count1 = 0
      lin_count2 = 0

*     **** get the number of HFX orbitals ****
      if (hfx_on) then
         !hfx_on = .true.
         do ms=1,ispin
           norbs(ms) = ne(ms)
           if(.not.BA_alloc_get(mt_int,norbs(ms),
     >       'orbital_list',orbital_list(2,ms),orbital_list(1,ms)))
     >        call errquit('pspw_init_HFX:out of heap memory',0,MA_ERR)
          
            do n1=1,norbs(ms)
               int_mb(orbital_list(1,ms)+n1-1) = n1 + (ms-1)*ne(1)
            end do
         end do

      else if (rtdb_ma_get(rtdb, 'pspw:HFX_up_orbitals', ma_type,
     >                norbs(1), orbital_list(2,1))) then

            if (.not.BA_get_index(orbital_list(2,1),orbital_list(1,1))) 
     >        call errquit(
     >       'pspw_init_HFX: ma_get_index failed for actlist',911,
     >        MA_ERR)

         if (rtdb_ma_get(rtdb, 'pspw:HFX_down_orbitals', ma_type,
     >                norbs(2), orbital_list(2,2))) then

           if (.not.BA_get_index(orbital_list(2,2),orbital_list(1,2))) 
     >        call errquit(
     >       'pspw_init_HFX: ma_get_index failed for actlist',911,
     >         MA_ERR)
         end if

         hfx_on = .true.

      end if


      if (hfx_on) then

         if (.not. btdb_get(rtdb,
     >                      'pspw:HFX_screening_radius',
     >                      mt_dbl,1,rcut)) 
     >       rcut = 8.0d0

         if (.not. btdb_get(rtdb,
     >                      'pspw:HFX_screening_power',
     >                      mt_dbl,1,pp)) 
     >       pp = 8.0d0

         if (.not. btdb_get(rtdb,
     >                      'pspw:HFX_screening_type',
     >                      mt_int,1,flag)) 
     >       flag = 0

         if (.not. btdb_get(rtdb,
     >                      'pspw:HFX_relax',
     >                      mt_log,1,relaxed)) 
     >       relaxed = .true.

         if (.not. btdb_get(rtdb,
     >                      'pspw:HFX_solver_type',
     >                      mt_int,1,solver_type)) then

            if (control_version().eq.3) solver_type = 1
            if (control_version().eq.4) solver_type = 2
         end if

         if (.not. btdb_get(rtdb,
     >                      'pspw:HFX_parameter',
     >                       mt_dbl,1,HFX_parameter))
     >       HFX_parameter = 1.0d0

         if (.not. btdb_get(rtdb,
     >                      'pspw:HFX_print_orbital_contribution',
     >                       mt_log,1,orb_contribution))
     >       orb_contribution = .false.

*        **** get localization parameters ****
         if (.not. btdb_get(rtdb,'pspw:HFX_localize',
     >                      mt_log,1,localize_on)) 
     >      localize_on = .false.

         if (.not. btdb_get(rtdb,'pspw:HFX_localize2',
     >                      mt_log,1,localize2_on)) 
     >      localize2_on = .false.
         if (localize2_on) localize_on = .true.

         if (.not. btdb_get(rtdb,'pspw:HFX_localize2:maxiter',
     >                      mt_dbl,1,maxiter)) 
     >      maxiter = 5

         if (.not. btdb_get(rtdb,'pspw:HFX_localize:prune_algorithm',
     >                      mt_int,1,prune_algorithm)) 
     >      prune_algorithm = 0

         if (localize_on) then
            if (.not. btdb_get(rtdb,'pspw:HFX_localize_rin',
     >                         mt_dbl,1,rin_thresh)) then
               rin_thresh = 6.5d0
               if (prune_algorithm.eq.0) rin_thresh = 5.0d-3
            end if
            if (.not. btdb_get(rtdb,'pspw:HFX_localize_rout',
     >                         mt_dbl,1,rout_thresh)) then
               rout_thresh = 7.5d0
               if (prune_algorithm.eq.0) rout_thresh = 1.0d-2
            end if
         end if

 

*        **** initialize coulomb_screened ****
         if (solver_type.eq.1) then
              call coulomb_screened_init(flag,rcut,pp)
              if (localize_on.and.has_smallgrid) then
                 call coulomb_screened_small_init(flag,rcut,pp)
                 call D3dB_r_1to3_init()
                 call D3dB_r_3to1_init()
              end if

*        **** initialize free-space coulomb if necessary ****
         else
            if (control_version().eq.3) then
               call D3dB_nx(1,n1)
               call D3dB_ny(1,n2)
               call D3dB_nz(1,n3)
               mapping = control_mapping()
               call D3dB_Init(2,2*n1,2*n2,2*n3,mapping)
               call coulomb2_init()
            end if

         end if

*        **** initialize orb_contribution ****
         do ms=1,ispin
           value = BA_alloc_get(mt_dbl,norbs(ms),
     >                'ehfx_orb',ehfx_orb(2,ms),ehfx_orb(1,ms))
           if (.not. value)   
     >       call errquit('pspw_init_HFX: out of heap memory',1, MA_ERR)
         end do

      end if


c     **** define Lin_projectors  ****
      if (lin_on) then
        call Dneall_neq(neq)
        neqall = neq(1)+neq(2)
        value = BA_alloc_get(mt_dbl,neqall*n2ft3d,
     >                      'Lin_projectors',
     >                       Lin_projectors(2),
     >                       Lin_projectors(1))
        if (.not. value)   
     >    call errquit('pspw_init_HFX: out of heap memory',3,MA_ERR)
      end if


c     **** define extra psi and Hpsi  ****
      call Parallel2d_np_j(npj)
      call Parallel2d_taskid_j(taskid_j)   
      replicated = (npj.gt.1)
      if (hfx_on.and.replicated) then

         if (.not.btdb_get(rtdb,'pspw:HFX_butter',mt_log,1,butterfly))
     >      butterfly = .false.

        call Dneall_neq(neq)
        neqall = neq(1)+neq(2)

        if (butterfly) then
           nrsize = n2ft3d*Butter_size(taskid_j,npj,
     >                                  int_mb(Dneall_na_ptr(1)))
        else
           nrsize = (ne(1)+ne(2))*n2ft3d
        end if

        value = BA_alloc_get(mt_dbl,nrsize,
     >                      'psi_r_replicated',
     >                       psi_r_replicated(2),
     >                       psi_r_replicated(1))
        value = value.and.
     >          BA_alloc_get(mt_dbl,nrsize,
     >                      'Hpsi_r_replicated',
     >                       Hpsi_r_replicated(2),
     >                       Hpsi_r_replicated(1))
        if (.not. value)   
     >    call errquit('pspw_init_HFX: out of heap memory',4,MA_ERR)
         call dcopy(nrsize,0.0d0,0,dbl_mb(psi_r_replicated(1)),1)
         call dcopy(nrsize,0.0d0,0,dbl_mb(Hpsi_r_replicated(1)),1)
      end if

      if (localize_on)  call pspw_hfx_localize_init()
   
      return
      end


*     *************************
*     *                       *
*     *     pspw_end_HFX      *
*     *                       *
*     *************************
      subroutine pspw_end_HFX()
      implicit none

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter(MASTER=0)
      logical value
      integer i,ms

*     **** external functions ****
      integer  control_version
      external control_version

      if ((norbs(1)+norbs(2)).gt.0) then

*       **** print out orbital contributions ****
        if (orb_contribution) then
           call Parallel_taskid(taskid)
           if (taskid.eq.MASTER) then
              write(6,487)
              write(6,488)
              do ms=1,ispin
              do i=1,norbs(ms)
                write(6,489) 
     >            ms,
     >            int_mb(orbital_list(1,ms)+i-1),
     >            dbl_mb(ehfx_orb(1,ms)+i-1)
              end do
              end do
           end if
  487   format(//,'== Orbital Contributions to HFX ==')
  488   format(/1x,'orbital',15x,
     >         'HF_Exchange') 
  489   format(1x,i3,i7,2x,e18.6)
        end if


*       **** deallocate memory ****
        if (localize_on) call pspw_hfx_localize_end()
        value = .true.
        do ms=1,ispin
          value = value.and.BA_free_heap(orbital_list(2,ms)) 
          value = value.and.BA_free_heap(ehfx_orb(2,ms))
        end do
c        if (fractional) 
c     >     value = value.and.BA_free_heap(amatrix(2))
c        if (.not. value)
c     >  call errquit('pspw_end_HFX:error freeing heap memory',0, MA_ERR)


*        **** end coulomb_screened ****
        if (solver_type.eq.1) then
          call coulomb_screened_end()
          if (localize_on.and.has_smallgrid) then
             call coulomb_screened_small_end()
             call D3dB_r_1to3_end()
             call D3dB_r_3to1_end()
          end if

*        **** end free-space coulomb if necessary ****
        else
           if (control_version().eq.3) then
              call coulomb2_end()
              call D3dB_end(2)
           end if
        end if

*        **** deallocate Lin_projectors ****
        if (lin_on) then
           value = value.and.BA_free_heap(Lin_projectors(2))
           if (.not. value)
     >     call errquit('pspw_end_HFX:error freeing heap memory',
     >                  0,MA_ERR)

        end if


*        **** deallocate replicated space if necessary ****
        if (replicated) then
           value = value.and.BA_free_heap(psi_r_replicated(2))
           value = value.and.BA_free_heap(Hpsi_r_replicated(2))
           if (.not. value)
     >     call errquit('pspw_end_HFX:error freeing heap memory',
     >                  1,MA_ERR)
        end if

      end if

      return
      end

*     **********************************
*     *                                *
*     *     pspw_print_HFX_counters    *
*     *                                *
*     **********************************
      subroutine pspw_print_HFX_counters(unit)
      implicit none
      integer unit

#include "bafdecls.fh"
#include "pspw_hfx.fh"
      
      if (hfx_on) then
         write(unit,501)
         write(unit,502) lin_count2 
         write(unit,503) lin_count1
         write(unit,*) 
         write(unit,*) 
      end if

      return
 501  format("== HFX Counters ==")
 502  format("  HFX Evalulations         =" I6)
 503  format("  Lin ACE HFX Evalulations =" I6)
      end

*     *************************
*     *                       *
*     *     pspw_print_HFX    *
*     *                       *
*     *************************
      subroutine pspw_print_HFX(unit)
      implicit none
      integer unit

#include "bafdecls.fh"
#include "pspw_hfx.fh"

*     **** local variables ****
      integer i,ms
      real*8   control_attenuation
      external control_attenuation

      if (hfx_on) then
        if (relaxed) then
          write(unit,1001)
        else
          write(unit,1002)
        end if
        if (ispin.eq.1) then
         write(unit,1003) (int_mb(orbital_list(1,1)+i-1),i=1,norbs(1))
        else
         write(unit,1004) (int_mb(orbital_list(1,1)+i-1),i=1,norbs(1))
         write(unit,1005) (int_mb(orbital_list(1,2)+i-1),i=1,norbs(2))
        end if

        if (solver_type.eq.1) then
          write(unit,1006)
           if (rcut.ge.0.0d0) write(unit,1008) rcut
           if (rcut.ge.0.0d0) write(unit,1009) pp
           if (rcut.ge.0.0d0) write(unit,1011) flag
           if ((rcut.ge.0.0d0).and.(flag.eq.2)) 
     >        write(unit,1012) control_attenuation()
        else
          write(unit,1007)
        end if
        if (hfx_parameter.ne.1.0d0) write(unit,1010) hfx_parameter
        if (localize_on) then
           if (localize2_on) then
              write(unit,1027)
              write(unit,1028) maxiter
           else
              write(unit,1020)
           end if
           if (prune_algorithm.eq.0) then
              write(unit,1021)
              write(unit,1022) rin_thresh
              write(unit,1023) rout_thresh
           else
              write(unit,1024)
              write(unit,1025) rin_thresh
              write(unit,1026) rout_thresh
           end if
        end if
        if (lin_on) write(unit,1029)
        write(unit,*)

      end if

      return
 1001 FORMAT(6x,"- HFX relaxed")
 1002 FORMAT(6x,"- HFX unrelaxed")
 1003 FORMAT(6x,"- HFX restricted orbitals :",10I5)
 1004 FORMAT(6x,"- HFX alpha orbitals:",10I5)
 1005 FORMAT(6x,"- HFX beta orbitals :",10I5)

 1006 FORMAT(6x,"- HFX screened coulomb solver")
 1007 FORMAT(6x,"- HFX free-space coulomb solver")
 1008 FORMAT(6x,"- HFX screening radius(pspw:HFX_screening_radius):",
     >       E10.3)
 1009 FORMAT(6x,"- HFX screening power (pspw:HFX_screening_power) :",
     >       E10.3)
 1010 FORMAT(6x,"- HFX scaling parameter (pspw:HFX_parameter)     :",
     >       E10.3)
 1011 FORMAT(6x,"- HFX screening type (pspw:HFX_screening_type)   :",
     >       I2)
 1012 FORMAT(6x,"- attenuation parameter (nwpw:attenuation)       :",
     >       E10.3)
 1020 FORMAT(6x,"- HFX localized orbital solver")
 1021 FORMAT(6x,"- HFX localized orbital absolute overlap pruning")
 1022 FORMAT(6x,"- Overlap threshold low (pspw:HFX_localize_rin)  :",
     >       E10.3)
 1023 FORMAT(6x,"- Overlap threshold high (pspw:HFX_localize_rout):",
     >       E10.3)
 1024 FORMAT(6x,"- HFX localized orbital distance pruning")
 1025 FORMAT(6x,"- Distance cutoff low  (pspw:HFX_localize_rin)   :",
     >       E10.3)
 1026 FORMAT(6x,"- Distance cutoff high (pspw:HFX_localize_rout)  :",
     >       E10.3)
 1027 FORMAT(6x,"- HFX Wannier orbital solver")
 1028 FORMAT(6x,"- HFX Wannier maxiter(pspw:HFX_localize2:maxiter):",
     >       I3)
 1029 FORMAT(6x,"- HFX Lin extrapolation")
      end


*     ****************************
*     *                    	 *
*     *     pspw_potential_HFX   *
*     *                          *
*     ****************************
      subroutine pspw_potential_HFX(ispin0,psi_r,Hpsi_r)
      implicit none
      integer    ispin0
      real*8     psi_r(*)
      real*8     Hpsi_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

      call nwpw_timing_start(33)


      !**** run Lin projector HFX operator ****
      if ((lin_run).and.(lin_on).and.(.not.lin_off)) then
         call pspw_Lin_proj_HFX(dbl_mb(Lin_projectors(1)),psi_r,Hpsi_r)
         lin_count1 = lin_count1 + 1

      !**** HFX operator ****
      else
         if (localize_on) then
            call pspw_hfx_localize_start(psi_r)
            call pspw_potential_HFX0(ispin0,dbl_mb(psiloc_r(1)),Hpsi_r)
            call pspw_hfx_localize_stop(Hpsi_r)
         else
            call pspw_potential_HFX0(ispin0,psi_r,Hpsi_r)
         end if

         !**** define Lin projectors ****
         if ((lin_on).and.(.not.lin_off)) then
            call pspw_Lin_proj_HFX_set(psi_r,Hpsi_r,
     >                                 dbl_mb(Lin_projectors(1)))
            lin_run = .true.
         end if
         lin_count2 = lin_count2 + 1
 
      end if
      call nwpw_timing_end(33)

      return
      end

*     ****************************
*     *                    	 *
*     *    pspw_potential_HFX0   *
*     *                          *
*     ****************************
      subroutine pspw_potential_HFX0(ispin0,psi_r,Hpsi_r)
      implicit none
      integer    ispin0
      real*8     psi_r(*)
      real*8     Hpsi_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

      integer istart,iend,jstart,jend,imodn,imodtask
      integer ms,l,q,n,indx1,indx2,Levels,neq(2)
      integer requests(5),reqcnt

      integer  Butter_Levels,Dneall_na_ptr
      external Butter_Levels,Dneall_na_ptr



*     ***** now do exchange as normal ****
      ehfx = 0.0d0
      phfx = 0.0d0
      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then

         if (replicated) then

*           **** butterfly algorithm ****
            if (butterfly) then
               call Dneall_neq(neq)
               Levels = Butter_Levels(npj)
               do ms=1,ispin
                  call dcopy(nrsize,0.0d0,0,
     >                       dbl_mb(Hpsi_r_replicated(1)),1)
                  call dcopy(neq(ms)*n2ft3d,
     >                       psi_r(1+(ms-1)*neq(1)*n2ft3d),1,
     >                       dbl_mb(psi_r_replicated(1)),1)

                  do l=0,Levels-1
                     call D1dB_Brdcst_step(l,
     >                       int_mb(Dneall_na_ptr(ms)),-1,
     >                       n2ft3d,
     >                       dbl_mb(psi_r_replicated(1)),
     >                       requests,reqcnt)

                     call Butter_indexes(l,taskid_j,npj,
     >                       int_mb(Dneall_na_ptr(ms)),
     >                       istart,iend,jstart,jend,
     >                       imodn,imodtask)
                     call pspw_potential_HFX_sub2(solver_type,
     >                                  istart,iend,
     >                                  jstart,jend,
     >                                  imodn,imodtask,
     >                                  n2ft3d,
     >                                  dbl_mb(psi_r_replicated(1)),
     >                                  dbl_mb(Hpsi_r_replicated(1)),
     >                                  ehfx)

                     call D1dB_WaitAll(requests,reqcnt)
                  end do

                  call Butter_indexes_L1(taskid_j,npj,
     >                       int_mb(Dneall_na_ptr(ms)),
     >                       istart,iend,jstart,jend,
     >                       imodn,imodtask)
                  if (jend.ge.jstart)
     >               call pspw_potential_HFX_sub2(solver_type,
     >                               istart,iend,
     >                               jstart,jend,
     >                               imodn,imodtask,
     >                               n2ft3d,
     >                               dbl_mb(psi_r_replicated(1)),
     >                               dbl_mb(Hpsi_r_replicated(1)),
     >                               ehfx)
                  call Butter_indexes_L2(taskid_j,npj,
     >                       int_mb(Dneall_na_ptr(ms)),
     >                       istart,iend,jstart,jend,
     >                       imodn,imodtask)
                  call pspw_potential_HFX_sub2(solver_type,
     >                               istart,iend,
     >                               jstart,jend,
     >                               imodn,imodtask,
     >                               n2ft3d,
     >                               dbl_mb(psi_r_replicated(1)),
     >                               dbl_mb(Hpsi_r_replicated(1)),
     >                               ehfx)

                  do l=Levels-1,0,-1
                     call D1dB_Reduce_step(l,
     >                       int_mb(Dneall_na_ptr(ms)),-1,
     >                       n2ft3d,
     >                       dbl_mb(Hpsi_r_replicated(1)),
     >                       dbl_mb(psi_r_replicated(1)))
                  end do
                  call daxpy_omp(neq(ms)*n2ft3d,hfx_parameter,
     >                       dbl_mb(Hpsi_r_replicated(1)),1,
     >                       Hpsi_r(1+(ms-1)*neq(1)*n2ft3d),1)
               end do


*              *** apply hfx_parameter ****
               ehfx = ehfx*hfx_parameter

               if (ispin.eq.1) ehfx = ehfx + ehfx
               call Parallel_SumAll(ehfx)
               phfx = 2.0d0*ehfx

*           **** reduceall algorithm ****
            else
            call dcopy(nrsize,0.0d0,0,dbl_mb(psi_r_replicated(1)),1)
            call dcopy(nrsize,0.0d0,0,dbl_mb(Hpsi_r_replicated(1)),1)
            do q=1,neqall
               call Dneall_qton(q,n)
               indx1 = (q-1)*n2ft3d + 1
               indx2 = psi_r_replicated(1)+(n-1)*n2ft3d
               call dcopy(n2ft3d,psi_r(indx1),1,dbl_mb(indx2),1)
            end do
            call D1dB_Vector_SumAll(nrsize,dbl_mb(psi_r_replicated(1)))
            call pspw_potential_HFX_sub(ispin0,
     >                                  dbl_mb(psi_r_replicated(1)),
     >                                  dbl_mb(Hpsi_r_replicated(1)))
            call D1dB_Vector_SumAll(nrsize,dbl_mb(Hpsi_r_replicated(1)))
            do q=1,neqall
               call Dneall_qton(q,n)
               indx1 = Hpsi_r_replicated(1)+(n-1)*n2ft3d
               indx2 = (q-1)*n2ft3d + 1
               call daxpy_omp(n2ft3d,1.0d0,dbl_mb(indx1),1,
     >                                     Hpsi_r(indx2),1)
            end do
            end if

         else
            call pspw_potential_HFX_sub(ispin0,psi_r,Hpsi_r)
         end if

      end if

      return
      end


*     *************************
*     *                       *
*     *     pspw_energy_HFX   *
*     *                       *
*     *************************
      subroutine pspw_energy_HFX(ispin0,psi_r,ehfx_out,phfx_out)
      implicit none
      integer ispin0
      real*8  psi_r(*)
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "pspw_hfx.fh"

      call nwpw_timing_start(33)
      if (((norbs(1)+norbs(2)).ne.0).and.(.not.relaxed)) then
         if (localize_on) then
            call pspw_hfx_localize_start(psi_r)
            call pspw_energy_HFX0(ispin0,dbl_mb(psiloc_r(1)),
     >                            ehfx_out,phfx_out)
         else
            call pspw_energy_HFX0(ispin0,psi_r,ehfx_out,phfx_out)
         end if
      else
         ehfx_out = ehfx
         phfx_out = phfx
      end if
      call nwpw_timing_end(33)

      return
      end

*     *************************
*     *                       *
*     *     pspw_energy_HFX0  *
*     *                       *
*     *************************
      subroutine pspw_energy_HFX0(ispin0,psi_r,ehfx_out,phfx_out)
      implicit none
      integer ispin0
      real*8  psi_r(*)
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

      integer q,n,indx1,indx2

      if (((norbs(1)+norbs(2)).ne.0).and.(.not.relaxed)) then

c     **** calculate HFX energy  ****
         if (replicated) then

            call dcopy(nrsize,0.0d0,0,dbl_mb(psi_r_replicated(1)),1)
            do q=1,neqall
               call Dneall_qton(q,n)
               indx1 = (q-1)*n2ft3d + 1
               indx2 = psi_r_replicated(1)+(n-1)*n2ft3d
               call dcopy(n2ft3d,psi_r(indx1),1,dbl_mb(indx2),1)
            end do
            call D1dB_Vector_SumAll(nrsize,dbl_mb(psi_r_replicated(1)))
            call pspw_energy_HFX_sub(ispin0,
     >                               dbl_mb(psi_r_replicated(1)),
     >                               ehfx_out,phfx_out)

         else

            call pspw_energy_HFX_sub(ispin0,psi_r,ehfx_out,phfx_out)

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
*     *     pspw_potential_HFX_orb   *
*     *                              *
*     ********************************
      subroutine pspw_potential_HFX_orb(ms,psi_r,
     >                                  orb_r,Horb_r)
      implicit none
      integer    ms
      real*8     psi_r(*)
      real*8     orb_r(*)
      real*8     Horb_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n1,n2,n3
      integer dn(2),vij(2),tmp1(2),index2
      real*8  scal1,scal2,dv,eh,ph

*     **** external functions ****
      real*8   lattice_omega,coulomb_screened_e
      external lattice_omega,coulomb_screened_e


      call nwpw_timing_start(33)
      if ((norbs(ms).ne.0).and.relaxed) then
        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_dbl,(n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'vij_hfx',vij(2),vij(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(vij(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2

        do j=1,norbs(ms)
           index2 = (int_mb(orbital_list(1,ms)+j-1)-1)*n2ft3d + 1

*          **** generate dnij for Vij  ****
           call D3dB_rr_Mul(1,psi_r(index2),orb_r,dbl_mb(dn(1)))
c           call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
           call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))

*          ***** screened coulomb solver ****
           if (solver_type.eq.1) then
c             call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
c     >                                dbl_mb(dn(1)))
             call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))
             call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
             call Pack_c_pack(0,dbl_mb(dn(1)))

*            **** get Ecoul energy ****
             eh = coulomb_screened_e(dbl_mb(dn(1)))

*            **** generate Vcoul ****
             call coulomb_screened_v(dbl_mb(dn(1)),dbl_mb(vij(1)))
             call Pack_c_unpack(0,dbl_mb(vij(1)))
             !call D3dB_cr_fft3b(1,dbl_mb(vij(1)))
             call D3dB_cr_pfft3b(1,0,dbl_mb(vij(1)))

*          ***** free-space coulomb solver ****
           else
              call coulomb2_v(dbl_mb(dn(1)),dbl_mb(vij(1)))
              call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(vij(1)),eh)
              eh = 0.5d0*eh*dv
           end if

*          **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
           eh = eh*hfx_parameter
c           call D3dB_r_SMul(1,hfx_parameter,
c     >                      dbl_mb(vij(1)),
c     >                      dbl_mb(vij(1)))
           call D3dB_r_SMul1(1,hfx_parameter,dbl_mb(vij(1)))
           if (ispin.eq.1) eh = eh + eh
           ph = 2.0d0*eh

*          **** generate (Vij)*psi_r ***
           call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                        psi_r(index2),
     >                        dbl_mb(tmp1(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))

*          **** add -(Vij)*psi_r to Hpsi_r ***
           call D3dB_rr_Sub(1,Horb_r,
     >                        dbl_mb(tmp1(1)),
     >                        Horb_r)
        end do

        value = value.and.BA_pop_stack(tmp1(2))
        value = value.and.BA_pop_stack(vij(2))
        value = value.and.BA_pop_stack(dn(2))
        if (.not. value) 
     >    call errquit('pspw_potential_HFX_orb:popping stack memory',0,
     &       MA_ERR)
      end if
      call nwpw_timing_end(33)
      return
      end



*     *************************
*     *                       *
*     *     pspw_HFX          *
*     *                       *
*     *************************
      logical function pspw_HFX()
      implicit none

#include "pspw_hfx.fh"

      pspw_HFX = hfx_on
      return
      end

*     *******************************
*     *                             *
*     *     pspw_HFX_virutal_set    *
*     *                             *
*     *******************************
      subroutine pspw_HFX_virtual_set(hfxvset)
      implicit none
      logical hfxvset

#include "pspw_hfx.fh"

      hfx_virtual_on = hfxvset
      return
      end

*     *************************
*     *                       *
*     *     pspw_HFX_virtual  *
*     *                       *
*     *************************
      logical function pspw_HFX_virtual()
      implicit none

#include "pspw_hfx.fh"

      pspw_HFX_virtual = hfx_virtual_on
      return
      end



*     *************************
*     *                       *
*     *   pspw_HFX_localize   *
*     *                       *
*     *************************
      logical function pspw_HFX_localize()
      implicit none

#include "pspw_hfx.fh"

      pspw_HFX_localize = localize_on
      return
      end

*     *************************
*     *                       *
*     *   pspw_HFX_localize2  *
*     *                       *
*     *************************
      logical function pspw_HFX_localize2()
      implicit none

#include "pspw_hfx.fh"

      pspw_HFX_localize2 = localize2_on
      return
      end

*     *************************
*     *                       *
*     *   pspw_HFX_relaxed    *
*     *                       *
*     *************************
      logical function pspw_HFX_relaxed()
      implicit none

#include "bafdecls.fh"
#include "pspw_hfx.fh"

      pspw_hfx_relaxed = relaxed
      return
      end






*     *****************************
*     *                           *
*     *     pspw_energy_euv_HFX   *
*     *                           *
*     *****************************
      subroutine pspw_energy_euv_HFX(ispin0,psi_r,stress)
      implicit none
      integer ispin0
      real*8  psi_r(*)
      real*8 stress(3,3)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

      integer q,n,indx1,indx2

      call nwpw_timing_start(33)
      if (((norbs(1)+norbs(2)).ne.0)) then

*     ***** localize psi_r here, and generate absoverlap and rotation matrix ****

c     **** calculate HFX stress  ****
         if (replicated) then
            call dcopy(nrsize,0.0d0,0,dbl_mb(psi_r_replicated(1)),1)
            do q=1,neqall
               call Dneall_qton(q,n)
               indx1 = (q-1)*n2ft3d + 1
               indx2 = psi_r_replicated(1)+(n-1)*n2ft3d
               call dcopy(n2ft3d,psi_r(indx1),1,dbl_mb(indx2),1)
            end do
            call D1dB_Vector_SumAll(nrsize,dbl_mb(psi_r_replicated(1)))
            call pspw_energy_euv_HFX_sub(ispin0,
     >                               dbl_mb(psi_r_replicated(1)),
     >                               stress)
         else
            call pspw_energy_euv_HFX_sub(ispin0,psi_r,stress)
         end if

      !*** nothing to do ***
      else
         call dcopy(9,0.0d0,0,stress,1)
      end if
      call nwpw_timing_end(33)

      return
      end

c***************** sub/replicated routines *****************************

*     ********************************
*     *                    	     *
*     *     pspw_potential_HFX_sub   *
*     *                              *
*     ********************************
      subroutine pspw_potential_HFX_sub(ispin0,psi_r,Hpsi_r)
      implicit none
      integer    ispin0
      real*8     psi_r(*)
      real*8     Hpsi_r(*)

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
      real*8  scal1,scal2,dv,eh,ph,ss,teh
      integer center(3)
      real*8  rcenter(3)
      integer taskid,icount
      real*8 cpu0,cpu1
      integer ktaskjid,kcompute(2)

*     **** external functions ****
      real*8   lattice_omega,icoulomb_screened_e
      real*8   icoulomb_screened_small_e
      logical  D3dB_rc_pfft3_queue_filled,D3dB_cr_pfft3_queue_filled
      logical  pspw_hfx_localize_closeenough
      real*8   pspw_hfx_localize_switchr
      external lattice_omega,icoulomb_screened_e
      external icoulomb_screened_small_e
      external D3dB_rc_pfft3_queue_filled,D3dB_cr_pfft3_queue_filled
      external pspw_hfx_localize_closeenough
      external pspw_hfx_localize_switchr

      call Parallel2d_taskid_i(taskid)
      icount = 0
      call current_second(cpu0)

      ehfx = 0.0d0
      phfx = 0.0d0
      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then
        !call D3dB_nx(1,n1)
        !call D3dB_ny(1,n2)
        !call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_dbl,(n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'vij_hfx',vij(2),vij(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'tmp2_hfx',tmp2(2),tmp2(1))
        NN = norbs(1)*(norbs(1)+1)/2
        value = value.and.
     >          BA_push_get(mt_int,NN,'kcmpute',kcompute(2),kcompute(1))
        if (.not.value) 
     >     call errquit('pspw_potential_HFX_sub:out of stack memory',
     >                  0,MA_ERR)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(vij(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp2(1)),1)

        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2

        if (localize_on.and.has_smallgrid) then
           call D3dB_nx(3,n1)
           call D3dB_ny(3,n2)
           call D3dB_nz(3,n3)
           scal1 = 1.0d0/dble(n1*n2*n3)
        end if

*       **** compute kcompute ****
        ktaskjid = 0

*       ***** screened coulomb solver ****
        if (solver_type.eq.1) then
        do ms=1,ispin0
           if (norbs(ms).eq.0) go to 898
           call dcopy(norbs(ms),0.0d0,0,dbl_mb(ehfx_orb(1,ms)),1)
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
     >            (.not.D3dB_rc_pfft3_queue_filled())) then

                 !if (mod(k1,npj).eq.taskid_j) then
                 if (int_mb(kcompute(1)+k1-1).eq.taskid_j) then
                 if (pspw_hfx_localize_closeenough(
     >                   int_mb(orbital_list(1,ms)+i1-1),
     >                   int_mb(orbital_list(1,ms)+j1-1))) then

                    index1 =(int_mb(orbital_list(1,ms)+i1-1)-1)*n2ft3d+1
                    index2 =(int_mb(orbital_list(1,ms)+j1-1)-1)*n2ft3d+1

*                   **** generate dnij for Vij  ****
                    call D3dB_rr_Mul(1,psi_r(index2),
     >                                 psi_r(index1),dbl_mb(dn(1)))

                    call D3dB_r_SMul1(1,scal2*scal1,dbl_mb(dn(1)))
                    call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))

c                    if (taskid.eq.0) then
c                       call current_second(cpu1)
c                       icount = icount + 1
c                       write(*,*) "into hfx=",i1,j1,
c     >                            " count,rate,NN=",icount,
c     >                            dble(icount)/(cpu1-cpu0),NN
c                    end if

                    if (localize_on.and.has_smallgrid) then
                       call pspw_hfx_localize_center_ovlp(
     >                          int_mb(orbital_list(1,ms)+i1-1),
     >                          int_mb(orbital_list(1,ms)+j1-1),center)
                       call D3dB_rc_pfft3f_queuein_center(2,center,
     >                                                    dbl_mb(dn(1)))
                    else
                       call D3dB_rc_pfft3f_queuein(0,dbl_mb(dn(1)))
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
 
              if ((     ((D3dB_rc_pfft3_queue_filled()).or.(k1.gt.NN))
     >            .and.(k2.le.NN)).and.
     >            (.not.D3dB_cr_pfft3_queue_filled())) then

                 !if (mod(k2,npj).eq.taskid_j) then
                 if (int_mb(kcompute(1)+k2-1).eq.taskid_j) then
                 if (pspw_hfx_localize_closeenough(
     >                   int_mb(orbital_list(1,ms)+i2-1),
     >                   int_mb(orbital_list(1,ms)+j2-1))) then

                    ss = pspw_hfx_localize_switchr(
     >                   int_mb(orbital_list(1,ms)+i2-1),
     >                   int_mb(orbital_list(1,ms)+j2-1))

                    !eh = coulomb_screened_e(dbl_mb(dn(1)))
                    if (localize_on.and.has_smallgrid) then
                       call D3dB_rc_pfft3f_queueout_center(2,
     >                                             dbl_mb(dn(1)))
                       eh = icoulomb_screened_small_e(dbl_mb(dn(1)))
                       call coulomb_screened_small_v(dbl_mb(dn(1)),
     >                                               dbl_mb(vij(1)))
                    else
                       call D3dB_rc_pfft3f_queueout(0,dbl_mb(dn(1)))
                       eh = icoulomb_screened_e(dbl_mb(dn(1)))
                       call coulomb_screened_v(dbl_mb(dn(1)),
     >                                      dbl_mb(vij(1)))
                    end if

*                   **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
                    eh = eh*hfx_parameter*ss
                    if (ispin0.eq.1) eh = eh + eh
                    ph = 2.0d0*eh
                    ehfx = ehfx - eh
                    phfx = phfx - ph
                    dbl_mb(ehfx_orb(1,ms)+i2-1) 
     >               = dbl_mb(ehfx_orb(1,ms)+i2-1) - eh
                    if (i2.ne.j2) then
                       ehfx = ehfx - eh
                       phfx = phfx - ph
                       dbl_mb(ehfx_orb(1,ms)+i2-1) 
     >                  = dbl_mb(ehfx_orb(1,ms)+i2-1) - eh
                    end if

                    if (localize_on.and.has_smallgrid) then
                       call pspw_hfx_localize_center_ovlp(
     >                          int_mb(orbital_list(1,ms)+i2-1),
     >                          int_mb(orbital_list(1,ms)+j2-1),center)
                       call D3dB_cr_pfft3b_queuein_center(2,center,
     >                                                   dbl_mb(vij(1)))
                    else
                       call D3dB_cr_pfft3b_queuein(0,dbl_mb(vij(1)))
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

              if ((D3dB_cr_pfft3_queue_filled()).or.(k2.gt.NN)) then

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
                       call D3dB_cr_pfft3b_queueout_center(2,
     >                                              dbl_mb(vij(1)))
                    else
                       call D3dB_cr_pfft3b_queueout(0,dbl_mb(vij(1)))
                    end if


*                   **** apply hfx_parameter ****
                    call D3dB_r_SMul1(1,hfx_parameter*ss,
     >                                dbl_mb(vij(1)))

*                   **** generate (Vij)*psi_r ***
                    call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                                 psi_r(index2),
     >                                 dbl_mb(tmp1(1)))
                    call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))

*                   **** add -(Vij)*psi_r to Hpsi_r ***
c                    call D3dB_rr_Sub(1,Hpsi_r(index1),
c     >                                 dbl_mb(tmp1(1)),
c     >                                 Hpsi_r(index1))
                    call D3dB_rr_Sub2(1,dbl_mb(tmp1(1)),Hpsi_r(index1))

                    !**** include off diagonal terms ****
                    if (i3.ne.j3) then
*                      **** generate (Vij)*psi_r ***
                       call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                                 psi_r(index1),
     >                                 dbl_mb(tmp2(1)))
                       call D3dB_r_Zero_Ends(1,dbl_mb(tmp2(1)))

*                      **** add -(Vij)*psi_r to Hpsi_r ***
c                       call D3dB_rr_Sub(1,Hpsi_r(index2),
c     >                                 dbl_mb(tmp2(1)),
c     >                                 Hpsi_r(index2))
                       call D3dB_rr_Sub2(1,dbl_mb(tmp2(1)),
     >                                     Hpsi_r(index2))
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
           !call D1dB_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
           call Parallel_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))

  898      continue
       end do !**** ms *****

*       ***** free-space coulomb solver ****
        else
        k1 = 1
        do ms=1,ispin0
        do i=1,norbs(ms)
         dbl_mb(ehfx_orb(1,ms)+i-1) = 0.0d0
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
              call D3dB_rr_Mul(1,psi_r(index2),psi_r(index1),
     >                         dbl_mb(dn(1)))
c              call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
              call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
              call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))
   
              call coulomb2_v(dbl_mb(dn(1)),dbl_mb(vij(1)))
              call D3dB_rr_idot(1,dbl_mb(dn(1)),dbl_mb(vij(1)),eh)
              eh = 0.5d0*eh*dv

*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
              eh = eh*hfx_parameter*ss
c              call D3dB_r_SMul(1,hfx_parameter,
c     >                         dbl_mb(vij(1)),
c     >                         dbl_mb(vij(1)))
              call D3dB_r_SMul1(1,hfx_parameter*ss,dbl_mb(vij(1)))
              if (ispin0.eq.1) eh = eh + eh
              ph = 2.0d0*eh


              ehfx = ehfx - eh
              phfx = phfx - ph
              dbl_mb(ehfx_orb(1,ms)+i-1) =dbl_mb(ehfx_orb(1,ms)+i-1)-eh

*             **** generate (Vij)*psi_r ***
              call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                           psi_r(index2),
     >                           dbl_mb(tmp1(1)))
              call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))


*             **** add -(Vij)*psi_r to Hpsi_r ***
c              call D3dB_rr_Sub(1,Hpsi_r(index1),
c     >                           dbl_mb(tmp1(1)),
c     >                           Hpsi_r(index1))
              call D3dB_rr_Sub2(1,dbl_mb(tmp1(1)),Hpsi_r(index1))

              !**** include off diagonal terms ****
              if (i.ne.j) then
                 ehfx = ehfx - eh
                 phfx = phfx - ph
                 dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)
     >                                      - eh
*                **** generate (Vij)*psi_r ***
                 call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                           psi_r(index1),
     >                           dbl_mb(tmp2(1)))
                 call D3dB_r_Zero_Ends(1,dbl_mb(tmp2(1)))

*                **** add -(Vij)*psi_r to Hpsi_r ***
c                 call D3dB_rr_Sub(1,Hpsi_r(index2),
c     >                           dbl_mb(tmp2(1)),
c     >                           Hpsi_r(index2))
                 call D3dB_rr_Sub2(1,dbl_mb(tmp2(1)),
     >                               Hpsi_r(index2))
              end if
           end if
           end if
           k1 = k1 + 1
         end do
        end do
        !call D1dB_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
        call Parallel_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
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
*     *     pspw_energy_HFX_sub   *
*     *                           *
*     *****************************
      subroutine pspw_energy_HFX_sub(ispin0,psi_r,ehfx_out,phfx_out)
      implicit none
      integer ispin0
      real*8  psi_r(*)
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n1,n2,n3,ms,k1
      integer dn(2),tmp1(2),index1,index2
      real*8  scal1,scal2,dv,eh,ph,ss

*     **** external functions ****
      real*8   lattice_omega,coulomb_screened_e
      external lattice_omega,coulomb_screened_e
      logical  pspw_hfx_localize_closeenough
      external pspw_hfx_localize_closeenough
      real*8   pspw_hfx_localize_switchr
      external pspw_hfx_localize_switchr

      
      if (((norbs(1)+norbs(2)).ne.0).and.(.not.relaxed)) then
        ehfx = 0.0d0
        phfx = 0.0d0

        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_dbl,(2*n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2

        k1 = 1
        do ms=1,ispin
        do i=1,norbs(ms)
         dbl_mb(ehfx_orb(1,ms)+i-1) = 0.0d0
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
              call D3dB_rr_Mul(1,psi_r(index1),psi_r(index2),
     >                         dbl_mb(dn(1)))
c              call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
              call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
              call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))

*             ***** screened coulomb solver ****
              if (solver_type.eq.1) then

*               **** generate dng ****
c                call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
c     >                                   dbl_mb(dn(1)))
                call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))
                call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
                call Pack_c_pack(0,dbl_mb(dn(1)))

*               **** get Ecoul energy ****
                eh = coulomb_screened_e(dbl_mb(dn(1)))
       
*             ***** free-space coulomb solver ****
              else
                 call coulomb2_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))
                 call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(tmp1(1)),eh)
                 eh = 0.5d0*eh*dv
              end if

*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
              eh = eh*hfx_parameter*ss
              if (ispin0.eq.1) eh = eh + eh
              ph = 2.0d0*eh

              ehfx = ehfx - eh
              phfx = phfx - ph
              dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)-eh

              !**** include off diagonal terms ****
              if (i.ne.j) then
                 ehfx = ehfx - eh
                 phfx = phfx - ph
                 dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)
     >                                      - eh
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


*     *********************************
*     *                               *
*     *     pspw_energy_euv_HFX_sub   *
*     *                               *
*     *********************************
      subroutine pspw_energy_euv_HFX_sub(ispin0,psi_r,stress)
      implicit none
      integer ispin0
      real*8  psi_r(*)
      real*8 stress(3,3)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n1,n2,n3,ms,u,v,k1
      integer dn(2),index1,index2
      real*8  scal1,scal2,dv,eh,ph,ss
      real*8  tstress(3,3)

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      logical  pspw_hfx_localize_closeenough
      external pspw_hfx_localize_closeenough
      real*8   pspw_hfx_localize_switchr

      
      call dcopy(9,0.0d0,0,stress,1)
      if (((norbs(1)+norbs(2)).ne.0)) then
        call coulomb_screened_euv_init(flag,rcut,pp)

        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_dbl,(2*n2ft3d),'dn_hfx',dn(2),dn(1))
        if (.not. value) call errquit('out of stack memory',0,MA_ERR)
        call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2

        k1 = 1
        do ms=1,ispin
        do i=1,norbs(ms)
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
              call D3dB_rr_Mul(1,psi_r(index1),psi_r(index2),
     >                         dbl_mb(dn(1)))
c              call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
              call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
              call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))

*             ***** screened coulomb solver ****
              if (solver_type.eq.1) then

*               **** generate dng ****
c                call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
c     >                                   dbl_mb(dn(1)))
                call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))
                call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
                call Pack_c_pack(0,dbl_mb(dn(1)))

*               **** get Ecoul energy ****
                call coulomb_screened_euv(dbl_mb(dn(1)),tstress)
                if (ispin.eq.1) call dscal(9,2.0d0,tstress,1)

*               **** apply the hfx_parameter ****
                call dscal(9,hfx_parameter*ss,tstress,1)

*             ***** free-space coulomb solver ****
              else
                write(*,*) "ERROR free-space coulomb solver called"
                call errquit('error: not periodic boundary conditions',
     >                       0,0)
              end if

              do v=1,3
              do u=1,3
                 stress(u,v) = stress(u,v) - tstress(u,v)
              end do
              end do

              !**** include off diagonal terms ****
              if (i.ne.j) then
                do v=1,3
                do u=1,3
                   stress(u,v) = stress(u,v) - tstress(u,v)
                end do
                end do
              end if

            end if
            end if
            k1 = k1 + 1

         end do
        end do
        end do

        call coulomb_screened_euv_end()
        value = BA_pop_stack(dn(2))
        if (.not. value) 
     >     call errquit('pspw_energy_euv_HFX_sub:popping stack memory',
     >                  0,MA_ERR)

        call D1dB_Vector_SumAll(9,stress)
      end if

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
      subroutine pspw_potential_HFX_orb_sub(ms,psi_r,
     >                                      orb_r,Horb_r)
      implicit none
      integer    ms
      real*8     psi_r(*)
      real*8     orb_r(*)
      real*8     Horb_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer j,n1,n2,n3,q2,p2
      integer dn(2),vij(2),tmp1(2),tmp2(2),index2
      real*8  scal1,scal2,dv,eh,ph

*     **** external functions ****
      real*8   lattice_omega,coulomb_screened_e
      external lattice_omega,coulomb_screened_e

      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then
        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_dbl,(n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'vij_hfx',vij(2),vij(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'tmp2_hfx',tmp2(2),tmp2(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(vij(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp2(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2

        do j=1,norbs(ms)
           call Dneall_ntoqp(int_mb(orbital_list(1,ms)+j-1),q2,p2)
           index2 = (q2-1)*n2ft3d + 1

           if (p2.eq.taskid_j) then
*             **** generate dnij for Vij  ****
              call D3dB_rr_Mul(1,psi_r(index2),orb_r,dbl_mb(dn(1)))
              call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
              call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))
   
*             ***** screened coulomb solver ****
              if (solver_type.eq.1) then
                call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))
                call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
                call Pack_c_pack(0,dbl_mb(dn(1)))

*               **** get Ecoul energy ****
                eh = coulomb_screened_e(dbl_mb(dn(1)))

*               **** generate Vcoul ****
                call coulomb_screened_v(dbl_mb(dn(1)),dbl_mb(vij(1)))
                call Pack_c_unpack(0,dbl_mb(vij(1)))
                call D3dB_cr_pfft3b(1,0,dbl_mb(vij(1)))

*             ***** free-space coulomb solver ****
              else
                 call coulomb2_v(dbl_mb(dn(1)),dbl_mb(vij(1)))
                 call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(vij(1)),eh)
                 eh = 0.5d0*eh*dv
              end if

*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
              eh = eh*hfx_parameter
              call D3dB_r_SMul1(1,hfx_parameter,dbl_mb(vij(1)))
              if (ispin.eq.1) eh = eh + eh
              ph = 2.0d0*eh


*             **** generate (Vij)*psi_r ***
              call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                           psi_r(index2),
     >                           dbl_mb(tmp1(1)))
              call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))

*             **** add -(Vij)*psi_r to Hpsi_r ***
              call D3dB_rr_Sub2(1,dbl_mb(tmp1(1)),Horb_r)
           end if
        end do
        value =           BA_pop_stack(tmp2(2))
        value = value.and.BA_pop_stack(tmp1(2))
        value = value.and.BA_pop_stack(vij(2))
        value = value.and.BA_pop_stack(dn(2))
        if (.not. value) 
     >   call errquit('pspw_potential_HFX_orb_sub:popping stack memory',
     >                  0,MA_ERR)

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

*   Exit - Hpsi_r: wavefunction gradients in realspace.
*          ehfx: running sum of exchange energy.

      subroutine pspw_potential_HFX_sub2(solver_type,
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
      real*8  psi_r(n2ft3d,*)
      real*8  Hpsi_r(n2ft3d,*)
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
      real*8  scal1,scal2,dv,eh,ph

*     **** external functions ****
      real*8   lattice_omega,icoulomb_screened_e
      external lattice_omega,icoulomb_screened_e
      logical  D3dB_rc_pfft3_queue_filled,D3dB_cr_pfft3_queue_filled
      external D3dB_rc_pfft3_queue_filled,D3dB_cr_pfft3_queue_filled

      call Parallel2d_taskid_j(taskid_j)

      special = ((istart.eq.jstart).and.(iend.eq.jend))

      call D3dB_nx(1,n1)
      call D3dB_ny(1,n2)
      call D3dB_nz(1,n3)
      value = BA_push_get(mt_dbl,(n2ft3d),'dn_hfx',dn(2),dn(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'vij_hfx',vij(2),vij(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('pspw_potential_HFX_sub: out of stack',0,MA_ERR)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(vij(1)),1)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)

      scal1 = 1.0d0/dble(n1*n2*n3)
      scal2 = 1.0d0/lattice_omega()
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
                 call D3dB_rr_Mul(1,psi_r(1,j1),
     >                              psi_r(1,i1),dbl_mb(dn(1)))
                 call D3dB_r_SMul1(1,scal2*scal1,dbl_mb(dn(1)))
                 call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))
                 call D3dB_rc_pfft3f_queuein(0,dbl_mb(dn(1)))

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
           if (     ((D3dB_rc_pfft3_queue_filled()).or.(k1.gt.NN))
     >         .and.(k2.le.NN)) then

              if (mod(k2,imodn).eq.imodtask) then
                 call D3dB_rc_pfft3f_queueout(0,dbl_mb(dn(1)))

*                **** generate Vcoul ****
                 eh = icoulomb_screened_e(dbl_mb(dn(1)))
                 call coulomb_screened_v(dbl_mb(dn(1)),
     >                                   dbl_mb(vij(1)))


*                **** calculcate ph ****
                 ehfx = ehfx - eh

*                **** include transpose ***
                 if ((i2.ne.j2).or.(.not.special)) then
                    ehfx = ehfx - eh
                 end if

                 call D3dB_cr_pfft3b_queuein(0,dbl_mb(vij(1)))
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
           if ((D3dB_cr_pfft3_queue_filled()).or.(k2.gt.NN)) then

              if (mod(k3,imodn).eq.imodtask) then
                 call D3dB_cr_pfft3b_queueout(0,dbl_mb(vij(1)))

*                **** generate (Vij)*psi_r ***
                 call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                              psi_r(1,j3),
     >                              dbl_mb(tmp1(1)))
                 call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))

*                **** add -(Vij)*psi_r to Hpsi_r ***
                 call D3dB_rr_Sub2(1,dbl_mb(tmp1(1)),Hpsi_r(1,i3))

                 !**** include transpose ****
                 if ((i3.ne.j3).or.(.not.special)) then

*                   **** generate (Vij)*psi_r ***
                    call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                              psi_r(1,i3),
     >                              dbl_mb(tmp1(1)))
                    call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))

*                   **** add -(Vij)*psi_r to Hpsi_r ***
                    call D3dB_rr_Sub2(1,dbl_mb(tmp1(1)),
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
               call D3dB_rr_Mul(1,psi_r(1,j1),psi_r(1,i1),
     >                          dbl_mb(dn(1)))
               call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
               call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))
   
               call coulomb2_v(dbl_mb(dn(1)),dbl_mb(vij(1)))
               call D3dB_rr_idot(1,dbl_mb(dn(1)),dbl_mb(vij(1)),eh)
               eh = 0.5d0*eh*dv

*              **** calculcate ph ****
               ehfx = ehfx - eh

*              **** generate (Vij)*psi_r ***
               call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                            psi_r(1,j1),
     >                            dbl_mb(tmp1(1)))
               call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))

*              **** add -(Vij)*psi_r to Hpsi_r ***
               call D3dB_rr_Sub2(1,dbl_mb(tmp1(1)),Hpsi_r(1,i1))
 
               !**** include transpose terms ****
               if ((i1.ne.j1).or.(.not.special)) then
                  ehfx = ehfx - eh
 
*                 **** generate (Vij)*psi_r ***
                  call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                            psi_r(1,i1),
     >                            dbl_mb(tmp1(1)))
                  call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))
 
*                 **** add -(Vij)*psi_r to Hpsi_r ***
                  call D3dB_rr_Sub2(1,dbl_mb(tmp1(1)),
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
     >   'pspw_potential_HFX_sub2:popping stack memory',0,MA_ERR)

      return
      end



*     *******************************************
*     *                    	                *
*     *     pspw_potential_HFX_orb_replicated   *
*     *                                         *
*     *******************************************
      subroutine pspw_potential_HFX_orb_replicated(ms,psi_r,
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
            call dcopy(n2ft3d,0.0d0,0,dbl_mb(Hpsi_r_replicated(1)),1)
            call pspw_potential_HFX_orb_sub(ms,psi_r,orb_r,
     >                                  dbl_mb(Hpsi_r_replicated(1)))
            call D1dB_Vector_SumAll(n2ft3d,dbl_mb(Hpsi_r_replicated(1)))
            call daxpy_omp(n2ft3d,1.0d0,dbl_mb(Hpsi_r_replicated(1)),1,
     >                 Horb_r,1)

         else
            call pspw_potential_HFX_orb_sub(ms,psi_r,orb_r,Horb_r)
         end if

      end if
      call nwpw_timing_end(33)
      return
      end






*    ***** routines below used for pspw_et calculations ****


*     *************************
*     *                       *
*     *     pspw_energy_HFX2  *
*     *                       *
*     *************************
      subroutine pspw_energy_HFX2(ispin0,psi1_r,psi2_r,
     >                            ehfx_out,phfx_out)
      implicit none
      integer ispin0
      real*8  psi1_r(*),psi2_r(*)
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
            call dcopy(nrsize,0.0d0,0,dbl_mb(psi_r_replicated(1)),1)
            call dcopy(nrsize,0.0d0,0,dbl_mb(Hpsi_r_replicated(1)),1)
            do q=1,neqall
               call Dneall_qton(q,n)
               indx1 = (q-1)*n2ft3d + 1
               indx2 = psi_r_replicated(1)+(n-1)*n2ft3d
               call dcopy(n2ft3d,psi1_r(indx1),1,dbl_mb(indx2),1)

               indx2 = Hpsi_r_replicated(1)+(n-1)*n2ft3d
               call dcopy(n2ft3d,psi2_r(indx1),1,dbl_mb(indx2),1)
            end do
            call D1dB_Vector_SumAll(nrsize,dbl_mb(psi_r_replicated(1)))
            call D1dB_Vector_SumAll(nrsize,dbl_mb(Hpsi_r_replicated(1)))
            call pspw_energy_HFX2_sub(ispin0,
     >                               dbl_mb(psi_r_replicated(1)),
     >                               dbl_mb(Hpsi_r_replicated(1)),
     >                               ehfx_out,phfx_out)

         else
            call pspw_energy_HFX2_sub(ispin0,psi1_r,psi2_r,
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
*     *     pspw_energy_HFX2_sub  *
*     *                           *
*     *****************************
      subroutine pspw_energy_HFX2_sub(ispin0,psi1_r,psi2_r,
     >                               ehfx_out,phfx_out)
      implicit none
      integer ispin0
      real*8  psi1_r(*)
      real*8  psi2_r(*)
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n1,n2,n3,ms,k1
      integer dn(2),tmp1(2),index1,index2
      real*8  scal1,scal2,dv,eh,ph

*     **** external functions ****
      real*8   lattice_omega,coulomb_screened_e
      external lattice_omega,coulomb_screened_e

      integer taskid

      call Parallel_taskid(taskid)

      
      if ((norbs(1)+norbs(2)).ne.0) then
        ehfx = 0.0d0
        phfx = 0.0d0

        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        !call D3dB_n2ft3d(1,n2ft3d)
        value = BA_push_get(mt_dbl,(2*n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)
        call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2

        k1 = 1
        do ms=1,ispin
        do i=1,norbs(ms)
         dbl_mb(ehfx_orb(1,ms)+i-1) = 0.0d0
         do j=1,i

            if (mod(k1,npj).eq.taskid_j) then
              index1 = (int_mb(orbital_list(1,ms)+i-1)-1)*n2ft3d + 1
              index2 = (int_mb(orbital_list(1,ms)+j-1)-1)*n2ft3d + 1

*             **** generate dnij ****
              call D3dB_rr_Mul(1,psi1_r(index1),psi2_r(index2),
     >                         dbl_mb(dn(1)))
              call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
              call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))

*             **** generate  dnji****
              call D3dB_rr_Mul(1,psi1_r(index2),psi2_r(index1),
     >                         dbl_mb(dn(1)+n2ft3d))
              call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)+n2ft3d))
              call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)+n2ft3d))

*             ***** screened coulomb solver ****
              if (solver_type.eq.1) then

*               **** generate dng ****
                call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))
                call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
                call Pack_c_pack(0,dbl_mb(dn(1)))

                call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)+n2ft3d))
                call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)+n2ft3d))
                call Pack_c_pack(0,dbl_mb(dn(1)+n2ft3d))

*               **** get Ecoul energy ****
                call coulomb_screened_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))
                call Pack_cc_dot(0,dbl_mb(dn(1)+n2ft3d),
     >                             dbl_mb(tmp1(1)),eh)
                eh = 0.5d0*eh*lattice_omega()
       
*             ***** free-space coulomb solver ****
              else
                 call coulomb2_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))
                 call D3dB_rr_dot(1,dbl_mb(dn(1)+n2ft3d),
     >                            dbl_mb(tmp1(1)),eh)
                 eh = 0.5d0*eh*dv

                 !if (taskid.eq.0) write(*,*) " - ms,i,j,ehfx=",ms,i,j,eh

              end if

*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
              eh = eh*hfx_parameter
              if (ispin0.eq.1) eh = eh + eh
              ph = 2.0d0*eh

              ehfx = ehfx - eh
              phfx = phfx - ph
              dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)-eh

              !**** include off diagonal terms ****
              if (i.ne.j) then
                 ehfx = ehfx - eh
                 phfx = phfx - ph
                 dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)
     >                                      - eh
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
     >     call errquit('pspw_energy_HFX2_sub:popping stack memory',0,
     &       MA_ERR)

        call D1dB_SumAll(ehfx)
        call D1dB_SumAll(phfx)
      end if
      ehfx_out = ehfx
      phfx_out = phfx

      return
      end


*     *****************************
*     *                           *
*     *   pspw_Lin_proj_HFX_set   *
*     *                           *
*     *****************************

      subroutine pspw_Lin_proj_HFX_set(psi_r,Hpsi_r,Lin_prj)
      implicit none
      real*8  psi_r(*)
      real*8  Hpsi_r(*)
      real*8  Lin_prj(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      integer hml(2),neq(2),nsize,n1,n2,n3
      real*8 scal1

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      call D3dB_nx(1,n1)
      call D3dB_ny(1,n2)
      call D3dB_nz(1,n3)
      scal1 = 1.0d0/dble(n1*n2*n3)

      if (.not.Dneall_m_push_get(0,hml))
     >   call errquit("pspw_Lin_Proj_HFX_set:out of stack",0,MA_ERR)

      call Dneall_neq(neq)
      nsize = (neq(1)+neq(2))*n2ft3d

      call Parallel_shared_vector_copy(.false.,nsize,Hpsi_r,Lin_prj)
c      call Dneall_ggm_sym_Multiply(0,psi_r,Hpsi_r,n2ft3d,dbl_mb(hml(1)))
c      call Dneall_m_scal(0,scal1,dbl_mb(hml(1)))
      call Dneall_ggm_Multiply(0,psi_r,Hpsi_r,-scal1,n2ft3d,
     >                         dbl_mb(hml(1)),0.0d0)

      call Dneall_m_cholesky(0,dbl_mb(hml(1)))
      call Dneall_mg_forwardsolve(0,dbl_mb(hml(1)),n2ft3d,Lin_prj)

      if (.not.Dneall_m_pop_stack(hml))
     >   call errquit("pspw_Lin_Proj_HFX_set:pop stack",0,MA_ERR)
      return
      end 


*     *****************************
*     *                           *
*     *   pspw_Lin_proj_HFX       *
*     *                           *
*     *****************************

      subroutine pspw_Lin_proj_HFX(Lin_prj,psi_r,Hpsi_r)
      implicit none
      real*8  psi_r(*)
      real*8  Lin_prj(*)
      real*8  Hpsi_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      integer hml(2),neq(2),n1,n2,n3
      real*8 scal1

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack
      real*8   Dneall_m_sqr_trace
      external Dneall_m_sqr_trace

      call D3dB_nx(1,n1)
      call D3dB_ny(1,n2)
      call D3dB_nz(1,n3)
      scal1 = 1.0d0/dble(n1*n2*n3)

      if (.not.Dneall_m_push_get(0,hml))
     >   call errquit("pspw_Lin_Proj_HFX:out of stack",0,MA_ERR)

      call Dneall_ggm_Multiply(0,Lin_prj,psi_r,scal1,n2ft3d,
     >                         dbl_mb(hml(1)),0.0d0)
      call Dneall_gmg_Multiply(0,Lin_prj,n2ft3d,
     >                         dbl_mb(hml(1)),-1.0d0,
     >                         Hpsi_r,0.0d0)

      phfx = Dneall_m_sqr_trace(0,dbl_mb(hml(1)))
      phfx = -phfx
      if (ispin.eq.1) phfx = phfx + phfx
      ehfx = 0.5d0*phfx

      if (.not.Dneall_m_pop_stack(hml))
     >   call errquit("pspw_Lin_Proj_HFX:pop stack",0,MA_ERR)
      return
      end

*     *****************************
*     *                           *
*     *   pspw_Lin_HFX_reset      *
*     *                           *
*     *****************************
      subroutine pspw_Lin_HFX_reset()
      implicit none
#include "pspw_hfx.fh"

!$OMP MASTER
      lin_run = .false.
!$OMP END MASTER
      return
      end

*     *****************************
*     *                           *
*     *   pspw_Lin_HFX_off        *
*     *                           *
*     *****************************
      subroutine pspw_Lin_HFX_off(bvalue)
      implicit none
      logical bvalue

#include "pspw_hfx.fh"

      !lin_off = .true.
!$OMP MASTER
      lin_off = bvalue
!$OMP END MASTER
      return
      end

*     *****************************
*     *                           *
*     *   pspw_Lin_HFX_on         *
*     *                           *
*     *****************************
      logical function pspw_Lin_HFX_on()
      implicit none
#include "pspw_hfx.fh"
      pspw_Lin_HFX_on = lin_on
      return
      end


