*
* $Id: band_hfx.F,v 1.24 2009-03-30 14:16:06 bylaska Exp $
*

*     *************************
*     *                       *
*     *     band_init_HFX     *
*     *                       *
*     *************************
      subroutine band_init_HFX(rtdb,nbrill0,ispin0,ne)
      implicit none
      integer rtdb
      integer nbrill0,ispin0
      integer ne(2)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "band_hfx.fh"

*     **** local variables ****
      logical value
      integer ma_type
      integer n1,n2,n3,mapping,ms,neq(2)

*     **** external functions ****
      integer  control_version,control_mapping,Butter_size
      external control_version,control_mapping,Butter_size

      nbrillioun = nbrill0
      ispin = ispin0
      norbs(1) = 0
      norbs(2) = 0
      ehfx = 0.0d0
      hfx_on = .false.
      call C3dB_n2ft3d(1,n2ft3d)

      if (.not.rtdb_get(rtdb,'band:HFX',mt_log,1,hfx_on))
     >   hfx_on = .false.


*     **** get the number of HFX orbitals ****
      if (hfx_on) then
         norbs(1) = ne(1)
         norbs(2) = ne(2)

         if (.not. rtdb_get(rtdb,
     >                      'band:HFX_screening_radius',
     >                      mt_dbl,1,rcut)) 
     >       rcut = 8.0d0

         if (.not. rtdb_get(rtdb,
     >                      'band:HFX_screening_power',
     >                      mt_dbl,1,pp)) 
     >       pp = 8.0d0

         if (.not. rtdb_get(rtdb,
     >                      'band:HFX_screening_type',
     >                      mt_int,1,flag)) 
     >       flag = 0

         if (.not. rtdb_get(rtdb,
     >                      'band:HFX_relax',
     >                      mt_log,1,relaxed)) 
     >       relaxed = .true.

         if (.not. rtdb_get(rtdb,
     >                      'band:HFX_solver_type',
     >                      mt_int,1,solver_type)) then

            solver_type = 1
         end if

         if (.not. rtdb_get(rtdb,
     >                      'band:HFX_parameter',
     >                       mt_dbl,1,HFX_parameter))
     >       HFX_parameter = 1.0d0

         if (.not. rtdb_get(rtdb,
     >                      'band:HFX_print_orbital_contribution',
     >                       mt_log,1,orb_contribution))
     >       orb_contribution = .false.
 

*        **** initialize coulomb_screened ****
c         call coulomb_screened_init(flag,rcut,pp)

*      
        
         

*        **** initialize orb_contribution ****
         do ms=1,ispin
           value = MA_alloc_get(mt_dbl,norbs(ms),
     >                'ehfx_orb',ehfx_orb(2,ms),ehfx_orb(1,ms))
           if (.not. value)   
     >       call errquit('band_init_HFX: out of heap memory',1, MA_ERR)
         end do

      end if


c     **** define extra psi and Hpsi  ****
   
      return
      end


*     *************************
*     *                       *
*     *     band_end_HFX      *
*     *                       *
*     *************************
      subroutine band_end_HFX()
      implicit none

#include "mafdecls.fh"
#include "band_hfx.fh"
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
     >            ms,i,
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
        value = .true.
        do ms=1,ispin
          value = value.and.MA_free_heap(ehfx_orb(2,ms))
        end do

*        **** end coulomb_screened ****
c         call coulomb_screened_end()

      end if

      return
      end

*     *************************
*     *                       *
*     *     band_print_HFX    *
*     *                       *
*     *************************
      subroutine band_print_HFX(unit)
      implicit none
      integer unit

#include "mafdecls.fh"
#include "band_hfx.fh"

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
        write(unit,1006)
        if (rcut.ge.0.0d0) write(unit,1008) rcut
        if (rcut.ge.0.0d0) write(unit,1009) pp
        if (rcut.ge.0.0d0) write(unit,1011) flag
        if ((rcut.ge.0.0d0).and.(flag.eq.2)) 
     >     write(unit,1012) control_attenuation()
        if (hfx_parameter.ne.1.0d0) write(unit,1010) hfx_parameter
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
 1008 FORMAT(6x,"- HFX screening radius(band:HFX_screening_radius):",
     >       E10.3)
 1009 FORMAT(6x,"- HFX screening power (band:HFX_screening_power) :",
     >       E10.3)
 1010 FORMAT(6x,"- HFX scaling parameter (band:HFX_parameter)     :",
     >       E10.3)
 1011 FORMAT(6x,"- HFX screening type (band:HFX_screening_type)   :",
     >       I2)
 1012 FORMAT(6x,"- attenuation parameter (nwpw:attenuation)       :",
     >       E10.3)
      end



*     ****************************
*     *                    	 *
*     *     band_potential_HFX   *
*     *                          *
*     ****************************
      subroutine band_potential_HFX(ispin0,psi_r_tag,Hpsi_r_tag)
      implicit none
      integer    ispin0
      integer    psi_r_tag
      integer    Hpsi_r_tag

#include "mafdecls.fh"
#include "band_hfx.fh"
#include "errquit.fh"

      integer istart,iend,jstart,jend,imodn,imodtask
      integer ms,l,q,n,indx1,indx2,Levels,neq(2)
      integer requests(5),reqcnt

      integer  Butter_Levels,Dneall_na_ptr
      external Butter_Levels,Dneall_na_ptr

      call nwpw_timing_start(33)
      ehfx = 0.0d0
      phfx = 0.0d0
      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then

c      call band_potential_HFX_sub(ispin0,psi_r_tag,Hpsi_r_tag)

      end if
      call nwpw_timing_end(33)

      return
      end


*     *************************
*     *                       *
*     *     band_energy_HFX   *
*     *                       *
*     *************************
      subroutine band_energy_HFX(ispin0,psi_r_tag,ehfx_out,phfx_out)
      implicit none
      integer ispin0
      integer psi_r_tag
      real*8 ehfx_out
      real*8 phfx_out

#include "mafdecls.fh"
#include "band_hfx.fh"
#include "errquit.fh"

      integer q,n,indx1,indx2

      call nwpw_timing_start(33)
c     **** calculate HFX energy  ****
      if (((norbs(1)+norbs(2)).ne.0).and.(.not.relaxed)) then
c        call band_energy_HFX_sub(ispin0,psi_r_tag,ehfx_out,phfx_out)

c     **** nothing to do ****
      else
         ehfx_out = ehfx
         phfx_out = phfx
      end if
      call nwpw_timing_end(33)

      return
      end


*     *************************
*     *                       *
*     *     band_HFX          *
*     *                       *
*     *************************
      logical function band_HFX()
      implicit none

#include "band_hfx.fh"

      band_HFX= hfx_on
      return
      end

*     *************************
*     *                       *
*     *   band_HFX_relaxed    *
*     *                       *
*     *************************
      logical function band_HFX_relaxed()
      implicit none

#include "mafdecls.fh"
#include "band_hfx.fh"

      band_hfx_relaxed = relaxed
      return
      end







c***************** sub/replicated routines *****************************

*     ********************************
*     *                    	     *
*     *     band_potential_HFX_sub   *
*     *                              *
*     ********************************
      subroutine band_potential_HFX_sub(ispin0,psi_r_tag,Hpsi_r_tag)
      implicit none
      integer    ispin0
      integer    psi_r_tag
      integer    Hpsi_r_tag

#include "mafdecls.fh"
#include "band_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value,done
      integer i,j,n1,n2,n3,ms
      integer dn(2),vij(2),tmp1(2),tmp2(2),index1,index2
      integer i1,j1,k1,NN
      integer i2,j2,k2
      integer i3,j3,k3
      real*8  scal1,scal2,dv,eh,ph

*     **** external functions ****
      real*8   lattice_omega,icoulomb_screened_e
      logical   D3dB_rc_pfft3_queue_filled,D3dB_cr_pfft3_queue_filled
      external lattice_omega,icoulomb_screened_e
      external  D3dB_rc_pfft3_queue_filled,D3dB_cr_pfft3_queue_filled


      ehfx = 0.0d0
      phfx = 0.0d0
      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then
c        call D3dB_nx(1,n1)
c        call D3dB_ny(1,n2)
c        call D3dB_nz(1,n3)
c        !call D3dB_n2ft3d(1,n2ft3d)
c        value = MA_push_get(mt_dbl,(n2ft3d),'dn_hfx',dn(2),dn(1))
c        value = value.and.
c     >          MA_push_get(mt_dbl,(n2ft3d),'vij_hfx',vij(2),vij(1))
c        value = value.and.
c     >          MA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
c        value = value.and.
c     >          MA_push_get(mt_dbl,(n2ft3d),'tmp2_hfx',tmp2(2),tmp2(1))
c        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
c
c        scal1 = 1.0d0/dble(n1*n2*n3)
c        scal2 = 1.0d0/lattice_omega()
c        dv = scal1/scal2
c
c*       ***** screened coulomb solver ****
c        do ms=1,ispin0
c           call dcopy(norbs(ms),0.0d0,0,dbl_mb(ehfx_orb(1,ms)),1)
c           NN = norbs(ms)*(norbs(ms)+1)/2
c           i1 = 1
c           j1 = 1
c           k1 = 1
c           i2 = 1
c           j2 = 1
c           k2 = 1
c           i3 = 1
c           j3 = 1
c           k3 = 1
c           done = .false.
c           do while (.not.done)
c
c              if (k1.le.NN) then
c
c                 if (mod(k1,npj).eq.taskid_j) then
c                    index1 =(int_mb(orbital_list(1,ms)+i1-1)-1)*n2ft3d+1
c                    index2 =(int_mb(orbital_list(1,ms)+j1-1)-1)*n2ft3d+1
c
c*                   **** generate dnij for Vij  ****
c                    call D3dB_rr_Mul(1,psi_r(index2),
c     >                                 psi_r(index1),dbl_mb(dn(1)))
cc                    call D3dB_r_SMul(1,scal2*scal1,dbl_mb(dn(1)),
cc     >                                             dbl_mb(dn(1)))
c                    call D3dB_r_SMul1(1,scal2*scal1,dbl_mb(dn(1)))
c                    call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))
c
c                    call D3dB_rc_pfft3f_queuein(0,dbl_mb(dn(1)))
c                 end if
c
c                 k1 = k1 + 1
c                 j1 = j1 + 1
c                 if (j1.gt.i1) then
c                    j1 = 1
c                    i1 = i1 + 1
c                 end if
c              end if
c 
c              if (     ((D3dB_rc_pfft3_queue_filled()).or.(k1.gt.NN))
c     >            .and.(k2.le.NN)) then
c
c                 if (mod(k2,npj).eq.taskid_j) then
c                    call D3dB_rc_pfft3f_queueout(0,dbl_mb(dn(1)))
c
c*                   **** generate Vcoul ****
c                    !eh = coulomb_screened_e(dbl_mb(dn(1)))
c                    eh = icoulomb_screened_e(dbl_mb(dn(1)))
c                    call coulomb_screened_v(dbl_mb(dn(1)),
c     >                                      dbl_mb(vij(1)))
c
c
c*                   **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
c                    eh = eh*hfx_parameter
c                    if (ispin0.eq.1) eh = eh + eh
c                    ph = 2.0d0*eh
c                    ehfx = ehfx - eh
c                    phfx = phfx - ph
c                    dbl_mb(ehfx_orb(1,ms)+i2-1) 
c     >               = dbl_mb(ehfx_orb(1,ms)+i2-1) - eh
c                    if (i2.ne.j2) then
c                       ehfx = ehfx - eh
c                       phfx = phfx - ph
c                       dbl_mb(ehfx_orb(1,ms)+i2-1) 
c     >                  = dbl_mb(ehfx_orb(1,ms)+i2-1) - eh
c                    end if
c
c                    call D3dB_cr_pfft3b_queuein(0,dbl_mb(vij(1)))
c                 end if
c
c                 k2 = k2 + 1
c                 j2 = j2 + 1
c                 if (j2.gt.i2) then
c                    j2 = 1
c                    i2 = i2 + 1
c                 end if
c              end if
c
c              if ((D3dB_cr_pfft3_queue_filled()).or.(k2.gt.NN)) then
c
c                 if (mod(k3,npj).eq.taskid_j) then
c                    index1 =(int_mb(orbital_list(1,ms)+i3-1)-1)*n2ft3d+1
c                    index2 =(int_mb(orbital_list(1,ms)+j3-1)-1)*n2ft3d+1
c
c                    call D3dB_cr_pfft3b_queueout(0,dbl_mb(vij(1)))
c
c*                   **** apply hfx_parameter ****
cc                    call D3dB_r_SMul(1,hfx_parameter,
cc     >                                 dbl_mb(vij(1)),
cc     >                                 dbl_mb(vij(1)))
c                    call D3dB_r_SMul1(1,hfx_parameter,dbl_mb(vij(1)))
c
c*                   **** generate (Vij)*psi_r ***
c                    call D3dB_rr_Mul(1,dbl_mb(vij(1)),
c     >                                 psi_r(index2),
c     >                                 dbl_mb(tmp1(1)))
c                    call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))
c
c*                   **** add -(Vij)*psi_r to Hpsi_r ***
cc                    call D3dB_rr_Sub(1,Hpsi_r(index1),
cc     >                                 dbl_mb(tmp1(1)),
cc     >                                 Hpsi_r(index1))
c                    call D3dB_rr_Sub2(1,dbl_mb(tmp1(1)),Hpsi_r(index1))
c
c                    !**** include off diagonal terms ****
c                    if (i3.ne.j3) then
c*                      **** generate (Vij)*psi_r ***
c                       call D3dB_rr_Mul(1,dbl_mb(vij(1)),
c     >                                 psi_r(index1),
c     >                                 dbl_mb(tmp2(1)))
c                       call D3dB_r_Zero_Ends(1,dbl_mb(tmp2(1)))
c
c*                      **** add -(Vij)*psi_r to Hpsi_r ***
cc                       call D3dB_rr_Sub(1,Hpsi_r(index2),
cc     >                                 dbl_mb(tmp2(1)),
cc     >                                 Hpsi_r(index2))
c                       call D3dB_rr_Sub2(1,dbl_mb(tmp2(1)),
c     >                                     Hpsi_r(index2))
c                    end if
c                 end if
c
c                 k3 = k3 + 1
c                 j3 = j3 + 1
c                 if (j3.gt.i3) then
c                    j3 = 1
c                    i3 = i3 + 1
c                 end if
c              end if
c              done = (k1.gt.NN).and.(k2.gt.NN).and.(k3.gt.NN)
c           end do !**** while ****
c           !call D1dB_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
c           call Parallel_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
c       end do !**** ms *****
c
c
c        value =           MA_pop_stack(tmp2(2))
c        value = value.and.MA_pop_stack(tmp1(2))
c        value = value.and.MA_pop_stack(vij(2))
c        value = value.and.MA_pop_stack(dn(2))
c        if (.not. value) 
c     >    call errquit('band_potential_HFX:popping stack memory',0,
c     &       MA_ERR)

         call Parallel_SumAll(ehfx)
         call Parallel_SumAll(phfx)
      end if

      return
      end


*     *****************************
*     *                           *
*     *     band_energy_HFX_sub   *
*     *                           *
*     *****************************
      subroutine band_energy_HFX_sub(ispin0,psi_r_tag,ehfx_out,phfx_out)
      implicit none
      integer ispin0
      integer psi_r_tag
      real*8 ehfx_out
      real*8 phfx_out

#include "mafdecls.fh"
#include "band_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n1,n2,n3,ms,k1
      integer dn(2),tmp1(2),index1,index2
      real*8  scal1,scal2,dv,eh,ph

*     **** external functions ****
      real*8   lattice_omega,coulomb_screened_e
      external lattice_omega,coulomb_screened_e

      
      if (((norbs(1)+norbs(2)).ne.0).and.(.not.relaxed)) then
        ehfx = 0.0d0
        phfx = 0.0d0

c        call D3dB_nx(1,n1)
c        call D3dB_ny(1,n2)
c        call D3dB_nz(1,n3)
c        !call D3dB_n2ft3d(1,n2ft3d)
c        value = MA_push_get(mt_dbl,(2*n2ft3d),'dn_hfx',dn(2),dn(1))
c        value = value.and.
c     >          MA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
c        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
c
c        scal1 = 1.0d0/dble(n1*n2*n3)
c        scal2 = 1.0d0/lattice_omega()
c        dv = scal1/scal2
c
c        k1 = 1
c        do ms=1,ispin
c        do i=1,norbs(ms)
c         dbl_mb(ehfx_orb(1,ms)+i-1) = 0.0d0
c         do j=1,i
c
c            if (mod(k1,npj).eq.taskid_j) then
c              index1 = (int_mb(orbital_list(1,ms)+i-1)-1)*n2ft3d + 1
c              index2 = (int_mb(orbital_list(1,ms)+j-1)-1)*n2ft3d + 1
c
c*             **** generate dnij ****
c              call D3dB_rr_Mul(1,psi_r(index1),psi_r(index2),
c     >                         dbl_mb(dn(1)))
cc              call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
c              call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
c              call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))
c
c*             ***** screened coulomb solver ****
c              if (solver_type.eq.1) then
c
c*               **** generate dng ****
cc                call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
cc     >                                   dbl_mb(dn(1)))
c                call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))
c                call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
c                call Pack_c_pack(0,dbl_mb(dn(1)))
c
c*               **** get Ecoul energy ****
c                eh = coulomb_screened_e(dbl_mb(dn(1)))
c       
c*             ***** free-space coulomb solver ****
c              else
c                 call coulomb2_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))
c                 call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(tmp1(1)),eh)
c                 eh = 0.5d0*eh*dv
c              end if
c
c*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
c              eh = eh*hfx_parameter
c              if (ispin0.eq.1) eh = eh + eh
c              ph = 2.0d0*eh
c
c              ehfx = ehfx - eh
c              phfx = phfx - ph
c              dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)-eh
c
c              !**** include off diagonal terms ****
c              if (i.ne.j) then
c                 ehfx = ehfx - eh
c                 phfx = phfx - ph
c                 dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)
c     >                                      - eh
c              end if
c
c            end if
c            k1 = k1 + 1
c
c         end do
c        end do
c        call D1dB_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
c        end do
c
c        value =           MA_pop_stack(tmp1(2))
c        value = value.and.MA_pop_stack(dn(2))
c        if (.not. value) 
c     >     call errquit('band_energy_HFX_sub:popping stack memory',0,
c     &       MA_ERR)
c
c        call D1dB_SumAll(ehfx)
c        call D1dB_SumAll(phfx)
      end if

      ehfx_out = ehfx
      phfx_out = phfx

      return
      end


