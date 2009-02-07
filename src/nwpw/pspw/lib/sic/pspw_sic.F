*
* $Id: pspw_sic.F,v 1.35 2009-02-07 03:50:59 bylaska Exp $
*

*     *************************
*     *                       *
*     *     pspw_init_SIC     *
*     *                       *
*     *************************
      subroutine pspw_init_SIC(rtdb,ne)
      implicit none
      integer rtdb
      integer ne(2)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      logical value,sic_on
      integer ma_type
      integer n1,n2,n3,mapping

*     **** external functions ****
      integer  control_version,control_mapping
      external control_version,control_mapping

      norbs = 0
      e_sic = 0.0d0
      p_sic = 0.0d0
      sic_on = .false.
      if (.not.rtdb_get(rtdb,'pspw:SIC_all',mt_log,1,sicall)) 
     >   sicall=.false.

*     **** get the number of SIC orbitals ****
      if (sicall) then
         sic_on = .true.
         norbs = ne(1)+ne(2)
         if(.not.MA_alloc_get(mt_int,norbs,
     >       'orbital_list',orbital_list(2),orbital_list(1)))
     >        call errquit('pspw_init_SIC:out of heap memory',0,MA_ERR)
          
         do n1=1,norbs
            int_mb(orbital_list(1)+n1-1) = n1
         end do

      else if (rtdb_ma_get(rtdb, 'pspw:SIC_orbitals', ma_type,
     >                norbs, orbital_list(2))) then

         if (.not.MA_get_index(orbital_list(2),orbital_list(1))) 
     >      call errquit(
     >       'pspw_init_SIC: ma_get_index failed for actlist',911,
     &       MA_ERR)
         sic_on = .true.
      end if

      if (sic_on) then

         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_screening_radius',
     >                      mt_dbl,1,rcut)) 
     >       rcut = 8.0d0

         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_screening_power',
     >                      mt_dbl,1,pp)) 
     >       pp = 8.0d0

         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_screening_type',
     >                      mt_int,1,flag)) 
     >       flag = 0

         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_relax',
     >                      mt_log,1,relaxed)) 
     >       relaxed = .true.

         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_solver_type',
     >                      mt_int,1,solver_type))  then
           if (control_version().eq.3) solver_type = 1
           if (control_version().eq.4) solver_type = 2
         end if

         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_xc_parameter',
     >                       mt_dbl,1,xc_parameter))
     >       xc_parameter = 1.0d0
         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_h_parameter',
     >                       mt_dbl,1,h_parameter))
     >       h_parameter = 1.0d0
         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_print_orbital_contribution',
     >                       mt_log,1,orb_contribution))
     >       orb_contribution = .false.
 


*        **** initialize coulomb_screened ****
         if (solver_type.eq.1) then
              call coulomb_screened_init(flag,rcut,pp)

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
         if (orb_contribution) then
           value = MA_alloc_get(mt_dbl,norbs,
     >                'eh_orb',eh_orb(2),eh_orb(1))
           value = value.and.
     >             MA_alloc_get(mt_dbl,norbs,
     >                'xc_orb',xc_orb(2),xc_orb(1))
           if (.not. value)   
     >       call errquit('pspw_init_SIC: out of heap memory',1, MA_ERR)
         end if

      else
         norbs = 0
      end if

   
      return
      end


*     *************************
*     *                       *
*     *     pspw_end_SIC      *
*     *                       *
*     *************************
      subroutine pspw_end_SIC()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter(MASTER=0)
      logical value
      integer i

*     **** external functions ****
      integer  control_version
      external control_version

      if (norbs.gt.0) then


*       **** print out orbital contributions ****
        if (orb_contribution) then
           call Parallel_taskid(taskid)
           if (taskid.eq.MASTER) then
              write(6,487)
              write(6,488)
              do i=1,norbs
                write(6,489) 
     >            int_mb(orbital_list(1)+i-1),
     >            (dbl_mb(eh_orb(1)+i-1)+dbl_mb(xc_orb(1)+i-1)),
     >            dbl_mb(eh_orb(1)+i-1),
     >            dbl_mb(xc_orb(1)+i-1)
              end do
           end if
  487   format(//,'== Orbital Contributions to SIC ==')
  488   format(/1x,'orbital',15x,
     >         'Total',11x, 'Hartree',10x,'exc-corr') 	
  489   format(1x,i7,2x,e18.6,e18.6,e18.6)
        end if



*       **** deallocate memory ****
        value =  MA_free_heap(orbital_list(2)) 
        if (orb_contribution) then
          value = value.and.MA_free_heap(eh_orb(2))
          value = value.and.MA_free_heap(xc_orb(2))
        end if

        if (.not. value)
     >  call errquit('pspw_end_SIC:error freeing heap memory',0, MA_ERR)

*        **** end coulomb_screened ****
        if (solver_type.eq.1) then
          call coulomb_screened_end()

*        **** end free-space coulomb if necessary ****
        else
           if (control_version().eq.3) then
              call coulomb2_end()
              call D3dB_end(2)
           end if
        end if

      end if


      return
      end

*     *************************
*     *                       *
*     *     pspw_print_SIC    *
*     *                       *
*     *************************
      subroutine pspw_print_SIC(unit)
      implicit none
      integer unit

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      integer i

      if (norbs.ne.0) then
        if (relaxed) then
          write(unit,1001)
        else
          write(unit,1002)
        end if
        write(unit,1003) (int_mb(orbital_list(1)+i-1),i=1,norbs)

        if (solver_type.eq.1) then
          write(unit,1004)
           if (rcut.ge.0.0d0) write(unit,1006) rcut
           if (rcut.ge.0.0d0) write(unit,1007) pp
           if (rcut.ge.0.0d0) write(unit,1010) flag
        else
          write(unit,1005)
        end if
        if (h_parameter.ne.1.0d0) write(unit,1008) h_parameter
        if (xc_parameter.ne.1.0d0)write(unit,1009) xc_parameter
        write(unit,*)

      end if

      return
 1001 FORMAT(6x,"- SIC relaxed")
 1002 FORMAT(6x,"- SIC unrelaxed")
 1003 FORMAT(6x,"- SIC orbitals:",10I5)
 1004 FORMAT(6x,"- SIC screened coulomb solver")
 1005 FORMAT(6x,"- SIC free-space coulomb solver")
 1006 FORMAT(6x,"- SIC screening radius(pspw:SIC_screening_radius):",
     >       E10.3)
 1007 FORMAT(6x,"- SIC screening power (pspw:SIC_screening_power) :",
     >       E10.3)
 1008 FORMAT(6x,"- SIC h_parameter     (pspw:SIC_h_parameter)     :",
     >       E10.3)
 1009 FORMAT(6x,"- SIC xc_parameter    (pspw:SIC_xc_parameter)    :",
     >       E10.3)
 1010 FORMAT(6x,"- SIC screening type (pspw:SIC_screening_type)   :",
     >       I2)
      end

*     ****************************
*     *                    	 *
*     *     pspw_potential_SIC   *
*     *                          *
*     ****************************
      subroutine pspw_potential_SIC(ispin,psi_r,Hpsi_r)
      implicit none
#include "errquit.fh"
      integer    ispin
      real*8     psi_r(*)
      real*8     Hpsi_r(*)

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      logical value
      integer i,n2ft3d,n1,n2,n3,npack,gga
      integer dn(2),xcp(2),xce(2),tmp1(2),index
      !integer index2
      real*8  scal1,scal2,dv,exc,pxc,eh,ph

*     **** external functions ****
      integer  control_gga
      real*8   lattice_omega,coulomb_screened_e
      external control_gga
      external lattice_omega,coulomb_screened_e


      eh_sic = 0.0d0
      ph_sic = 0.0d0
      ex_sic = 0.0d0
      px_sic = 0.0d0
      if ((norbs.ne.0).and.relaxed) then

        gga = control_gga()
        call Pack_npack(1,npack)
        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        call D3dB_n2ft3d(1,n2ft3d)
        value = MA_push_get(mt_dbl,(2*n2ft3d),'dn_SIC',dn(2),dn(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(2*n2ft3d),'xcpSC',xcp(2),xcp(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(2*n2ft3d),'xceSC',xce(2),xce(1))
        value = value.and.
     >             MA_push_get(mt_dbl,(n2ft3d),'tmp1C',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2
        do i=1,norbs

           index  = (int_mb(orbital_list(1)+i-1)-1)*n2ft3d + 1
           !index2 = (int_mb(orbital_list(1)+i-1)-1)*npack  + 1

*          **** generate dn for Vxc and VH ****
           call D3dB_rr_Sqr(1,psi_r(index),dbl_mb(dn(1)))
c           call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
           call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))


*          **** generate exchange-correlation potential ****
           if (gga.eq.0) then
             call vxc(n2ft3d,2,dbl_mb(dn(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)),
     >                      dbl_mb(tmp1(1)))
           else
             call v_bwexc(gga,n2ft3d,2,dbl_mb(dn(1)),
     >                    1.0d0,1.0d0,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)))
           end if

*          **** apply the xc_parameter ****
c           call D3dB_r_SMul(1,xc_parameter,
c     >                      dbl_mb(xcp(1)),
c     >                      dbl_mb(xcp(1)))
c           call D3dB_r_SMul(1,xc_parameter,
c     >                      dbl_mb(xce(1)),
c     >                      dbl_mb(xce(1)))
           call D3dB_r_SMul1(1,xc_parameter,dbl_mb(xcp(1)))
           call D3dB_r_SMul1(1,xc_parameter,dbl_mb(xce(1)))

*          **** get exchange-correlation energy ****
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xce(1)),exc)
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xcp(1)),pxc)
           if (ispin.eq.1) exc= exc + exc
           if (ispin.eq.1) pxc= pxc + pxc
           exc = exc*dv
           pxc = pxc*dv

           ex_sic = ex_sic - exc
           px_sic = px_sic - pxc


*          ***** screened coulomb solver ****
           if (solver_type.eq.1) then

*            **** generate dng ****
c            if (ispin.eq.1) 
c    >         call D3dB_rr_Sum(1,dbl_mb(dn(1)),
c    >                            dbl_mb(dn(1)),
c    >                            dbl_mb(dn(1)))

c             call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
c     >                              dbl_mb(dn(1)))
             call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))

             !call D3dB_rc_fft3f(1,dbl_mb(dn(1)))
             call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
             call Pack_c_pack(0,dbl_mb(dn(1)))

*            **** get Ecoul energy ****
             eh = coulomb_screened_e(dbl_mb(dn(1)))
             if (ispin.eq.1) eh = eh + eh
             ph = 2.0d0*eh

*            **** generate Vcoul ****
             call coulomb_screened_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))
             call Pack_c_unpack(0,dbl_mb(tmp1(1)))
             !call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
             call D3dB_cr_pfft3b(1,0,dbl_mb(tmp1(1)))



*          ***** free-space coulomb solver ****
           else
c             if (ispin.eq.1) 
c     >         call D3dB_rr_Sum(1,dbl_mb(dn(1)),
c     >                            dbl_mb(dn(1)),
c     >                            dbl_mb(dn(1)))

*             ** perform shift here? ***
              call coulomb2_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))
              call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(tmp1(1)),eh)
*             ** perform reverse shift here? ***

              eh = 0.5d0*eh*dv
              if (ispin.eq.1) eh = eh + eh
              ph = 2.0d0*eh
           end if

*          **** apply the h_parameter ****
           eh = eh*h_parameter
           ph = ph*h_parameter
c           call D3dB_r_SMul(1,h_parameter,
c     >                      dbl_mb(tmp1(1)),
c     >                      dbl_mb(tmp1(1)))
           call D3dB_r_SMul1(1,h_parameter,dbl_mb(tmp1(1)))

           eh_sic = eh_sic - eh
           ph_sic = ph_sic - ph

*          **** generate (Vcoul+Vxc)*psi_r ***
c           call D3dB_rr_Sum(1,dbl_mb(tmp1(1)),
c     >                        dbl_mb(xcp(1)),
c     >                        dbl_mb(xcp(1)))
           call D3dB_rr_Sum2(1,dbl_mb(tmp1(1)),dbl_mb(xcp(1)))

           call D3dB_rr_Mul(1,dbl_mb(xcp(1)),
     >                        psi_r(index),
     >                        dbl_mb(tmp1(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))

*          **** add -(Vcoul+Vxc)*psi_r to Hpsi_r ***
c           call D3dB_rr_Sub(1,Hpsi_r(index),
c     >                        dbl_mb(tmp1(1)),
c     >                        Hpsi_r(index))
           call D3dB_rr_Sub2(1,dbl_mb(tmp1(1)),Hpsi_r(index))

*          **** fourier transform (Vcoul+Vxc)*psi_r and add to Hpsi_k ***
c          call D3dB_rc_fft3f(1,dbl_mb(tmp1(1)))
c          call Pack_c_pack(1,dbl_mb(tmp1(1)))
c          call Pack_cc_daxpy(1,(scal1),
c    >                        dbl_mb(tmp1(1)),
c    >                        Hpsi_k(index2))

*          **** keep track of orbital contributions ****
           if (orb_contribution) then
             dbl_mb(eh_orb(1)+i-1) = -eh
             dbl_mb(xc_orb(1)+i-1) = -exc

           end if
        end do

        value = MA_pop_stack(tmp1(2))
        value = value.and.MA_pop_stack(xce(2))
        value = value.and.MA_pop_stack(xcp(2))
        value = value.and.MA_pop_stack(dn(2))
        if (.not. value) 
     >    call errquit('pspw_potential_SIC:popping stack memory',0,
     &       MA_ERR)
      end if
      return
      end


*     *************************
*     *                       *
*     *     pspw_energy_SIC   *
*     *                       *
*     *************************
      subroutine pspw_energy_SIC(ispin,psi_r,ehsic,phsic,exsic,pxsic)
      implicit none
#include "errquit.fh"
      integer ispin
      real*8  psi_r(*)
      real*8 ehsic,phsic
      real*8 exsic,pxsic

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      logical value
      integer i,n2ft3d,n1,n2,n3,npack,gga
      integer dn(2),xcp(2),xce(2),tmp1(2),index
      !integer index2
      real*8  scal1,scal2,dv,exc,pxc,eh,ph

*     **** external functions ****
      integer  control_gga
      real*8   lattice_omega,coulomb_screened_e
      external control_gga
      external lattice_omega,coulomb_screened_e

      ehsic = eh_sic
      phsic = ph_sic
      exsic = ex_sic
      pxsic = px_sic
      if ((norbs.ne.0).and.(.not.relaxed)) then
        ehsic = 0.0d0
        phsic = 0.0d0
        exsic = 0.0d0
        pxsic = 0.0d0

        gga = control_gga()
        call Pack_npack(1,npack)
        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        call D3dB_n2ft3d(1,n2ft3d)
        value = MA_push_get(mt_dbl,(2*n2ft3d),'dn_SIC',dn(2),dn(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(2*n2ft3d),'xcpSC',xcp(2),xcp(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(2*n2ft3d),'xceSC',xce(2),xce(1))
        value = value.and.
     >             MA_push_get(mt_dbl,(n2ft3d),'tmp1C',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2
        do i=1,norbs

           index  = (int_mb(orbital_list(1)+i-1)-1)*n2ft3d + 1
           !index2 = (int_mb(orbital_list(1)+i-1)-1)*npack  + 1

*          **** generate dn ****
           call D3dB_rr_Sqr(1,psi_r(index),dbl_mb(dn(1)))
c           call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
           call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))


*          **** generate exchange-correlation potential ****
           if (gga.eq.0) then
             call vxc(n2ft3d,2,dbl_mb(dn(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)),
     >                      dbl_mb(tmp1(1)))
           else
             call v_bwexc(gga,n2ft3d,2,dbl_mb(dn(1)),
     >                    1.0d0,1.0d0,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)))
           end if

*          **** apply the xc_parameter ****
c           call D3dB_r_SMul(1,xc_parameter,
c     >                      dbl_mb(xcp(1)),
c     >                      dbl_mb(xcp(1)))
c           call D3dB_r_SMul(1,xc_parameter,
c     >                      dbl_mb(xce(1)),
c     >                      dbl_mb(xce(1)))
           call D3dB_r_SMul1(1,xc_parameter,dbl_mb(xcp(1)))
           call D3dB_r_SMul1(1,xc_parameter,dbl_mb(xce(1)))

*          **** get exchange-correlation energy ****
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xce(1)),exc)
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xcp(1)),pxc)
           if (ispin.eq.1) exc= exc + exc
           if (ispin.eq.1) pxc= pxc + pxc
           exc = exc*dv
           pxc = pxc*dv

           exsic = exsic - exc
           pxsic = pxsic - pxc


*          ***** screened coulomb solver ****
           if (solver_type.eq.1) then

*            **** generate dng ****
c             call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
c     >                                dbl_mb(dn(1)))
             call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))
             !call D3dB_rc_fft3f(1,dbl_mb(dn(1)))
             call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
             call Pack_c_pack(0,dbl_mb(dn(1)))

*            **** get Ecoul energy ****
             eh = coulomb_screened_e(dbl_mb(dn(1)))
             eh = eh*h_parameter
             if (ispin.eq.1) eh = eh + eh
             ph = 2.0d0*eh

*          ***** free-space coulomb solver ****
           else
              call coulomb2_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))

              call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(tmp1(1)),eh)
              eh = 0.5d0*eh*dv
              eh = eh*h_parameter
              if (ispin.eq.1) eh = eh + eh
              ph = 2.0d0*eh
           end if

           ehsic = ehsic - eh
           phsic = phsic - ph

*          **** keep track of orbital contributions ****
           if (orb_contribution) then
             dbl_mb(eh_orb(1)+i-1) = -eh
             dbl_mb(xc_orb(1)+i-1) = -exc
           end if

          end do

        value = MA_pop_stack(tmp1(2))
        value = value.and.MA_pop_stack(xce(2))
        value = value.and.MA_pop_stack(xcp(2))
        value = value.and.MA_pop_stack(dn(2))
        if (.not. value) 
     >     call errquit('pspw_energy_SIC:popping stack memory',0,
     &       MA_ERR)

      end if


      return
      end

*     *************************
*     *                       *
*     *     pspw_SIC          *
*     *                       *
*     *************************
      logical function pspw_SIC()
      implicit none

#include "pspw_sic.fh"

      pspw_SIC= (norbs.gt.0)
      return
      end

*     *************************
*     *                       *
*     *   pspw_SIC_relaxed    *
*     *                       *
*     *************************
      logical function pspw_SIC_relaxed()
      implicit none

#include "mafdecls.fh"
#include "pspw_sic.fh"

      pspw_SIC_relaxed = relaxed
      return
      end


*     ********************************
*     *                      	     *
*     *     pspw_potential_SIC_OEP   *
*     *                              *
*     ********************************

      subroutine pspw_potential_SIC_OEP(ispin,ne,dntotal,psi_r,v_oep)
      implicit none
#include "errquit.fh"
      integer    ispin,ne(2)
      real*8     dntotal(*)
      real*8     psi_r(*)
      real*8     v_oep(*)

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      logical value
      integer i,j,n2ft3d,nx,ny,nz,npack,n1(2),n2(2),gga
      integer index,index1,indext,ms
      integer indexi,indexj
      integer dn(2),xcp(2),xce(2),tmp1(2),vsic(2),rhoi(2)
      integer M(2),b(2),ipv(2),info
      real*8  scal1,scal2,dv,exc,pxc,eh,ph,xmax(2)

*     **** external functions ****
      integer  control_gga
      real*8   lattice_omega,coulomb_screened_e
      external control_gga
      external lattice_omega,coulomb_screened_e


      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      
        gga = control_gga()
        call Pack_npack(1,npack)
        call D3dB_nx(1,nx)
        call D3dB_ny(1,ny)
        call D3dB_nz(1,nz)
        call D3dB_n2ft3d(1,n2ft3d)
        value = MA_push_get(mt_dbl,(2*n2ft3d),'dn_SIC',dn(2),dn(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(2*n2ft3d),'xcpSC',xcp(2),xcp(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(2*n2ft3d),'xceSC',xce(2),xce(1))
        value = value.and.
     >             MA_push_get(mt_dbl,(n2ft3d),'tmp1C',tmp1(2),tmp1(1))
        value = value.and.
     >             MA_push_get(mt_dbl,((ne(1)+ne(2))*n2ft3d),
     >                         'rhoi',rhoi(2),rhoi(1))
        value = value.and.
     >             MA_push_get(mt_dbl,((ne(1)+ne(2))*n2ft3d),
     >                         'vsic',vsic(2),vsic(1))
        value = value.and.
     >             MA_push_get(mt_dbl,((ne(1)*ne(1))*2),
     >                         'M',M(2),M(1))
        value = value.and.
     >             MA_push_get(mt_dbl,((ne(1)+ne(2))),
     >                         'b',b(2),b(1))
        if (.not. value) 
     >    call errquit('pspw_potential_SIC_OEP:out of stack memory',0,0)
        call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)


        scal1 = 1.0d0/dble(nx*ny*nz)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2
        do i=1,ne(1)+ne(2)

           index  = (i-1)*n2ft3d + 1

*          **** generate dn for Vxc and VH ****
           call D3dB_rr_Sqr(1,psi_r(index),dbl_mb(rhoi(1)+(i-1)*n2ft3d))
c           call D3dB_r_SMul(1,scal2,
c     >                      dbl_mb(rhoi(1)+(i-1)*n2ft3d),
c     >                      dbl_mb(rhoi(1)+(i-1)*n2ft3d))
           call D3dB_r_SMul1(1,scal2,dbl_mb(rhoi(1)+(i-1)*n2ft3d))
           call D3dB_r_Zero_Ends(1,dbl_mb(rhoi(1)+(i-1)*n2ft3d))
           call dcopy(n2ft3d,dbl_mb(rhoi(1)+(i-1)*n2ft3d),1,
     >                        dbl_mb(dn(1)),1)


*          **** generate exchange-correlation potential ****
           if (gga.eq.0) then
             call vxc(n2ft3d,2,
     >                      dbl_mb(dn(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)),
     >                      dbl_mb(tmp1(1)))
           else
             call v_bwexc(gga,n2ft3d,2,dbl_mb(dn(1)),
     >                    1.0d0,1.0d0,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)))
           end if

*          **** apply the xc_parameter ****
c           call D3dB_r_SMul(1,xc_parameter,
c     >                      dbl_mb(xcp(1)),
c     >                      dbl_mb(xcp(1)))
c           call D3dB_r_SMul(1,xc_parameter,
c     >                      dbl_mb(xce(1)),
c     >                      dbl_mb(xce(1)))
           call D3dB_r_SMul1(1,xc_parameter,dbl_mb(xcp(1)))
           call D3dB_r_SMul1(1,xc_parameter,dbl_mb(xce(1)))


*          ***** screened coulomb solver ****
           if (solver_type.eq.1) then

*            **** generate dng ****

             call D3dB_r_SMul(1,scal1,dbl_mb(rhoi(1)+(i-1)*n2ft3d),
     >                              dbl_mb(dn(1)))
             !call D3dB_rc_fft3f(1,dbl_mb(dn(1)))
             call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
             call Pack_c_pack(0,dbl_mb(dn(1)))

*            **** generate Vcoul ****
             call coulomb_screened_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))
             call Pack_c_unpack(0,dbl_mb(tmp1(1)))
             !call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
             call D3dB_cr_pfft3b(1,0,dbl_mb(tmp1(1)))


*          ***** free-space coulomb solver ****
           else
              call coulomb2_v(dbl_mb(rhoi(1)+(i-1)*n2ft3d),
     >                        dbl_mb(tmp1(1)))
              
           end if

*          **** apply the h_parameter ****
c           call D3dB_r_SMul(1,h_parameter,
c     >                      dbl_mb(tmp1(1)),
c     >                      dbl_mb(tmp1(1)))
           call D3dB_r_SMul1(1,h_parameter,dbl_mb(tmp1(1)))


*          **** generate vsic(i) = (Vcoul+Vxc) ***
           call D3dB_rr_Sum(1,dbl_mb(tmp1(1)),
     >                        dbl_mb(xcp(1)),
     >                        dbl_mb(vsic(1)+(i-1)*n2ft3d))

c           call D3dB_r_SMul(1,(-1.0d0),
c     >                        dbl_mb(vsic(1)+(i-1)*n2ft3d),
c     >                        dbl_mb(vsic(1)+(i-1)*n2ft3d))
           call D3dB_r_SMul1(1,(-1.0d0),dbl_mb(vsic(1)+(i-1)*n2ft3d))


        end do

        !**** define Vs  ****
        call dcopy(ispin*n2ft3d,0.0d0,0,v_oep,1)

        do ms=1,ispin
        do i=n1(ms),n2(ms)

           call D3dB_rr_Mul(1,dbl_mb(rhoi(1)+(i-1)*n2ft3d),
     >                        dbl_mb(vsic(1)+(i-1)*n2ft3d),
     >                        dbl_mb(tmp1(1)))
c           call D3dB_rr_Sum(1,v_oep((ms-1)*n2ft3d+1),
c     >                        dbl_mb(tmp1(1)),
c     >                        v_oep((ms-1)*n2ft3d+1))
           call D3dB_rr_Sum2(1,dbl_mb(tmp1(1)),
     >                        v_oep((ms-1)*n2ft3d+1))
        end do
        end do

        do ms=1,ispin
c          call D3dB_rr_Divide(1,  v_oep((ms-1)*n2ft3d+1),
c     >                          dntotal((ms-1)*n2ft3d+1),
c     >                            v_oep((ms-1)*n2ft3d+1))
          call D3dB_rr_Divide2(1,dntotal((ms-1)*n2ft3d+1),
     >                             v_oep((ms-1)*n2ft3d+1))
        end do

        !**** define (-M)  ****
        do ms=1,ispin
        do i=1,ne(ms)
        do j=i,ne(ms)
           indexi = (i-1 + (ms-1)*ne(1))*n2ft3d
           indexj = (j-1 + (ms-1)*ne(1))*n2ft3d
           call D3dB_rr_Mul(1,dbl_mb(rhoi(1)+indexi),
     >                        dbl_mb(rhoi(1)+indexj),
     >                        dbl_mb(tmp1(1)))

c           call D3dB_rr_Divide(1,dbl_mb(tmp1(1)),
c     >                          dntotal((ms-1)*n2ft3d+1),
c     >                           dbl_mb(tmp1(1)))
           call D3dB_rr_Divide2(1,dntotal((ms-1)*n2ft3d+1),
     >                           dbl_mb(tmp1(1)))

           call D3dB_r_dsum(1,dbl_mb(tmp1(1)),eh)
           index  = (i-1) + (j-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
           indext = (j-1) + (i-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
           dbl_mb(M(1)+index)  = -eh*dv
           dbl_mb(M(1)+indext) = -eh*dv
        
        end do
        end do
        end do


        !**** M = I-M ****
        do ms=1,ispin
        do i=1,ne(ms)

           index = (i-1) + (i-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
           dbl_mb(M(1)+index)  = 1.0d0 + dbl_mb(M(1)+index)
        end do
        end do



        !**** define b=(vs-vsic)  ****
        do ms=1,ispin
        do i=n1(ms),n2(ms)

          call D3dB_rr_dot(1,dbl_mb(rhoi(1)+(i-1)*n2ft3d),
     >                       dbl_mb(vsic(1)+(i-1)*n2ft3d),
     >                       exc)
          call D3dB_rr_dot(1,dbl_mb(rhoi(1)+(i-1)*n2ft3d),
     >                       v_oep((ms-1)*n2ft3d+1),
     >                       pxc)
          index = (i-1) 
          dbl_mb(b(1)+index) = (pxc-exc)*dv
        end do
        end do

*       ***** Check M and b *****
        write(*,*)
        write(*,*)
        do ms=1,ispin
        do j=1,ne(ms)
          pxc = 0.0d0
          do i=1,ne(ms)
           index = (i-1) + (j-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
           pxc = pxc + dbl_mb(M(1)+index)
          end do
          write(*,*) "Sum M:",ms,j,pxc
        end do
        end do
        write(*,*)
        do ms=1,ispin
          pxc = 0.0d0
          do i=n1(ms),n2(ms)
           pxc = pxc + dbl_mb(b(1)+i-1)
          end do
          write(*,*) "Sum b:",ms,pxc
        end do
        write(*,*)
        write(*,*)


        !*** solve Ax = b ***
        value = MA_push_get(mt_int,(ne(1)+ne(2)),
     >                      'ipiv',ipv(2),ipv(1))
        if (.not. value) 
     >    call errquit('pspw_potential_SIC_OEP:out of stack memory',1,0)
        do ms=1,ispin
          call DGESV(ne(ms),1, 
     >             dbl_mb(M(1)+ (ms-1)*ne(1)*ne(1)), ne(ms),
     >             int_mb(ipv(1)),
     >             dbl_mb(b(1)+(ms-1)*ne(1)),ne(ms),
     >             info)
        end do
        value = MA_pop_stack(ipv(2))
        if (.not. value) 
     >    call errquit(
     >    'pspw_potential_SIC_OEP:popping stack memory',2,0)

        do ms=1,ispin
        do i=1,ne(ms)
          write(*,*) "x:",ms,i,dbl_mb(b(1)+i-1+(ms-1)*ne(1))
        end do
        write(*,*)
        end do

        !*** determine xmax ***
        do ms=1,ispin
         xmax(ms) = dbl_mb(b(1)+ (ms-1)*ne(1))
         do i=2,ne(ms)
           if (dbl_mb(b(1)+i-1+(ms-1)*ne(1)).gt.xmax(ms))
     >       xmax(ms) = dbl_mb(b(1)+i-1+(ms-1)*ne(1))
         end do
        end do


        !***** tmp = (1/rho)* Sum(i=2,N) rhoi * (xi - xmax) *****
        do ms=1,ispin
          call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)
          do i=n1(ms),n2(ms)
           index  = (i-1) 
           index1 = (ms-1)*ne(1)
           eh     = dbl_mb(b(1)+index) - xmax(ms)
           call D3dB_rr_daxpy(1,eh,
     >                        dbl_mb(rhoi(1)+(i-1)*n2ft3d),
     >                        dbl_mb(tmp1(1)))

          end do
c          call D3dB_rr_Divide(1,dbl_mb(tmp1(1)),
c     >                        dntotal((ms-1)*n2ft3d+1),
c     >                          dbl_mb(tmp1(1)))
          call D3dB_rr_Divide2(1,dntotal((ms-1)*n2ft3d+1),
     >                          dbl_mb(tmp1(1)))

c          call D3dB_rr_Sum(1,v_oep((ms-1)*n2ft3d+1),
c     >                       dbl_mb(tmp1(1)),
c     >                       v_oep((ms-1)*n2ft3d+1))
          call D3dB_rr_Sum2(1,dbl_mb(tmp1(1)),
     >                       v_oep((ms-1)*n2ft3d+1))
        end do

      
        value =           MA_pop_stack(b(2))
        value = value.and.MA_pop_stack(M(2))
        value = value.and.MA_pop_stack(vsic(2))
        value = value.and.MA_pop_stack(rhoi(2))
        value = value.and.MA_pop_stack(tmp1(2))
        value = value.and.MA_pop_stack(xce(2))
        value = value.and.MA_pop_stack(xcp(2))
        value = value.and.MA_pop_stack(dn(2))
        if (.not. value) 
     >    call errquit(
     >     'pspw_potential_SIC_OEP:popping stack memory',3,0)

      return
      end



*     *************************
*     *                       *
*     *     pspw_SIC_euv      *
*     *                       *
*     *************************
      subroutine pspw_SIC_euv(ispin,psi_r,stress)
      implicit none
      integer ispin
      real*8  psi_r(*)
      real*8 stress(3,3)

#include "mafdecls.fh"
#include "pspw_sic.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,n2ft3d,n1,n2,n3,npack,u,v,gga
      integer dn(2),xcp(2),xce(2),tmp1(2),index
      !integer index2
      real*8  scal1,scal2,dv,exc,pxc,eh,ph
      real*8 hm(3,3),scal,pi,tstress(3,3)

*     **** external functions ****
      integer  control_gga
      real*8   lattice_omega,lattice_unitg
      external control_gga
      external lattice_omega,lattice_unitg

      call dcopy(9,0.0d0,0,stress,1)
      if ((norbs.ne.0)) then
        call coulomb_screened_euv_init(flag,rcut,pp)

*       *** define hm ****
        pi   = 4.0d0*datan(1.0d0)
        scal = 1.0d0/(2.0d0*pi)
        do v=1,3
        do u=1,3
           hm(u,v) = scal*lattice_unitg(u,v)
        end do
        end do

        gga = control_gga()
        call Pack_npack(1,npack)
        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        call D3dB_n2ft3d(1,n2ft3d)
        value = MA_push_get(mt_dbl,(2*n2ft3d),'dn_SIC',dn(2),dn(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(2*n2ft3d),'xcpSC',xcp(2),xcp(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(2*n2ft3d),'xceSC',xce(2),xce(1))
        value = value.and.
     >             MA_push_get(mt_dbl,(n2ft3d),'tmp1C',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
        call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2
        do i=1,norbs

           index  = (int_mb(orbital_list(1)+i-1)-1)*n2ft3d + 1
           !index2 = (int_mb(orbital_list(1)+i-1)-1)*npack  + 1

*          **** generate dn ****
           call D3dB_rr_Sqr(1,psi_r(index),dbl_mb(dn(1)))
c           call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
           call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))


*          **** generate exchange-correlation potential ****
           if (gga.eq.0) then
             call vxc(n2ft3d,2,dbl_mb(dn(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)),
     >                      dbl_mb(tmp1(1)))
           else
             call v_bwexc(gga,n2ft3d,2,dbl_mb(dn(1)),
     >                    1.0d0,1.0d0,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)))
           end if


*          **** get exchange-correlation energy ****
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xce(1)),exc)
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xcp(1)),pxc)
           if (ispin.eq.1) exc= exc + exc
           if (ispin.eq.1) pxc= pxc + pxc
           exc = exc*dv
           pxc = pxc*dv


*          **** LDA part ****
           do v=1,3
           do u=1,3
              stress(u,v) = stress(u,v) - xc_parameter*(exc-pxc)*hm(u,v)
           end do
           end do


*          **** PBE96 GGA part ****
           if ((gga.eq.1).or.(gga.eq.2)) then
            call v_bwexc_euv(gga,n2ft3d,ispin,dbl_mb(dn(1)),
     >                       1.0d0,1.0d0,tstress)
            do v=1,3
            do u=1,3
               stress(u,v) = stress(u,v) - xc_parameter*tstress(u,v)
            end do
            end do
           end if


*          ***** screened coulomb solver ****
           if (solver_type.eq.1) then

*            **** generate dng ****
c             call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
c     >                                dbl_mb(dn(1)))
             call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))
             !call D3dB_rc_fft3f(1,dbl_mb(dn(1)))
             call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
             call Pack_c_pack(0,dbl_mb(dn(1)))

*            **** get Ecoul energy ****
             call coulomb_screened_euv(dbl_mb(dn(1)),tstress)
             if (ispin.eq.1) call dscal(9,2.0d0,tstress,1)

*            **** apply the h_parameter ****
             call dscal(9,h_parameter,tstress,1)

*          ***** free-space coulomb solver ****
           else
             write(*,*) "ERROR free-space coulomb solver called"
             call errquit('error: not periodic boundary conditions',0,0)
           end if
           
           do v=1,3
           do u=1,3
              stress(u,v) = stress(u,v) - tstress(u,v)
           end do
           end do

        end do

        call coulomb_screened_euv_end()
        value = MA_pop_stack(tmp1(2))
        value = value.and.MA_pop_stack(xce(2))
        value = value.and.MA_pop_stack(xcp(2))
        value = value.and.MA_pop_stack(dn(2))
        if (.not. value) 
     >     call errquit('pspw_SIC_euv:popping stack memory',0,
     &       MA_ERR)

      end if


      return
      end

