*     *************************
*     *                       *
*     *     pspw_init_SIC     *
*     *                       *
*     *************************
      subroutine pspw_init_SIC(ne,rtdb)
      implicit none
      integer ne(2)
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      integer ma_type

      norbs = 0
      e_sic = 0.0d0
      p_sic = 0.0d0

*     **** get the number of SIC orbitals ****
      if (rtdb_ma_get(rtdb, 'pspw:SIC_orbitals', ma_type,
     >                norbs, orbital_list(2))) then

         if (.not.MA_get_index(orbital_list(2),orbital_list(1))) 
     >      call errquit(
     >       'pspw_init_SIC: ma_get_index failed for actlist',911)

         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_screening_radius',
     >                      mt_dbl,1,rcut)) 
     >       rcut = -1.0d0

         if (.not. rtdb_get(rtdb,
     >                      'pspw:SIC_relax',
     >                      mt_log,1,relaxed)) 
     >       relaxed = .false.

*        **** initialize coulomb_screened ****
         if (rcut.lt.0.0d0) then
            call coulomb_screened_init(99999.0d9)
         else
            call coulomb_screened_init(rcut)
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

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      logical value

      if (norbs.gt.0) then
        value =  MA_free_heap(orbital_list(2)) 
        if (.not. value)
     >  call errquit('pspw_end_SIC:error freeing heap memory',0)

        call coulomb_screened_end()
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
          write(unit,*) "      - SIC relaxed"
        else
          write(unit,*) "      - SIC unrelaxed"
        end if
        write(unit,*) "      - SIC orbitals:     ",
     >      (int_mb(orbital_list(1)+i-1),i=1,norbs)
        if (rcut.ge.0.0d0)
     >  write(unit,'(A,E8.3)') "       - SIC screening radius: ",rcut
        write(unit,*)
      end if

      return
      end

*     ****************************
*     *                    	     *
*     *     pspw_potential_SIC   *
*     *                          *
*     ****************************
      subroutine pspw_potential_SIC(ispin,psi_r,Hpsi_k)
      implicit none
      integer    ispin
      real*8     psi_r(*)
      complex*16 Hpsi_k(*)

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      logical value
      integer i,n2ft3d,n1,n2,n3,npack
      integer dn(2),xcp(2),xce(2),tmp1(2),index,index2
      real*8  scal1,scal2,dv,exc,pxc,eh,ph

*     **** external functions ****
      integer  control_gga
      real*8   lattice_omega,coulomb_screened_e
      external control_gga
      external lattice_omega,coulomb_screened_e


      e_sic = 0.0d0
      p_sic = 0.0d0
      if ((norbs.ne.0).and.relaxed) then
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
        if (.not. value) call errquit('out of stack memory',0)
        call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2
        do i=1,norbs

           index  = (int_mb(orbital_list(1)+i-1)-1)*n2ft3d + 1
           index2 = (int_mb(orbital_list(1)+i-1)-1)*npack  + 1

*          **** generate dn ****
           call D3dB_rr_Sqr(1,psi_r(index),dbl_mb(dn(1)))
           call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))


*          **** generate exchange-correlation potential ****
           if (control_gga().eq.0) then
             call vxc(n2ft3d,ispin,dbl_mb(dn(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)),
     >                      dbl_mb(tmp1(1)))
           else
             call v_bwexc(n2ft3d,ispin,dbl_mb(dn(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)))
           end if

*          **** get exchange-correlation energy ****
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xce(1)),exc)
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xcp(1)),pxc)
           if (ispin.eq.1) exc= exc + exc
           if (ispin.eq.1) pxc= pxc + pxc
           exc = exc*dv
           pxc = pxc*dv

           e_sic = e_sic - exc
           p_sic = p_sic - pxc


*          **** generate dng ****
           if (ispin.eq.1) 
     >       call D3dB_rr_Sum(1,dbl_mb(dn(1)),
     >                        dbl_mb(dn(1)),
     >                        dbl_mb(dn(1)))
           call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
     >                              dbl_mb(dn(1)))
           call D3dB_rc_fft3f(1,dbl_mb(dn(1)))
           call Pack_c_pack(0,dbl_mb(dn(1)))

*          **** get Ecoul energy ****
           eh = coulomb_screened_e(dbl_mb(dn(1)))
           ph = 2.0d0*eh
           e_sic = e_sic + eh
           p_sic = p_sic + ph


*          **** generate Vcoul ****
           call coulomb_screened_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))
           call Pack_c_unpack(0,dbl_mb(tmp1(1)))
           call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))


*          **** generate (Vcoul-Vxc)*psi_r ***
           call D3dB_rr_Sub(1,dbl_mb(tmp1(1)),
     >                        dbl_mb(xcp(1)),
     >                        dbl_mb(xcp(1)))
           call D3dB_rr_Mul(1,dbl_mb(xcp(1)),
     >                        psi_r(index),
     >                        dbl_mb(tmp1(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))

*          **** add to Hpsi_k ***
           call D3dB_rc_fft3f(1,dbl_mb(tmp1(1)))
           call Pack_c_pack(1,dbl_mb(tmp1(1)))
           call Pack_cc_daxpy(1,(-scal1),
     >                        dbl_mb(tmp1(1)),
     >                        Hpsi_k(index2))

        end do

        value = MA_pop_stack(tmp1(2))
        value = value.and.MA_pop_stack(xce(2))
        value = value.and.MA_pop_stack(xcp(2))
        value = value.and.MA_pop_stack(dn(2))
        if (.not. value) call errquit('popping stack memory',0)
      end if
      return
      end

*     *************************
*     *                       *
*     *     pspw_energy_SIC   *
*     *                       *
*     *************************
      subroutine pspw_energy_SIC(ispin,psi_r,esic,psic)
      implicit none
      integer ispin
      real*8  psi_r(*)
      real*8 esic,psic

#include "mafdecls.fh"
#include "pspw_sic.fh"

*     **** local variables ****
      logical value
      integer i,n2ft3d,n1,n2,n3,npack
      integer dn(2),xcp(2),xce(2),tmp1(2),index,index2
      real*8  scal1,scal2,dv,exc,pxc,eh,ph

*     **** external functions ****
      integer  control_gga
      real*8   lattice_omega,coulomb_screened_e
      external control_gga
      external lattice_omega,coulomb_screened_e

      if ((norbs.ne.0).and.(.not.relaxed)) then
        e_sic = 0.0d0
        p_sic = 0.0d0
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
        if (.not. value) call errquit('out of stack memory',0)
        call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dn(1)),1)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2
        do i=1,norbs

           index  = (int_mb(orbital_list(1)+i-1)-1)*n2ft3d + 1
           index2 = (int_mb(orbital_list(1)+i-1)-1)*npack  + 1

*          **** generate dn ****
           call D3dB_rr_Sqr(1,psi_r(index),dbl_mb(dn(1)))
           call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))


*          **** generate exchange-correlation potential ****
           if (control_gga().eq.0) then
             call vxc(n2ft3d,ispin,dbl_mb(dn(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)),
     >                      dbl_mb(tmp1(1)))
           else
             call v_bwexc(n2ft3d,ispin,dbl_mb(dn(1)),
     >                      dbl_mb(xcp(1)),
     >                      dbl_mb(xce(1)))
           end if

*          **** get exchange-correlation energy ****
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xce(1)),exc)
           call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(xcp(1)),pxc)
           if (ispin.eq.1) exc= exc + exc
           if (ispin.eq.1) pxc= pxc + pxc
           exc = exc*dv
           pxc = pxc*dv

           e_sic = e_sic - exc
           p_sic = p_sic - pxc


*          **** generate dng ****
           if (ispin.eq.1)
     >       call D3dB_rr_Sum(1,dbl_mb(dn(1)),
     >                        dbl_mb(dn(1)),
     >                        dbl_mb(dn(1)))
           call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
     >                              dbl_mb(dn(1)))
           call D3dB_rc_fft3f(1,dbl_mb(dn(1)))
           call Pack_c_pack(0,dbl_mb(dn(1)))

*          **** get Ecoul energy ****
           eh = coulomb_screened_e(dbl_mb(dn(1)))
           ph = 2.0d0*eh
           e_sic = e_sic + eh
           p_sic = p_sic + ph

          end do

        value = MA_pop_stack(tmp1(2))
        value = value.and.MA_pop_stack(xce(2))
        value = value.and.MA_pop_stack(xcp(2))
        value = value.and.MA_pop_stack(dn(2))
        if (.not. value) call errquit('popping stack memory',0)

      end if

      esic = e_sic
      psic = p_sic

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

#include "pspw_sic.fh"

      pspw_SIC_relaxed = relaxed
      return
      end
