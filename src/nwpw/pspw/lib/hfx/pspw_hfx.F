*
* $Id: pspw_hfx.F,v 1.1 2004-11-05 23:15:01 bylaska Exp $
*

*     *************************
*     *                       *
*     *     pspw_init_HFX     *
*     *                       *
*     *************************
      subroutine pspw_init_HFX(rtdb,ne)
      implicit none
      integer rtdb
      integer ne(2)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "pspw_hfx.fh"

*     **** local variables ****
      logical value
      integer ma_type
      integer n1,n2,n3,mapping

*     **** external functions ****
      integer  control_version,control_mapping
      external control_version,control_mapping

      norbs = 0
      ehfx = 0.0d0
      hfx_on = .false.


*     **** get the number of HFX orbitals ****
      if (rtdb_get(rtdb,'pspw:HFX',mt_log,1,value)) then
         hfx_on = .true.
         norbs = ne(1)+ne(2)
         if(.not.MA_alloc_get(mt_int,norbs,
     >       'orbital_list',orbital_list(2),orbital_list(1)))
     >        call errquit('pspw_init_HFX:out of heap memory',0,MA_ERR)
          
         do n1=1,norbs
            int_mb(orbital_list(1)+n1-1) = n1
         end do

      else if (rtdb_ma_get(rtdb, 'pspw:HFX_orbitals', ma_type,
     >                norbs, orbital_list(2))) then

         if (.not.MA_get_index(orbital_list(2),orbital_list(1))) 
     >      call errquit(
     >       'pspw_init_HFX: ma_get_index failed for actlist',911,
     &       MA_ERR)
         hfx_on = .true.
      end if


      if (hfx_on) then

         if (.not. rtdb_get(rtdb,
     >                      'pspw:HFX_screening_radius',
     >                      mt_dbl,1,rcut)) 
     >       rcut = 8.0d0

         if (.not. rtdb_get(rtdb,
     >                      'pspw:HFX_screening_power',
     >                      mt_dbl,1,pp)) 
     >       pp = 8.0d0

         if (.not. rtdb_get(rtdb,
     >                      'pspw:HFX_relax',
     >                      mt_log,1,relaxed)) 
     >       relaxed = .false.

         if (.not. rtdb_get(rtdb,
     >                      'pspw:HFX_solver_type',
     >                      mt_int,1,solver_type)) 
     >       solver_type = 1

         if (.not. rtdb_get(rtdb,
     >                      'pspw:HFX_parameter',
     >                       mt_dbl,1,HFX_parameter))
     >       HFX_parameter = 1.0d0

         if (.not. rtdb_get(rtdb,
     >                      'pspw:HFX_print_orbital_contribution',
     >                       mt_log,1,orb_contribution))
     >       orb_contribution = .false.
 

*        **** initialize coulomb_screened ****
         if (solver_type.eq.1) then
              call coulomb_screened_init(rcut,pp)

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
     >                'ehfx_orb',ehfx_orb(2),ehfx_orb(1))
           if (.not. value)   
     >       call errquit('pspw_init_HFX: out of heap memory',1, MA_ERR)
         end if

      else
         norbs = 0
      end if

   
      return
      end


*     *************************
*     *                       *
*     *     pspw_end_HFX      *
*     *                       *
*     *************************
      subroutine pspw_end_HFX()
      implicit none

#include "mafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

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
     >            dbl_mb(ehfx_orb(1)+i-1)
              end do
           end if
  487   format(//,'== Orbital Contributions to HFX ==')
  488   format(/1x,'orbital',15x,
     >         'HF_Exchange') 	
  489   format(1x,i7,2x,e18.6)
        end if


*       **** deallocate memory ****
        value =  MA_free_heap(orbital_list(2)) 
        if (orb_contribution) then
          value = value.and.MA_free_heap(ehfx_orb(2))
        end if
        if (.not. value)
     >  call errquit('pspw_end_HFX:error freeing heap memory',0, MA_ERR)


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
*     *     pspw_print_HFX    *
*     *                       *
*     *************************
      subroutine pspw_print_HFX(unit)
      implicit none
      integer unit

#include "mafdecls.fh"
#include "pspw_hfx.fh"

*     **** local variables ****
      integer i

      if (hfx_on) then
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
        else
          write(unit,1005)
        end if
        if (hfx_parameter.ne.1.0d0) write(unit,1008) hfx_parameter
        write(unit,*)

      end if

      return
 1001 FORMAT(6x,"- HFX relaxed")
 1002 FORMAT(6x,"- HFX unrelaxed")
 1003 FORMAT(6x,"- HFX orbitals:",10I5)
 1004 FORMAT(6x,"- HFX screened coulomb solver")
 1005 FORMAT(6x,"- HFX free-space coulomb solver")
 1006 FORMAT(6x,"- HFX screening radius(pspw:HFX_screening_radius):",
     >       E10.3)
 1007 FORMAT(6x,"- HFX screening power (pspw:HFX_screening_power) :",
     >       E10.3)
 1008 FORMAT(6x,"- HFX scaling parameter (pspw:HFX_parameter)     :",
     >       E10.3)
      end

*     ****************************
*     *                    	 *
*     *     pspw_potential_HFX   *
*     *                          *
*     ****************************
      subroutine pspw_potential_HFX(ispin,psi_r,Hpsi_r)
      implicit none
      integer    ispin
      real*8     psi_r(*)
      real*8     Hpsi_r(*)

#include "mafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n2ft3d,n1,n2,n3
      integer dn(2),vij(2),tmp1(2),tmp2(2),index1,index2
      real*8  scal1,scal2,dv,eh

*     **** external functions ****
      integer  control_gga
      real*8   lattice_omega,coulomb_screened_e
      external control_gga
      external lattice_omega,coulomb_screened_e


      ehfx = 0.0d0
      if ((norbs.ne.0).and.relaxed) then
        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        call D3dB_n2ft3d(1,n2ft3d)
        value = MA_push_get(mt_dbl,(n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(n2ft3d),'vij_hfx',vij(2),vij(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(n2ft3d),'tmp2_hfx',tmp2(2),tmp2(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2

        do i=1,norbs
        do j=1,i-1
           index1 = (int_mb(orbital_list(1)+i-1)-1)*n2ft3d + 1
           index2 = (int_mb(orbital_list(1)+j-1)-1)*n2ft3d + 1

*          **** generate dnij for Vij  ****
           call D3dB_rr_Mul(1,psi_r(index1),psi_r(index2),dbl_mb(dn(1)))
           call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))
   

*          ***** screened coulomb solver ****
           if (solver_type.eq.1) then

             call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
     >                              dbl_mb(dn(1)))
             call D3dB_rc_fft3f(1,dbl_mb(dn(1)))
             call Pack_c_pack(0,dbl_mb(dn(1)))

*            **** get Ecoul energy ****
             eh = coulomb_screened_e(dbl_mb(dn(1)))
             eh = 0.5d0*eh

*            **** generate Vcoul ****
             call coulomb_screened_v(dbl_mb(dn(1)),dbl_mb(vij(1)))
             call Pack_c_unpack(0,dbl_mb(vij(1)))
             call D3dB_cr_fft3b(1,dbl_mb(vij(1)))



*          ***** free-space coulomb solver ****
           else

              call coulomb2_v(dbl_mb(dn(1)),dbl_mb(vij(1)))
              call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(vij(1)),eh)

              eh = 0.5d0*eh*dv
           end if

*          **** apply the hfx_parameter ****
           eh = eh*hfx_parameter
           call D3dB_r_SMul(1,hfx_parameter,
     >                      dbl_mb(vij(1)),
     >                      dbl_mb(vij(1)))

           ehfx = ehfx - eh

*          **** generate (Vij)*psi_r ***
           call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                        psi_r(index1),
     >                        dbl_mb(tmp1(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))
           call D3dB_rr_Mul(1,dbl_mb(vij(1)),
     >                        psi_r(index2),
     >                        dbl_mb(tmp2(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(tmp2(1)))


*          **** add -(Vij)*psi_r to Hpsi_r ***
           call D3dB_rr_Sub(1,Hpsi_r(index2),
     >                        dbl_mb(tmp1(1)),
     >                        Hpsi_r(index2))

           call D3dB_rr_Sub(1,Hpsi_r(index1),
     >                        dbl_mb(tmp2(1)),
     >                        Hpsi_r(index1))


        end do
        end do

        value =           MA_pop_stack(tmp2(2))
        value = value.and.MA_pop_stack(tmp1(2))
        value = value.and.MA_pop_stack(vij(2))
        value = value.and.MA_pop_stack(dn(2))
        if (.not. value) 
     >    call errquit('pspw_potential_HFX:popping stack memory',0,
     &       MA_ERR)
      end if
      return
      end


*     *************************
*     *                       *
*     *     pspw_energy_HFX   *
*     *                       *
*     *************************
      subroutine pspw_energy_HFX(ispin,psi_r,ehfx_out)
      implicit none
      integer ispin
      real*8  psi_r(*)
      real*8 ehfx_out

#include "mafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n2ft3d,n1,n2,n3
      integer dn(2),vij(2),tmp1(2),index1,index2
      real*8  scal1,scal2,dv,eh

*     **** external functions ****
      integer  control_gga
      real*8   lattice_omega,coulomb_screened_e
      external control_gga
      external lattice_omega,coulomb_screened_e

      if ((norbs.ne.0).and.(.not.relaxed)) then
        ehfx = 0.0d0

        call D3dB_nx(1,n1)
        call D3dB_ny(1,n2)
        call D3dB_nz(1,n3)
        call D3dB_n2ft3d(1,n2ft3d)
        value = MA_push_get(mt_dbl,(2*n2ft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          MA_push_get(mt_dbl,(n2ft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2
        do i=1,norbs
        do j=1,i-1

           index1 = (int_mb(orbital_list(1)+i-1)-1)*n2ft3d + 1
           index2 = (int_mb(orbital_list(1)+j-1)-1)*n2ft3d + 1

*          **** generate dnij ****
           call D3dB_rr_MUL(1,psi_r(index1),psi_r(index2),dbl_mb(dn(1)))
           call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))


*          ***** screened coulomb solver ****
           if (solver_type.eq.1) then

*            **** generate dng ****
             call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
     >                                dbl_mb(dn(1)))
             call D3dB_rc_fft3f(1,dbl_mb(dn(1)))
             call Pack_c_pack(0,dbl_mb(dn(1)))

*            **** get Ecoul energy ****
             eh = coulomb_screened_e(dbl_mb(dn(1)))
             eh = 0.5d0*eh*hfx_parameter
       

*          ***** free-space coulomb solver ****
           else
              call coulomb2_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))

              call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(tmp1(1)),eh)
              eh = eh*dv
              eh = 0.5d0*eh*hfx_parameter
           end if

           ehfx = ehfx - eh

        end do
        end do

        value =           MA_pop_stack(tmp1(2))
        value = value.and.MA_pop_stack(dn(2))
        if (.not. value) 
     >     call errquit('pspw_energy_HFX:popping stack memory',0,
     &       MA_ERR)

      end if
      ehfx_out = ehfx


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

      pspw_HFX= hfx_on
      return
      end

*     *************************
*     *                       *
*     *   pspw_HFX_relaxed    *
*     *                       *
*     *************************
      logical function pspw_HFX_relaxed()
      implicit none

#include "mafdecls.fh"
#include "pspw_hfx.fh"

      pspw_hfx_relaxed = relaxed
      return
      end


