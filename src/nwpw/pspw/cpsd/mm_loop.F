
*
* $Id$
*

*  ************************************************************
*  *                MPI mmsd routine                          *
*  *                                                          *
*  *  This is a developing mmsdv1 parallel code wrtten in     *
*  *  Fortran and MPI.                                        *
*  *                                                          *
*  ************************************************************
  
      subroutine mm_loop(it_in,E,deltae,deltar)
      implicit none
      integer    it_in
      real*8     E(*)
      real*8     deltae,deltar

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer i,n,n1(2),n2(2),it,ms,ierr
      real*8  sum,Eold,eion
      real*8  dt
      real*8  e_lj, e_q, e_spring

*     **** MA local variables ****
      logical value
      integer natmx,fion(2),ftest(2)

*     ***** external functions ****
      integer  ion_nion
      real*8   control_time_step
      real*8   ewald_e
      external ion_nion
      external control_time_step
      external ewald_e
      integer  control_version
      external control_version
      real*8   ion_ion_e
      external ion_ion_e

*     ***** QM/MM external functions ****
      logical  pspw_charge_found,pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_charge_found,pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E


      call nwpw_timing_start(12)
*     **** allocate MA local variables ****
      natmx = ion_nion()
      value = BA_push_get(mt_dbl,(3*natmx),'fion',fion(2),fion(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(3*natmx),'ftest',ftest(2),ftest(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      call nwpw_timing_end(12)
     
      dt = control_time_step()

*     ******************************************
*     ****                                  ****
*     ****   Start of steepest descent loop ****
*     ****                                  ****
*     ******************************************
      do it=1,it_in

        call ion_shift()
        call dcopy(3*natmx,0.0d0,0,dbl_mb(fion(1)),1)

        if (control_version().eq.3) then
           call ewald_phafac()
           call ewald_f(dbl_mb(fion(1)))
        end if
        if (control_version().eq.4) call ion_ion_f(dbl_mb(fion(1)))

        call pspw_qmmm_fion(dbl_mb(fion(1)))

*       **** remove ion forces using ion_FixIon ****
        call ion_FixIon(dbl_mb(fion(1)))
        call ion_optimize_step(dbl_mb(fion(1)))

      end do

*     *************************************
*     ***** total energy calculation ******
*     *************************************
      call nwpw_timing_start(10)

*     **** get ewald energy ****
      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()
    
*     **** get free-space ion-ion energy ****
      if (control_version().eq.4) eion = ion_ion_e()

*     **** total energy ****
      Eold=E(1)
      E(1) = eion 
      E(2) = 0.0d0
      E(3) = 0.0d0
      E(4) = 0.0d0
      E(5) = eion
      E(6) = 0.0d0
      E(7) = 0.0d0
      E(8) = 0.0d0
      E(9) = 0.0d0
      E(10) = 0.0d0
      e_lj  = pspw_qmmm_LJ_E()
      e_q   = pspw_qmmm_Q_E()
      e_spring = pspw_qmmm_spring_E()
      E(1)  = E(1) + e_lj + e_q + e_spring
      E(11) = e_lj
      E(12) = e_q
      E(13) = e_spring

*     **** set convergence variables ****
      deltae = (E(1)-Eold)/(dt*dble(it_in))

*     *** deltar ***
      deltar = 0.0d0
      do i=1,ion_nion()
         sum = dsqrt( dbl_mb(fion(1)+(i-1)*3  )**2 
     >              + dbl_mb(fion(1)+(i-1)*3+1)**2 
     >              + dbl_mb(fion(1)+(i-1)*3+2)**2)
        if (sum.gt.deltar) deltar = sum
      end do
       
      call nwpw_timing_end(10)

*     **** dealocate MA local variables ****
      call nwpw_timing_start(12)
      value = BA_pop_stack(ftest(2))
      value = BA_pop_stack(fion(2))

      call nwpw_timing_end(12)
      return
      end
 
