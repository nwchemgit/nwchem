*     *************************
*     *                       *
*     *     pspw_Q_init	  *
*     *                       *
*     *************************
      subroutine pspw_Q_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"


***** LJparam common block ****
#include "pspw_Q.fh"

*     **** local variables
      logical value
      integer taskid
      integer MASTER
      parameter(MASTER=0)

      integer i,j
      integer nkatm_mm,nkatm_ion
      real*8  q
      character*30 rtdbname

*     **** external functions ****
      character*2 pspw_qmmm_atom,ion_atom
      integer     pspw_qmmm_nkatm,ion_nkatm
      integer     pspw_qmmm_katm,ion_katm
      real*8      psp_zv
      external    pspw_qmmm_atom,ion_atom
      external    pspw_qmmm_nkatm,ion_nkatm
      external    pspw_qmmm_katm,ion_katm
      external    psp_zv
      

      nkatm_mm = pspw_qmmm_nkatm()
      nkatm_ion = ion_nkatm()

*     **** allocate Electrostatic parameters ****
*     **** default is zero                   ****
      value = rtdb_parallel(.true.)
      value = value.and.
     >        MA_alloc_get(mt_dbl,nkatm_mm,
     >                 'q_mm',q_mm(2),q_mm(1))
      value = MA_alloc_get(mt_dbl,nkatm_ion,
     >                 'q_ion',q_ion(2),q_ion(1))
      if (.not. value) call errquit('out of heap memory',0)


*     **** read in MM Electrostatic parameters ***
      do i=1,nkatm_mm
         rtdbname = 'pspw_Q_param_mm:'//pspw_qmmm_atom(i)
         j = index(rtdbname,' ')-1
         value = rtdb_get(rtdb,rtdbname(1:j),
     >                  mt_dbl,1,q)
         if (.not. value) q=0.0d0

         dbl_mb(q_mm(1)+i-1) = q
      end do
      
*     **** read in Ion Electrostatic parameters      ****
*     **** default is core charge of pseudopotential ****
      do i=1,nkatm_ion
         rtdbname = 'pspw_Q_param_ion:'//ion_atom(i)
         j = index(rtdbname,' ')-1
         value = rtdb_get(rtdb,rtdbname(1:j),
     >                  mt_dbl,1,q)
         if (.not. value) q=psp_zv(i)

         dbl_mb(q_ion(1)+i-1) = q
      end do
      
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
          write(*,*)
          write(*,*) 'Ion Electrostatic Parameters:'
          do i=1,nkatm_ion
            write(*,'(A3,4x,A12,E10.3)') 
     >                  ion_atom(i),
     >                  'Q:',  dbl_mb(q_ion(1)  + i-1) 
          end do
         write(*,*)
         write(*,*) 'MM Electrostatic Parameters:'
         do i=1,nkatm_mm
            write(*,'(A3,4x,A12,E10.3)') 
     >                  pspw_qmmm_atom(i),
     >                  'Q:',  dbl_mb(q_mm(1)  + i-1) 
          end do
          write(*,*)

      end if

      return
      end

*     *************************
*     *                       *
*     *     pspw_Q_end		  *
*     *                       *
*     *************************
      subroutine pspw_Q_end()
      implicit none

#include "mafdecls.fh"

***** LJparam common block ****
#include "pspw_Q.fh"

      logical value

      value = MA_free_heap(q_mm(2))
      value = value.and.MA_free_heap(q_ion(2))
      if (.not.value) call errquit('error MA_free_heap',0)

      return
      end


*     ***************************
*     *							*
*     *    pspw_Q_param_Ion		*
*     *							*
*     ***************************
      subroutine pspw_Q_param_Ion(i,q)
      implicit none
      integer i
      real*8  q

#include "mafdecls.fh"

***** LJparam common block ****
#include "pspw_Q.fh"

*     **** local variables ****
      integer ka

*     **** external functions ****
      integer  ion_katm
      external ion_katm

      ka = ion_katm(i)
      q  = dbl_mb(q_ion(1)  + ka-1)
      return
      end

*     ***************************
*     *							*
*     *    pspw_Q_param_MM 	*
*     *							*
*     ***************************
      subroutine pspw_Q_param_MM(i,q)
      implicit none
      integer i
      real*8  q

#include "mafdecls.fh"

***** LJparam common block ****
#include "pspw_Q.fh"

*     **** local variables ****
      integer ka

*     **** external functions ****
      integer  pspw_qmmm_katm
      external pspw_qmmm_katm

      ka = pspw_qmmm_katm(i)
      q  = dbl_mb(q_mm(1)  + ka-1)
      return
      end

*     ***************************
*     *							*
*     *    pspw_Q_zv_MM 	*
*     *							*
*     ***************************
      real*8 function pspw_Q_zv_MM(ia)
      implicit none
      integer ia

#include "mafdecls.fh"

***** LJparam common block ****
#include "pspw_Q.fh"

      pspw_Q_zv_MM =  dbl_mb(q_mm(1) + ia-1)
      return
      end


      
*     ***************************
*     *							*
*     *   pspw_Q_Energy_MM  	*
*     *							*
*     ***************************
      real*8 function pspw_Q_Energy_MM()
      implicit none

*     **** local variables ****
      real*8  E,qi,qj
      real*8 ri(3),rj(3)
      integer i,j,nion

*     **** external functions ****
      logical  pspw_qmmm_different_frag
      integer  pspw_qmmm_nion
      real*8   pspw_qmmm_rion,Q_Electrostatic
      external pspw_qmmm_different_frag
      external pspw_qmmm_nion
      external pspw_qmmm_rion,Q_Electrostatic

      nion = pspw_qmmm_nion()
     
      E = 0.0d0
      do i=1,nion
         call pspw_Q_param_MM(i,qi)
         ri(1) = pspw_qmmm_rion(1,i)
         ri(2) = pspw_qmmm_rion(2,i)
         ri(3) = pspw_qmmm_rion(3,i)
         do j=i+1,nion
           if (pspw_qmmm_different_frag(i,j)) then
             call pspw_Q_param_MM(j,qj)
             rj(1) = pspw_qmmm_rion(1,j)
             rj(2) = pspw_qmmm_rion(2,j)
             rj(3) = pspw_qmmm_rion(3,j)
             E = E + Q_Electrostatic(ri,qi,rj,qj)
           end if
         end do
      end do


      pspw_Q_Energy_MM = E
      return
      end



*     ***************************
*     *							*
*     *   pspw_Q_Energy_IonMM 	*
*     *							*
*     ***************************
      real*8 function pspw_Q_Energy_IonMM()
      implicit none

*     **** local variables ****
      real*8  E,qi,qj
      real*8 ri(3),rj(3)
      integer i,j,nion_ion,nion_mm

*     **** external functions ****
      integer  pspw_qmmm_nion,ion_nion
      real*8   pspw_qmmm_rion,ion_rion,Q_Electrostatic
      external pspw_qmmm_nion,ion_nion
      external pspw_qmmm_rion,ion_rion,Q_Electrostatic

      nion_mm  = pspw_qmmm_nion()
      nion_ion = ion_nion()
     
      E = 0.0d0
      do i=1,nion_ion
         call pspw_Q_param_Ion(i,qi)
         ri(1) = ion_rion(1,i)
         ri(2) = ion_rion(2,i)
         ri(3) = ion_rion(3,i)
         do j=1,nion_mm
            call pspw_Q_param_MM(j,qj)
            rj(1) = pspw_qmmm_rion(1,j)
            rj(2) = pspw_qmmm_rion(2,j)
            rj(3) = pspw_qmmm_rion(3,j)
            E = E + Q_Electrostatic(ri,qi,rj,qj)
         end do
      end do


      pspw_Q_Energy_IonMM = E
      return
      end


*     ***************************
*     *							*
*     *   pspw_Q_Force_MM	 	*
*     *							*
*     ***************************
      subroutine pspw_Q_Force_MM(f_mm)
      implicit none
      real*8 f_mm(3,*)

*     **** local variables ****
      real*8  qi,qj
      real*8 ri(3),rj(3)
      integer i,j,nion

*     **** external functions ****
      logical  pspw_qmmm_different_frag
      integer  pspw_qmmm_nion
      real*8   pspw_qmmm_rion
      external pspw_qmmm_different_frag
      external pspw_qmmm_nion
      external pspw_qmmm_rion

      nion = pspw_qmmm_nion()

      do i=1,nion
         call pspw_Q_param_MM(i,qi)
         ri(1) = pspw_qmmm_rion(1,i)
         ri(2) = pspw_qmmm_rion(2,i)
         ri(3) = pspw_qmmm_rion(3,i)
         do j=i+1,nion
           if (pspw_qmmm_different_frag(i,j)) then
             call pspw_Q_param_MM(j,qj)
             rj(1) = pspw_qmmm_rion(1,j)
             rj(2) = pspw_qmmm_rion(2,j)
             rj(3) = pspw_qmmm_rion(3,j)
             call Q_Electrostatic_Force(ri,qi,f_mm(1,i),
     >                                  rj,qj,f_mm(1,j))
           end if  
         end do
      end do

      return
      end


*     ***************************
*     *							*
*     *   pspw_Q_Force_IonMM 	*
*     *							*
*     ***************************
      subroutine pspw_Q_Force_IonMM(fion,f_mm)
      implicit none
      real*8 fion(3,*)
      real*8 f_mm(3,*)

*     **** local variables ****
      real*8  qi,qj
      real*8 ri(3),rj(3)
      integer i,j,nion_ion,nion_mm

*     **** external functions ****
      integer  pspw_qmmm_nion,ion_nion
      real*8   pspw_qmmm_rion,ion_rion
      external pspw_qmmm_nion,ion_nion
      external pspw_qmmm_rion,ion_rion

      nion_mm  = pspw_qmmm_nion()
      nion_ion = ion_nion()
     
      do i=1,nion_ion
         call pspw_Q_param_Ion(i,qi)
         ri(1) = ion_rion(1,i)
         ri(2) = ion_rion(2,i)
         ri(3) = ion_rion(3,i)
         do j=1,nion_mm
            call pspw_Q_param_MM(j,qj)
            rj(1) = pspw_qmmm_rion(1,j)
            rj(2) = pspw_qmmm_rion(2,j)
            rj(3) = pspw_qmmm_rion(3,j)
            call Q_Electrostatic_Force(ri,qi,fion(1,i),
     >                                 rj,qj,f_mm(1,j))
         end do
      end do

      return
      end


*     **********************************
*     *	                               *
*     *          Q_Electrostatic	   *
*     *                                *
*     **********************************

*   This routine calculates the Electrostatic energy between two centers 
*
      real*8 function Q_Electrostatic(r1,q1,r2,q2)
      implicit none
      real*8 r1(3),q1
      real*8 r2(3),q2

*     **** local variables ****
      real*8 xo,yo,zo,d
     

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)
  
      d = dsqrt(xo*xo + yo*yo + zo*zo)

      Q_Electrostatic = q1*q2/d
      return 
      end

*     **********************************
*     *	                               *
*     *      Q_Electrostatic_Force 	   *
*     *                                *
*     **********************************

*   This routine calculates the Electrostaticn force between two centers
*
      subroutine Q_Electrostatic_Force(r1,q1,f1,r2,q2,f2)
      implicit none
      real*8 r1(3),q1
      real*8 f1(3)
      real*8 r2(3),q2
      real*8 f2(3)

*     **** local variables ****
      real*8 xo,yo,zo,d0,d2
      real*8 der
      real*8 fx,fy,fz
     

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)
  
      d2 = (xo*xo + yo*yo + zo*zo)
      d0 = dsqrt(d2)

      der = -q1*q2/d2
      fx = -(xo/d0)*der
      fy = -(yo/d0)*der
      fz = -(zo/d0)*der  

      f2(1) = f2(1) + fx
      f2(2) = f2(2) + fy
      f2(3) = f2(3) + fz

      f1(1) = f1(1) - fx
      f1(2) = f1(2) - fy
      f1(3) = f1(3) - fz
      return 
      end



