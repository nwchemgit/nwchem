c
c $Id: pspw_LJ.F,v 1.8 2002-02-28 00:34:21 bylaska Exp $
c

*     *************************
*     *                       *
*     *     pspw_LJ_init	  *
*     *                       *
*     *************************
      subroutine pspw_LJ_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"


***** LJparam common block ****
#include "pspw_LJ.fh"

*     **** local variables
      logical value
      integer taskid
      integer MASTER
      parameter(MASTER=0)

      integer i,j
      integer nkatm_mm,nkatm_ion
      real*8  sigma_epsilon(2)
      character*30 rtdbname

*     **** external functions ****
      character*2 pspw_qmmm_atom,ion_atom
      integer     pspw_qmmm_nkatm,ion_nkatm
      integer     pspw_qmmm_katm,ion_katm
      external    pspw_qmmm_atom,ion_atom
      external    pspw_qmmm_nkatm,ion_nkatm
      external    pspw_qmmm_katm,ion_katm
      

      nkatm_mm = pspw_qmmm_nkatm()
      nkatm_ion = ion_nkatm()

*     **** allocate LJ parameters ****
      value = rtdb_parallel(.true.)
      value = value.and.
     >        MA_alloc_get(mt_dbl,nkatm_mm,
     >                 'sigma_mm',sigma_mm(2),sigma_mm(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,nkatm_mm,
     >                 'epsilon_mm',epsilon_mm(2),epsilon_mm(1))
      value = MA_alloc_get(mt_dbl,nkatm_ion,
     >                 'sigma_ion',sigma_ion(2),sigma_ion(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,nkatm_mm,
     >                 'epsilon_ion',epsilon_ion(2),epsilon_ion(1))
      if (.not. value) call errquit('out of heap memory',0)


*     **** read in MM LJ parameters ***
      do i=1,nkatm_mm
         rtdbname = 'pspw_LJ_param_mm:'//pspw_qmmm_atom(i)
         j = index(rtdbname,' ')-1
         value = rtdb_get(rtdb,rtdbname(1:j),
     >                  mt_dbl,2,sigma_epsilon)
         if (.not. value) call dcopy(2,0.0d0,0,sigma_epsilon,1)

        dbl_mb(sigma_mm(1)+i-1)  =sigma_epsilon(1)/(0.529177d0)
        dbl_mb(epsilon_mm(1)+i-1)=sigma_epsilon(2)/(23.06d0*27.2116d0)
      end do
      
*     **** read in Ion LJ parameters ***
      do i=1,nkatm_ion
         rtdbname = 'pspw_LJ_param_ion:'//ion_atom(i)
         j = index(rtdbname,' ')-1
         value = rtdb_get(rtdb,rtdbname(1:j),
     >                  mt_dbl,2,sigma_epsilon)
         if (.not. value) call dcopy(2,0.0d0,0,sigma_epsilon,1)

        dbl_mb(sigma_ion(1)+i-1)  =sigma_epsilon(1)/(0.529177d0)
        dbl_mb(epsilon_ion(1)+i-1)=sigma_epsilon(2)/(23.06d0*27.2116d0)
      end do

      
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
          write(*,*)
          write(*,*) 'Ion LJ Parameters:'
          do i=1,nkatm_ion
            write(*,'(A3,4x,A12,E10.3,A9,A12,E10.3,A9)') 
     >         ion_atom(i),
     >         'sigma:',  
     >         dbl_mb(sigma_ion(1)+i-1)*(0.529177d0), 
     >         'Angstrom',
     >         'epsilon:',
     >         dbl_mb(epsilon_ion(1)+ i-1)*(23.06d0*27.2116d0),
     >         'kcal/mol'
          end do
         write(*,*)
         write(*,*) 'MM LJ Parameters:'
         do i=1,nkatm_mm
            write(*,'(A3,4x,A12,E10.3,A9,A12,E10.3,A9)') 
     >         pspw_qmmm_atom(i),
     >         'sigma:',  
     >         dbl_mb(sigma_mm(1)+i-1)*(0.529177d0), 
     >         'Angstrom',
     >         'epsilon:',
     >         dbl_mb(epsilon_mm(1)+i-1)*(23.06d0*27.2116d0),
     >         'kcal/mol'
          end do
          write(*,*)

      end if

      return
      end

*     *************************
*     *                       *
*     *     pspw_LJ_end		  *
*     *                       *
*     *************************
      subroutine pspw_LJ_end()
      implicit none

#include "mafdecls.fh"

***** LJparam common block ****
#include "pspw_LJ.fh"

      logical value

      value = MA_free_heap(sigma_mm(2))
      value = value.and.MA_free_heap(epsilon_mm(2))
      value = value.and.MA_free_heap(sigma_ion(2))
      value = value.and.MA_free_heap(epsilon_ion(2))
      if (.not.value) call errquit('error MA_free_heap',0)

      return
      end


*     ***************************
*     *							*
*     *    pspw_LJ_param_Ion	*
*     *							*
*     ***************************
      subroutine pspw_LJ_param_Ion(i,s,e)
      implicit none
      integer i
      real*8  s,e

#include "mafdecls.fh"

***** LJparam common block ****
#include "pspw_LJ.fh"

*     **** local variables ****
      integer ka

*     **** external functions ****
      integer  ion_katm
      external ion_katm

      ka = ion_katm(i)
      s  = dbl_mb(sigma_ion(1)  + ka-1)
      e  = dbl_mb(epsilon_ion(1)+ ka-1)
      return
      end

*     ***************************
*     *							*
*     *    pspw_LJ_param_MM 	*
*     *							*
*     ***************************
      subroutine pspw_LJ_param_MM(i,s,e)
      implicit none
      integer i
      real*8  s,e

#include "mafdecls.fh"

***** LJparam common block ****
#include "pspw_LJ.fh"

*     **** local variables ****
      integer ka

*     **** external functions ****
      integer  pspw_qmmm_katm
      external pspw_qmmm_katm

      ka = pspw_qmmm_katm(i)
      s  = dbl_mb(sigma_mm(1)  + ka-1)
      e  = dbl_mb(epsilon_mm(1)+ ka-1)
      return
      end


      
*     ***************************
*     *							*
*     *   pspw_LJ_Energy_MM  	*
*     *							*
*     ***************************
      real*8 function pspw_LJ_Energy_MM()
      implicit none

*     **** local variables ****
      real*8  E,si,ei,sj,ej
      real*8 ri(3),rj(3)
      integer i,j,nion

*     **** external functions ****
      logical  pspw_qmmm_different_frag
      integer  pspw_qmmm_nion
      real*8   pspw_qmmm_rion,LJ_Repulsion
      external pspw_qmmm_different_frag
      external pspw_qmmm_nion
      external pspw_qmmm_rion,LJ_Repulsion

      nion = pspw_qmmm_nion()
     
      E = 0.0d0
      do i=1,nion
         call pspw_LJ_param_MM(i,si,ei)
         if (ei.gt.1.0d-10) then
         ri(1) = pspw_qmmm_rion(1,i)
         ri(2) = pspw_qmmm_rion(2,i)
         ri(3) = pspw_qmmm_rion(3,i)
         do j=i+1,nion
           if (pspw_qmmm_different_frag(i,j)) then
             call pspw_LJ_param_MM(j,sj,ej)
             if (ej.gt.1.0d-10) then
             rj(1) = pspw_qmmm_rion(1,j)
             rj(2) = pspw_qmmm_rion(2,j)
             rj(3) = pspw_qmmm_rion(3,j)
             E = E + LJ_Repulsion(ri,si,ei,rj,sj,ej)
             end if
           end if
         end do
         end if
      end do


      pspw_LJ_Energy_MM = E
      return
      end



*     ***************************
*     *							*
*     *   pspw_LJ_Energy_IonMM 	*
*     *							*
*     ***************************
      real*8 function pspw_LJ_Energy_IonMM()
      implicit none

*     **** local variables ****
      real*8  E,si,ei,sj,ej
      real*8 ri(3),rj(3)
      integer i,j,nion_ion,nion_mm

*     **** external functions ****
      integer  pspw_qmmm_nion,ion_nion
      real*8   pspw_qmmm_rion,ion_rion,LJ_Repulsion
      external pspw_qmmm_nion,ion_nion
      external pspw_qmmm_rion,ion_rion,LJ_Repulsion

      nion_mm  = pspw_qmmm_nion()
      nion_ion = ion_nion()
     
      E = 0.0d0
      do i=1,nion_ion
         call pspw_LJ_param_Ion(i,si,ei)
         if (ei.gt.1.0d-10) then
         ri(1) = ion_rion(1,i)
         ri(2) = ion_rion(2,i)
         ri(3) = ion_rion(3,i)
         do j=1,nion_mm
            call pspw_LJ_param_MM(j,sj,ej)
            if (ej.gt.1.0d-10) then
            rj(1) = pspw_qmmm_rion(1,j)
            rj(2) = pspw_qmmm_rion(2,j)
            rj(3) = pspw_qmmm_rion(3,j)
            E = E + LJ_Repulsion(ri,si,ei,rj,sj,ej)
            end if
         end do
         end if
      end do


      pspw_LJ_Energy_IonMM = E
      return
      end


*     ***************************
*     *							*
*     *       pspw_LJ_Fmm	 	*
*     *							*
*     ***************************

*     This routine calculates the force on the MM ions from the
*   Lenard-Jones interaction with the MM atoms and the QM ions and 
*   then adds it to the MM ion force, f_mm.
*
*     Exit - f_mm: force on the MM ions
*
*     Uses - ion_nion,ion_rion
*            pspw_qmmm_nion,pspw_qmmm_rion
*            pspw_qmmm_different_frag
*
*     Author - Eric Bylaska
*
      subroutine pspw_LJ_Fmm(f_mm)
      implicit none
      real*8 f_mm(3,*)

*     **** local variables ****
      real*8  si,ei,sj,ej
      real*8 ri(3),rj(3),fion(3)
      integer i,j

*     **** external functions ****
      logical  pspw_qmmm_different_frag
      integer  pspw_qmmm_nion,ion_nion
      real*8   pspw_qmmm_rion,ion_rion
      external pspw_qmmm_different_frag
      external pspw_qmmm_nion,ion_nion
      external pspw_qmmm_rion,ion_rion


*     **** force from the MM atoms ****
      do i=1,pspw_qmmm_nion()
         call pspw_LJ_param_MM(i,si,ei)
         if (ei.gt.1.0d-10) then
         ri(1) = pspw_qmmm_rion(1,i)
         ri(2) = pspw_qmmm_rion(2,i)
         ri(3) = pspw_qmmm_rion(3,i)
         do j=i+1,pspw_qmmm_nion()
           if (pspw_qmmm_different_frag(i,j)) then
             call pspw_LJ_param_MM(j,sj,ej)
             if (ej.gt.1.0d-10) then
             rj(1) = pspw_qmmm_rion(1,j)
             rj(2) = pspw_qmmm_rion(2,j)
             rj(3) = pspw_qmmm_rion(3,j)
             call LJ_Ion_Repulsion_Force(ri,si,ei,f_mm(1,i),
     >                                   rj,sj,ej,f_mm(1,j))
             end if  
           end if  
         end do
         end if
      end do


*     **** force from the QM atoms ****
      call dcopy(3,0.0d0,0,fion,1)
      do i=1,ion_nion()
         call pspw_LJ_param_Ion(i,si,ei)
         if (ei.gt.1.0d-10) then
         ri(1) = ion_rion(1,i)
         ri(2) = ion_rion(2,i)
         ri(3) = ion_rion(3,i)
         do j=1,pspw_qmmm_nion()
            call pspw_LJ_param_MM(j,sj,ej)
            if (ej.gt.1.0d-10) then
            rj(1) = pspw_qmmm_rion(1,j)
            rj(2) = pspw_qmmm_rion(2,j)
            rj(3) = pspw_qmmm_rion(3,j)
            call LJ_Ion_Repulsion_Force(ri,si,ei,fion,
     >                                  rj,sj,ej,f_mm(1,j))
            end if
         end do
         end if
      end do

      return
      end


*     ***************************
*     *							*
*     *   pspw_LJ_Fion 			*
*     *							*
*     ***************************

*     This routine calculates the force on the QM ions from the 
*   Lenard-Jones interaction with the MM atoms and adds it
*   to the QM ion force, fion.
*
*     Exit - fion: force on the QM ions
*
*     Uses - ion_nion,ion_rion
*            pspw_qmmm_nion,pspw_qmmm_rion
*
*     Author - Eric Bylaska

      subroutine pspw_LJ_Fion(fion)
      implicit none
      real*8 fion(3,*)

*     **** local variables ****
      real*8  si,ei,sj,ej
      real*8 ri(3),rj(3),f_mm(3)
      integer i,j,nion_ion,nion_mm

*     **** external functions ****
      integer  pspw_qmmm_nion,ion_nion
      real*8   pspw_qmmm_rion,ion_rion
      external pspw_qmmm_nion,ion_nion
      external pspw_qmmm_rion,ion_rion

      nion_mm  = pspw_qmmm_nion()
      nion_ion = ion_nion()
     
      call dcopy(3,0.0d0,0,f_mm,1)
      do i=1,nion_ion
         call pspw_LJ_param_Ion(i,si,ei)
         if (ei.gt.1.0d-10) then
         ri(1) = ion_rion(1,i)
         ri(2) = ion_rion(2,i)
         ri(3) = ion_rion(3,i)
         do j=1,nion_mm
            call pspw_LJ_param_MM(j,sj,ej)
            if (ej.gt.1.0d-10) then
            rj(1) = pspw_qmmm_rion(1,j)
            rj(2) = pspw_qmmm_rion(2,j)
            rj(3) = pspw_qmmm_rion(3,j)
            call LJ_Ion_Repulsion_Force(ri,si,ei,fion(1,i),
     >                                  rj,sj,ej,f_mm)
            end if
         end do
         end if
      end do

      return
      end


*     **********************************
*     *	                               *
*     *          LJ_Repulsion		   *
*     *                                *
*     **********************************

*   This routine calculates the LJ Repulsion between two centers 
*
      real*8 function LJ_Repulsion(r1,s1,e1,r2,s2,e2)
      implicit none
      real*8 r1(3),s1,e1
      real*8 r2(3),s2,e2

*     **** local variables ****
      real*8 xo,yo,zo,d2
      real*8 ee,ss
     
      ee = dsqrt(e1*e2)
      ss = (s1*s2)

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)
  
      d2 = (xo*xo + yo*yo + zo*zo)
        

      LJ_Repulsion = 4.0d0*ee * ((ss/d2)**6 - (ss/d2)**3)    
      return 
      end

*     **********************************
*     *	                               *
*     *      LJ_Repulsion_Force   	   *
*     *                                *
*     **********************************

*   This routine calculates the LJ Repulsion force between two centers
*
      subroutine LJ_Ion_Repulsion_Force(r1,s1,e1,f1,r2,s2,e2,f2)
      implicit none
      real*8 r1(3),s1,e1
      real*8 f1(3)
      real*8 r2(3),s2,e2
      real*8 f2(3)

*     **** local variables ****
      real*8 xo,yo,zo,d2
      real*8 der
      real*8 ee,ss,fx,fy,fz
     
      ee = dsqrt(e1*e2)
      ss = (s1*s2)

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)
  
      d2 = (xo*xo + yo*yo + zo*zo)

      der = (4.0d0*ee/d2) * ( -12.0d0*(ss/d2)**6 
     >                       +  6.0d0*(ss/d2)**3 ) 
      fx = -(xo)*der
      fy = -(yo)*der
      fz = -(zo)*der  

      f2(1) = f2(1) + fx
      f2(2) = f2(2) + fy
      f2(3) = f2(3) + fz

      f1(1) = f1(1) - fx
      f1(2) = f1(2) - fy
      f1(3) = f1(3) - fz
      return 
      end



